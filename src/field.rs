use std::ops::{Add, Mul, Sub};

pub trait Field: Copy + Eq + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> {
    const ZERO: Self;
    const ONE: Self;

    fn inverse(self) -> Option<Self>;

    fn pow(self, mut n: usize) -> Self {
        let mut result = Self::ONE;
        let mut base = self;
        while n > 0 {
            if n & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            n >>= 1;
        }
        result
    }
}

// GF(2)[x] / (x^8 + x^4 + x^3 + x^2 + 1)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Gf256(pub u8);

impl Field for Gf256 {
    const ZERO: Self = Self(0);
    const ONE: Self = Self(1);

    fn inverse(self) -> Option<Self> {
        (self != Self::ZERO).then(|| self.pow(254))
    }
}

impl Add for Gf256 {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, rhs: Self) -> Self {
        Self(self.0 ^ rhs.0)
    }
}

impl Sub for Gf256 {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 ^ rhs.0)
    }
}

impl Mul for Gf256 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut a = self.0;
        let mut b = rhs.0;
        let mut result = 0;
        while b != 0 {
            if b & 1 != 0 {
                result ^= a;
            }
            let carry = a & 0x80 != 0;
            a <<= 1;
            if carry {
                a ^= 0x1d;
            }
            b >>= 1;
        }
        Self(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverses() {
        assert_eq!(Gf256::ZERO.inverse(), None);
        for n in 1..=255 {
            let a = Gf256(n);
            assert_eq!(a * a.inverse().unwrap(), Gf256::ONE);
        }
        assert_eq!(Gf256(0x80) * Gf256(2), Gf256(0x1d));
    }
}
