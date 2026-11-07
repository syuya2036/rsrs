use crate::field::{Field, Gf256};
use crate::polynomial::Polynomial;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    InvalidParity,
    InvalidLength,
    Uncorrectable,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoded_polynomial_has_generator_roots() {
        let rs = ReedSolomon::new(8).unwrap();
        let encoded = rs.encode(b"hello").unwrap();
        assert_eq!(&encoded[..5], b"hello");
        let polynomial = Polynomial::new(encoded.into_iter().map(Gf256).collect());
        for i in 0..8 {
            assert_eq!(polynomial.evaluate(Gf256(2).pow(i)), Gf256::ZERO);
        }
        assert_eq!(rs.encode(&[0; 12]).unwrap(), vec![0; 20]);
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidParity => "parity must be between 1 and 254",
            Self::InvalidLength => {
                "expected 1..=255-parity data bytes or parity+1..=255 codeword bytes"
            }
            Self::Uncorrectable => "could not correct the codeword",
        })
    }
}

impl std::error::Error for Error {}

pub struct ReedSolomon {
    parity: usize,
    generator: Polynomial<Gf256>,
}

impl ReedSolomon {
    pub fn new(parity: usize) -> Result<Self, Error> {
        if !(1..=254).contains(&parity) {
            return Err(Error::InvalidParity);
        }
        let mut generator = Polynomial::new(vec![Gf256::ONE]);
        for i in 0..parity {
            let factor = Polynomial::new(vec![Gf256::ONE, Gf256(2).pow(i)]);
            generator = generator.multiply(&factor);
        }
        Ok(Self { parity, generator })
    }

    pub fn encode(&self, data: &[u8]) -> Result<Vec<u8>, Error> {
        if data.is_empty() || data.len() > 255 - self.parity {
            return Err(Error::InvalidLength);
        }
        let mut coefficients: Vec<_> = data.iter().copied().map(Gf256).collect();
        coefficients.resize(data.len() + self.parity, Gf256::ZERO);
        let remainder = Polynomial::new(coefficients)
            .remainder(&self.generator)
            .unwrap();
        let mut result = data.to_vec();
        result.resize(data.len() + self.parity, 0);
        let start = result.len() - remainder.coefficients().len();
        for (i, c) in remainder.coefficients().iter().enumerate() {
            result[start + i] = c.0;
        }
        Ok(result)
    }
}
