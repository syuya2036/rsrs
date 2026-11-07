use crate::field::{Field, Gf256};
use crate::polynomial::Polynomial;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    InvalidParity,
    InvalidLength,
    Uncorrectable,
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

    pub fn decode(&self, codeword: &[u8]) -> Result<Vec<u8>, Error> {
        if codeword.len() <= self.parity || codeword.len() > 255 {
            return Err(Error::InvalidLength);
        }
        let mut word: Vec<_> = codeword.iter().copied().map(Gf256).collect();
        let syndromes = self.syndromes(&word);
        if syndromes.iter().any(|&s| s != Gf256::ZERO) {
            let locator = error_locator(&syndromes);
            let errors = locator.degree().ok_or(Error::Uncorrectable)?;
            if errors == 0 || errors * 2 > self.parity {
                return Err(Error::Uncorrectable);
            }

            // 位置 p にある誤りは x = alpha^(n-1-p) と対応する。
            let mut positions = Vec::new();
            let mut roots = Vec::new();
            for p in 0..word.len() {
                let x = Gf256(2).pow(word.len() - 1 - p);
                if locator.evaluate(x.inverse().unwrap()) == Gf256::ZERO {
                    positions.push(p);
                    roots.push(x);
                }
            }
            if positions.len() != errors {
                return Err(Error::Uncorrectable);
            }
            let magnitudes = error_magnitudes(&roots, &syndromes)?;
            for (&p, magnitude) in positions.iter().zip(magnitudes) {
                word[p] = word[p] - magnitude;
            }
            if self.syndromes(&word).iter().any(|&s| s != Gf256::ZERO) {
                return Err(Error::Uncorrectable);
            }
        }
        word.truncate(word.len() - self.parity);
        Ok(word.into_iter().map(|c| c.0).collect())
    }

    fn syndromes(&self, word: &[Gf256]) -> Vec<Gf256> {
        let polynomial = Polynomial::new(word.to_vec());
        (0..self.parity)
            .map(|i| polynomial.evaluate(Gf256(2).pow(i)))
            .collect()
    }
}

// Berlekamp–Massey。この計算中だけ係数を低い次数から持つ。
fn error_locator(syndromes: &[Gf256]) -> Polynomial<Gf256> {
    let mut current = vec![Gf256::ONE];
    let mut previous = vec![Gf256::ONE];
    let mut degree = 0;
    let mut shift = 1;
    let mut last = Gf256::ONE;

    for (n, &syndrome) in syndromes.iter().enumerate() {
        let mut discrepancy = syndrome;
        for i in 1..=degree {
            discrepancy = discrepancy + current[i] * syndromes[n - i];
        }
        if discrepancy == Gf256::ZERO {
            shift += 1;
            continue;
        }

        let saved = current.clone();
        let scale = discrepancy * last.inverse().unwrap();
        current.resize(current.len().max(previous.len() + shift), Gf256::ZERO);
        for (i, &coefficient) in previous.iter().enumerate() {
            current[i + shift] = current[i + shift] - scale * coefficient;
        }
        if 2 * degree <= n {
            degree = n + 1 - degree;
            previous = saved;
            last = discrepancy;
            shift = 1;
        } else {
            shift += 1;
        }
    }
    current.reverse();
    Polynomial::new(current)
}

fn error_magnitudes(roots: &[Gf256], syndromes: &[Gf256]) -> Result<Vec<Gf256>, Error> {
    let n = roots.len();
    // S_i = sum(e_j * X_j^i) を掃き出し法で解く。
    let mut matrix: Vec<Vec<_>> = (0..n)
        .map(|i| {
            let mut row: Vec<_> = roots.iter().map(|x| x.pow(i)).collect();
            row.push(syndromes[i]);
            row
        })
        .collect();

    for col in 0..n {
        let pivot = (col..n)
            .find(|&row| matrix[row][col] != Gf256::ZERO)
            .ok_or(Error::Uncorrectable)?;
        matrix.swap(col, pivot);
        let inverse = matrix[col][col].inverse().ok_or(Error::Uncorrectable)?;
        for value in &mut matrix[col][col..] {
            *value = *value * inverse;
        }
        let pivot_row = matrix[col].clone();
        for (row, values) in matrix.iter_mut().enumerate() {
            if row == col {
                continue;
            }
            let factor = values[col];
            for j in col..=n {
                values[j] = values[j] - factor * pivot_row[j];
            }
        }
    }
    Ok(matrix.into_iter().map(|row| row[n]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reedsolo_example() {
        // https://github.com/tomerfiliba-org/reedsolomon の使用例
        let rs = ReedSolomon::new(10).unwrap();
        let expected = b"hello world\xed\x25\x54\xc4\xfd\xfd\x89\xf3\xa8\xaa";
        assert_eq!(rs.encode(b"hello world").unwrap(), expected);
        assert_eq!(rs.decode(expected).unwrap(), b"hello world");
    }

    #[test]
    fn encoded_polynomial_has_generator_roots() {
        let rs = ReedSolomon::new(8).unwrap();
        let encoded = rs.encode(b"hello").unwrap();
        assert_eq!(&encoded[..5], b"hello");
        let polynomial = Polynomial::new(encoded.into_iter().map(Gf256).collect());
        for i in 0..8 {
            assert_eq!(polynomial.evaluate(Gf256(2).pow(i)), Gf256::ZERO);
        }
        let zeros = rs.encode(&[0; 12]).unwrap();
        assert_eq!(zeros, vec![0; 20]);
        assert_eq!(rs.decode(&zeros).unwrap(), vec![0; 12]);
    }

    #[test]
    fn correct_errors() {
        let mut seed = 42_u32;
        let mut next = || {
            seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
            seed
        };
        for parity in [1, 2, 3, 8, 16, 32, 254] {
            let rs = ReedSolomon::new(parity).unwrap();
            for len in [1, (255 - parity).min(11), 255 - parity] {
                let data: Vec<_> = (0..len).map(|_| next() as u8).collect();
                let encoded = rs.encode(&data).unwrap();
                assert_eq!(rs.decode(&encoded).unwrap(), data);
                let mut positions: Vec<_> = (0..encoded.len()).collect();
                for i in (1..positions.len()).rev() {
                    positions.swap(i, next() as usize % (i + 1));
                }
                let mut damaged = encoded;
                for &p in positions.iter().take(parity / 2) {
                    damaged[p] ^= (next() % 255 + 1) as u8;
                }
                assert_eq!(rs.decode(&damaged).unwrap(), data);
            }
        }
    }

    #[test]
    fn errors_with_zero_first_syndrome() {
        let rs = ReedSolomon::new(4).unwrap();
        let mut encoded = rs.encode(b"hello").unwrap();
        encoded[0] ^= 7;
        encoded[8] ^= 7;
        assert_eq!(rs.decode(&encoded).unwrap(), b"hello");

        let rs = ReedSolomon::new(2).unwrap();
        let mut encoded = rs.encode(b"hello").unwrap();
        encoded[0] ^= 7;
        encoded[1] ^= 7;
        assert_eq!(rs.decode(&encoded), Err(Error::Uncorrectable));
    }

    #[test]
    fn invalid_inputs() {
        assert!(matches!(ReedSolomon::new(0), Err(Error::InvalidParity)));
        assert!(matches!(ReedSolomon::new(255), Err(Error::InvalidParity)));
        let rs = ReedSolomon::new(8).unwrap();
        assert_eq!(rs.encode(&[]), Err(Error::InvalidLength));
        assert_eq!(rs.encode(&[0; 248]), Err(Error::InvalidLength));
        assert_eq!(rs.decode(&[0; 8]), Err(Error::InvalidLength));
        assert_eq!(rs.decode(&[0; 256]), Err(Error::InvalidLength));
    }
}
