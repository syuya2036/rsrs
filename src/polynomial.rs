use crate::field::Field;

// 係数は次数の高い順。ゼロ多項式は [0] に揃える。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial<F: Field> {
    coefficients: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    pub fn new(coefficients: Vec<F>) -> Self {
        let start = coefficients.iter().position(|&c| c != F::ZERO);
        let coefficients = match start {
            Some(i) => coefficients[i..].to_vec(),
            None => vec![F::ZERO],
        };
        Self { coefficients }
    }

    pub fn coefficients(&self) -> &[F] {
        &self.coefficients
    }

    pub fn degree(&self) -> Option<usize> {
        if self.coefficients == [F::ZERO] {
            None
        } else {
            Some(self.coefficients.len() - 1)
        }
    }

    pub fn evaluate(&self, x: F) -> F {
        self.coefficients.iter().fold(F::ZERO, |y, &c| y * x + c)
    }

    pub fn multiply(&self, rhs: &Self) -> Self {
        let mut result = vec![F::ZERO; self.coefficients.len() + rhs.coefficients.len() - 1];
        for (i, &a) in self.coefficients.iter().enumerate() {
            for (j, &b) in rhs.coefficients.iter().enumerate() {
                result[i + j] = result[i + j] + a * b;
            }
        }
        Self::new(result)
    }

    // ゼロ多項式で割る場合だけ None。
    pub fn remainder(&self, divisor: &Self) -> Option<Self> {
        let inverse = divisor.coefficients[0].inverse()?;
        let mut result = self.coefficients.clone();
        if result.len() >= divisor.coefficients.len() {
            for i in 0..=result.len() - divisor.coefficients.len() {
                let factor = result[i] * inverse;
                for (j, &c) in divisor.coefficients.iter().enumerate() {
                    result[i + j] = result[i + j] - factor * c;
                }
            }
        }
        Some(Self::new(result))
    }
}
