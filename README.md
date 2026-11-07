# rsrs

Reed Solomon Codeを手で書いてみる。

- `field.rs`: 体の演算と、GF(256) の実装
- `polynomial.rs`: 体上の多項式。係数は次数の高い順
- `reed_solomon.rs`: 符号化と復号。誤り位置はBerlekamp–Masseyと根の探索で求め、誤り量は連立方程式を解く

## サンプルコード

```rust
use rsrs::ReedSolomon;

let rs = ReedSolomon::new(8).unwrap(); // パリティ8バイト
let mut encoded = rs.encode(b"hello").unwrap();
encoded[0] ^= 42;
assert_eq!(rs.decode(&encoded).unwrap(), b"hello");
```

GF(256)の法多項式は`0x11dで、生成元は2

生成多項式の根は`alpha^0`から`alpha^(parity-1)`。
