use rsrs::ReedSolomon;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let rs = ReedSolomon::new(8)?;
    let mut encoded = rs.encode(b"hello reed solomon")?;
    println!("encoded: {encoded:?}");

    encoded[0] ^= 0x12;
    encoded[5] ^= 0x34;
    encoded[10] ^= 0x56;
    encoded[20] ^= 0x78;

    let decoded = rs.decode(&encoded)?;
    println!("decoded: {}", String::from_utf8(decoded)?);
    Ok(())
}
