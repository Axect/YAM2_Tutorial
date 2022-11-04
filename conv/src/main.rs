use peroxide::fuga::*;

fn main() -> Result<(), Box<dyn Error>> {
    let mut df = DataFrame::read_csv("../data/m2ccb.dat", ',')?;
    let mut dg = DataFrame::read_csv("../data/ttbar_m2ccb.dat", ',')?;

    df.as_types(vec![F64]);
    dg.as_types(vec![F64]);

    df.print();
    println!("");
    dg.print();

    df.write_nc("../data/m2ccb.nc")?;
    dg.write_nc("../data/ttbar_m2ccb.nc")?;

    Ok(())
}
