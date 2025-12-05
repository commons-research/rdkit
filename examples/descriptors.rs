use std::collections::HashMap;

use rdkit::{Properties, ROMol};

fn stereochem_status(properties: &HashMap<String, f64>) -> &'static str {
    let total = properties
        .get("NumAtomStereoCenters")
        .copied()
        .unwrap_or_default()
        .round() as i32;
    let unspecified = properties
        .get("NumUnspecifiedAtomStereoCenters")
        .copied()
        .unwrap_or_default()
        .round() as i32;

    match (total, unspecified) {
        (0, _) => "achiral — no stereocenters present",
        (_, n) if n == 0 => "fully specified stereochemistry",
        (t, u) if t == u => "unspecified stereochemistry",
        _ => "partially specified stereochemistry",
    }
}

fn main() {
    let smiles_samples = [
        ("C[C@H](F)Cl", "fully specified single center"),
        ("CC(F)(Cl)Br", "no stereochemistry assigned"),
        (
            "C[C@H](F)C(Br)(Cl)I",
            "one specified center and one unspecified",
        ),
        ("c1ccccc1", "achiral"),
    ];

    let properties = Properties::new();

    for (smiles, description) in smiles_samples {
        let mol = ROMol::from_smiles(smiles).expect("valid SMILES");
        let computed = properties.compute_properties(&mol);
        let total = computed
            .get("NumAtomStereoCenters")
            .copied()
            .unwrap_or_default();
        let unspecified = computed
            .get("NumUnspecifiedAtomStereoCenters")
            .copied()
            .unwrap_or_default();

        println!("SMILES: {smiles} ({description})");
        println!(
            "  stereochemistry: {} (total={total}, unspecified={unspecified})",
            stereochem_status(&computed)
        );
        println!(
            "  descriptors: amw={:.3}, CrippenClogP={:.3}",
            computed.get("amw").copied().unwrap_or_default(),
            computed.get("CrippenClogP").copied().unwrap_or_default()
        );
        println!();
    }
}
