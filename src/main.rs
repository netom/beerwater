use rand::Rng;
use rand_distr::{Distribution, Normal};
pub use std::fs::{self, File};
pub use std::io::{self, BufRead, BufReader};

fn nudge(eps: f32, qin: &Vec<f32>, qout: &mut Vec<f32>) {
    let normal = Normal::new(-1.0 * eps, 1.0 * eps).unwrap();
    let mut rng = rand::rng();

    for s in 0..qin.len() {
        qout[s] = f32::max(0.0, qin[s] + normal.sample(&mut rng));
    }
}

// Alkalinity based on ion concentrations
fn alkalinity(concentrations: &Vec<f32>, i_hco3: usize) -> f32 {
    return concentrations[i_hco3] * 50.0 / 61.0;
}

// Error function for concentrations
fn err(c1: &Vec<f32>, c2: &Vec<f32>) -> f32 {
    let mut sum: f32 = 0.0;
    for i in 0..c1.len() {
        // TODO: bounds?
        let diff = c1[i] - c2[i];
        sum += diff * diff;
    }
    return sum;
}

// Return ion concentrations based on salt quantities
fn conc(contributions: &Vec<Vec<f32>>, quantities: &Vec<f32>, concentrations: &mut Vec<f32>) {
    for i in 0..quantities.len() {
        concentrations[i] = 0.0;
        for s in 0..quantities.len() {
            concentrations[i] += contributions[s][i] * quantities[s];
        }
    }
}

fn main() {
    let file = File::open("ion_contributions.txt").unwrap();

    // "lines" will contain all non-empty lines split at whitespace sequences
    let lines: Vec<Vec<String>> = BufReader::new(file)
        .lines()
        .map(|line: Result<String, io::Error>| line.expect("Could not read"))
        .map(|line: String| String::from(line.trim()))
        .filter(|line: &String| !line.is_empty())
        .map(|line: String| line.split_whitespace().map(|s| String::from(s)).collect())
        .collect();

    let ions: Vec<String>;
    let mut salts: Vec<String> = Vec::new();
    let mut contributions: Vec<Vec<f32>> = Vec::new();

    ions = lines[0].clone();

    for i in 1..lines.len() {
        salts.push(lines[i][0].clone());
        contributions.push(
            lines[i]
                .iter()
                .skip(1)
                .map(|v| v.parse().expect("Not a number"))
                .collect(),
        );
    }

    drop(lines);

    // Size of step in each direction, g/l
    let eps: f32 = 0.0002;

    // Water quantity in litres
    let water_quantity: f32 = 25.0;

    // Target concentrations
    let target = vec![60.0, 10.0, 0.0, 70.0, 80.0, 00.0];

    let mut best_concentrations: Vec<f32> = vec![0.0; ions.len()];
    let mut try_concentrations: Vec<f32> = vec![0.0; ions.len()];

    // Initial random quantities, 0 - 1 g/l
    let mut try_quantities: Vec<f32> = rand::rng()
        .sample_iter(rand::distr::StandardUniform)
        .take(salts.len())
        .collect();
    let mut best_quantities: Vec<f32> = try_quantities.clone();

    conc(&contributions, &best_quantities, &mut best_concentrations);

    let mut best_err: f32 = err(&target, &best_concentrations);

    for i in 1..400001 {
        nudge(eps, &best_quantities, &mut try_quantities);
        conc(&contributions, &try_quantities, &mut try_concentrations);

        let try_err = err(&target, &try_concentrations);

        if try_err < best_err {
            best_err = try_err;
            best_concentrations.copy_from_slice(try_concentrations.as_slice());
            best_quantities.copy_from_slice(try_quantities.as_slice());
        }

        if i % 10000 == 0 {
            println!("ERR @{}: {}", i, best_err);
        }
    }
    println!("");
    println!("Target concentrations:");
    println!("");
    for i in 0..ions.len() {
        println!("{} {}", ions[i], target[i]);
    }

    //best_quantities = np.round(best_quantities, decimals = 3);
    //best_concentrations = np.sum(contributions * best_quantities, axis = 0);
    //conc(&best_quantities, &mut best_concentrations);

    println!("");
    println!("Achieved concentrations:");
    println!("");
    for i in 0..ions.len() {
        println!("{} {}", ions[i], best_concentrations[i]);
    }

    println!("");
    println!("Salt additions for {}l of water:", water_quantity);
    println!("");
    for s in 0..salts.len() {
        println!("{} {}", salts[s], best_quantities[s] * water_quantity);
    }
}
