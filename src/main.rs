use rand::random;
use rand_distr::{Distribution, Normal};

// Indexes of ions in the salt array
const I_CA: usize = 0;
const I_MG: usize = 1;
const I_NA: usize = 2;
const I_SO4: usize = 3;
const I_CL: usize = 4;
const I_HCO3: usize = 5;

// Ions

const NUM_IONS: usize = 6;
const ions: [&str; NUM_IONS] = ["Ca+2", "Mg+2", "Na+", "SO4-2", "Cl-", "HCO3-"];

// Salts and their ion contributions
// in ppm * g / l

const NUM_SALTS: usize = 6;

const salts: [&str; NUM_SALTS] = ["CaCO3", "NaHCO3", "CaSO4", "CaCl", "MgSO4", "NaCl"];

const contributions: [[f32; NUM_IONS]; NUM_SALTS] = [
    [198.7, 0.0, 0.0, 0.0, 0.0, 607.6],
    [0.0, 0.0, 283.8, 0.0, 0.0, 723.0],
    [232.8, 0.0, 0.0, 558.0, 0.0, 0.0],
    [272.5, 0.0, 0.0, 0.0, 480.7, 0.0],
    [0.0, 98.4, 0.0, 389.9, 0.0, 0.0],
    [0.0, 0.0, 393.7, 0.0, 605.7, 0.0],
];

// Size of step in each direction, g/l
const EPS: f32 = 0.0005;

fn nudge(qin: &[f32; NUM_SALTS], qout: &mut [f32; NUM_SALTS]) {
    let normal = Normal::new(-1.0 * EPS, 1.0 * EPS).unwrap();
    let mut rng = rand::thread_rng();

    for s in 0..NUM_SALTS {
        qout[s] = f32::max(0.0, qin[s] + normal.sample(&mut rng));
    }
}

// Alkalinity based on ion concentrations
fn alkalinity(concentrations: &[f32; NUM_IONS]) -> f32 {
    return concentrations[I_HCO3] * 50.0 / 61.0;
}

// Error function for concentrations
fn err(c1: &[f32; NUM_IONS], c2: &[f32; NUM_IONS]) -> f32 {
    let mut sum: f32 = 0.0;
    for i in 0..NUM_IONS {
        sum += f32::powf(c1[i] - c2[i], 2.0);
    }
    return sum;
}

// Return ion concentrations based on salt quantities
fn conc(quantities: &[f32; NUM_SALTS], concentrations: &mut [f32; NUM_IONS]) {
    for i in 0..NUM_IONS {
        concentrations[i] = 0.0;
        for s in 0..NUM_SALTS {
            concentrations[i] += contributions[s][i] * quantities[s];
        }
    }
}

// Water quantity in litres
const q: f32 = 30.0;

// Target concentrations
const target: [f32; NUM_IONS] = [54.0, 10.0, 10.0, 80.0, 80.0, 0.0];

fn main() {
    let mut best_concentrations: [f32; NUM_IONS] = [0.0; NUM_IONS];
    let mut try_concentrations: [f32; NUM_IONS] = [0.0; NUM_IONS];

    // Initial random quantities, 0 - 1 g/l
    let mut try_quantities: [f32; NUM_SALTS] = rand::random();
    let mut best_quantities: [f32; NUM_SALTS] = rand::random();
    conc(&best_quantities, &mut best_concentrations);
    let mut best_err: f32 = err(&target, &best_concentrations);

    for i in 1..500001 {
        nudge(&best_quantities, &mut try_quantities);
        conc(&try_quantities, &mut try_concentrations);

        let try_err = err(&target, &try_concentrations);

        if try_err < best_err {
            best_err = try_err;
            best_concentrations = try_concentrations;
            best_quantities = try_quantities;
        }

        if i % 10000 == 0 {
            println!("ERR @{}: {}", i, best_err);
        }
    }
    println!("");
    println!("Target concentrations:");
    println!("");
    for i in 0..NUM_IONS {
        println!("{} {}", ions[i], target[i]);
    }

    //best_quantities = np.round(best_quantities, decimals = 3);
    //best_concentrations = np.sum(contributions * best_quantities, axis = 0);
    //conc(&best_quantities, &mut best_concentrations);

    println!("");
    println!("Achieved concentrations:");
    println!("");
    for i in 0..NUM_IONS {
        println!("{} {}", ions[i], best_concentrations[i]);
    }

    println!("");
    println!("Salt additions for {}l of water:", q);
    println!("");
    for s in 0..NUM_SALTS {
        println!("{} {}", salts[s], best_quantities[s] * q);
    }
}
