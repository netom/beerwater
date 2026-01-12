use rand::Rng;
use rand_distr::{Distribution, Normal};
pub use std::{
    fs::{self, File},
    io::{self, BufRead, BufReader, Lines},
    process::exit,
};

#[derive(Clone, Debug)]
enum Constraint {
    Exact(f32),           // Single target value
    Range(f32, f32),      // Min and max (inclusive)
    Unconstrained,        // Any value is acceptable
}

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

// Error function for concentrations against constraints
fn err(constraints: &Vec<Constraint>, concentrations: &Vec<f32>) -> f32 {
    let mut sum: f32 = 0.0;
    for i in 0..constraints.len() {
        let error_contribution = match &constraints[i] {
            Constraint::Exact(target) => {
                let diff = concentrations[i] - target;
                diff * diff
            }
            Constraint::Range(min, max) => {
                if concentrations[i] < *min {
                    let diff = concentrations[i] - min;
                    diff * diff
                } else if concentrations[i] > *max {
                    let diff = concentrations[i] - max;
                    diff * diff
                } else {
                    0.0 // Within range, no error
                }
            }
            Constraint::Unconstrained => 0.0, // No constraint, no error
        };
        sum += error_contribution;
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

fn file_lines_or_exit(file_name: &str, file_description: &str) -> Lines<BufReader<File>> {
    match File::open(file_name) {
        Ok(file) => return BufReader::new(file).lines(),
        Err(io_error) => {
            println!("Could not open {file_description} file {file_name}: {io_error}");
            exit(1);
        }
    }
}

fn process_data_file_or_exit<F: FnMut(Option<Vec<&str>>) -> Option<Result<(), String>>>(
    line_number: &mut u64,
    lines: &mut Lines<BufReader<File>>,
    file_name: &str,
    file_description: &str,
    mut processor: F,
) {
    for line_result in lines.by_ref() {
        *line_number += 1;
        let line;
        match line_result {
            Ok(line_) => line = line_,
            Err(io_error) => {
                println!(
                    "Could not read {} file {} at position {}: {}.",
                    file_description, file_name, line_number, io_error
                );
                exit(1);
            }
        }
        let line: &str = line.trim();
        if line.is_empty() {
            continue;
        }
        let line_parts: Vec<&str> = line.split_whitespace().collect();
        match processor(Some(line_parts)) {
            None => continue,
            Some(Ok(())) => return,
            Some(Err(error_message)) => {
                println!(
                    "Error reading {file_description} file {file_name} at line {line_number}: {error_message}."
                );
                exit(1);
            }
        }
    }
    match processor(None) {
        None => {
            println!("Unexpected end of {file_description} file {file_name}.")
        }
        Some(Ok(())) => return,
        Some(Err(error_message)) => {
            println!(
                "Error reading {file_description} file {file_name} at line {line_number}: {error_message}."
            );
            exit(1);
        }
    }
}

fn main() {
    let ion_contributions_file_name = "ion_contributions.txt";
    let ion_contributions_file_description = "ion contribution data";

    let mut ion_contribution_lines = file_lines_or_exit(
        ion_contributions_file_name,
        ion_contributions_file_description,
    );

    let mut ion_contributions_line_number: u64 = 0;

    let mut maybe_hco3_index = None; // Index of the "HCO3-" ion, if exists
    let mut ions: Vec<String> = Vec::new();
    process_data_file_or_exit(
        &mut ion_contributions_line_number,
        &mut ion_contribution_lines,
        ion_contributions_file_name,
        ion_contributions_file_description,
        |maybe_words| {
            let words;
            match maybe_words {
                Some(some_words) => words = some_words, // We've got some words!
                None => return None,                    // Unexpected end of file
            }

            /* "words" contains the first non-empty,
             *  non-comment line, these are our ions */
            let mut ion_index = 0;
            for word in words {
                ions.push(String::from(word));
                if word.eq("HCO3-") {
                    maybe_hco3_index = Some(ion_index);
                }
                ion_index += 1;
            }
            return Some(Ok(())); // We're done.
        },
    );

    let mut salts: Vec<String> = Vec::new();
    let mut contributions: Vec<Vec<f32>> = Vec::new();
    process_data_file_or_exit(
        &mut ion_contributions_line_number,
        &mut ion_contribution_lines,
        ion_contributions_file_name,
        ion_contributions_file_description,
        |maybe_words| {
            let words;
            match maybe_words {
                Some(some_words) => words = some_words, // We've got some words!
                None => return Some(Ok(())), // End of file, done reading ion contributions
            }

            /* "words" now contain non-empty, non-comment lines,
             * these are our ion contributions of salts. */
            if words.len() != ions.len() + 1 {
                return Some(Err(
                    "word count should be the number of ions plus one".to_string()
                ));
            }

            let salt = words[0].to_string();
            let mut contributions_for_this_salt: Vec<f32> = Vec::new();

            salts.push(salt);

            let mut field_counter = 1;
            for contribution in words.iter().skip(1) {
                field_counter += 1;
                match contribution.parse() {
                    Ok(value) => {
                        if value < 0.0 {
                            return Some(Err(format!(
                                "ion contribution at field {field_counter} is negative"
                            )));
                        }
                        contributions_for_this_salt.push(value);
                    }
                    Err(parse_error) => {
                        return Some(Err(format!(
                            "parse error at field {field_counter}: {parse_error}"
                        )))
                    }
                }
            }
            contributions.push(contributions_for_this_salt);

            return None; // We'd like to have some more please.
        },
    );

    // Size of step in each direction, g/l
    let eps: f32 = 0.0002;

    // Water quantity in litres
    let water_quantity: f32 = 25.0;

    let target_file_name = "target.txt";
    let target_file_description = "target concentration data";

    let mut target_lines = file_lines_or_exit(target_file_name, target_file_description);

    let mut target_line_number: u64 = 0;

    let mut constraints = vec![Constraint::Unconstrained; ions.len()];

    process_data_file_or_exit(
        &mut target_line_number,
        &mut target_lines,
        target_file_name,
        target_file_description,
        |maybe_words| {
            let words;
            match maybe_words {
                Some(some_words) => words = some_words, // We've got some words!
                None => return Some(Ok(())), // End of file, done reading ion contributions
            }

            if words.len() < 2 {
                return Some(Err("the number of fields must be at least 2".to_string()));
            }

            let ion: &str = words[0];

            let constraint = if words[1] == "*" {
                // Unconstrained
                Constraint::Unconstrained
            } else if words.len() == 4 && words[2] == "-" {
                // Range format: "ion min - max"
                let min_value: f32 = match words[1].parse() {
                    Ok(value) => {
                        if value < 0.0 {
                            return Some(Err("minimum value is negative".to_string()));
                        }
                        value
                    }
                    Err(parse_error) => {
                        return Some(Err(format!(
                            "error parsing minimum value: {parse_error}"
                        )));
                    }
                };

                let max_value: f32 = match words[3].parse() {
                    Ok(value) => {
                        if value < 0.0 {
                            return Some(Err("maximum value is negative".to_string()));
                        }
                        value
                    }
                    Err(parse_error) => {
                        return Some(Err(format!(
                            "error parsing maximum value: {parse_error}"
                        )));
                    }
                };

                if min_value > max_value {
                    return Some(Err("minimum value is greater than maximum value".to_string()));
                }

                Constraint::Range(min_value, max_value)
            } else if words.len() == 2 {
                // Exact value format: "ion value"
                let target_concentration: f32 = match words[1].parse() {
                    Ok(value) => {
                        if value < 0.0 {
                            return Some(Err("target concentration is negative".to_string()));
                        }
                        value
                    }
                    Err(parse_error) => {
                        return Some(Err(format!(
                            "error parsing target concentration: {parse_error}"
                        )));
                    }
                };
                Constraint::Exact(target_concentration)
            } else {
                return Some(Err("invalid constraint format (expected 'ion value', 'ion min - max', or 'ion *')".to_string()));
            };

            let ion_index;
            match ions.iter().position(|v| v == ion) {
                Some(value) => ion_index = value,
                None => return Some(Err(format!("unkown ion: {ion}"))),
            }

            constraints[ion_index] = constraint;

            return None;
        },
    );

    let mut best_concentrations: Vec<f32> = vec![0.0; ions.len()];
    let mut try_concentrations: Vec<f32> = vec![0.0; ions.len()];

    // Initial random quantities, 0 - 1 g/l
    let mut try_quantities: Vec<f32> = rand::rng()
        .sample_iter(rand::distr::StandardUniform)
        .take(salts.len())
        .collect();
    let mut best_quantities: Vec<f32> = try_quantities.clone();

    conc(&contributions, &best_quantities, &mut best_concentrations);

    let mut best_err: f32 = err(&constraints, &best_concentrations);

    for i in 1..500001 {
        nudge(eps, &best_quantities, &mut try_quantities);
        conc(&contributions, &try_quantities, &mut try_concentrations);

        let try_err = err(&constraints, &try_concentrations);

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
    println!("Target constraints:");
    println!("");
    for i in 0..ions.len() {
        let constraint_str = match &constraints[i] {
            Constraint::Exact(value) => format!("{}", value),
            Constraint::Range(min, max) => format!("{} - {}", min, max),
            Constraint::Unconstrained => "*".to_string(),
        };
        println!("{} {}", ions[i], constraint_str);
    }

    println!("");
    println!("Achieved concentrations:");
    println!("");
    for i in 0..ions.len() {
        let achieved = best_concentrations[i];
        let status = match &constraints[i] {
            Constraint::Exact(target) => {
                if (achieved - target).abs() < 1.0 {
                    "✓"
                } else {
                    "✗"
                }
            }
            Constraint::Range(min, max) => {
                if achieved >= *min && achieved <= *max {
                    "✓"
                } else {
                    "✗"
                }
            }
            Constraint::Unconstrained => "✓",
        };
        println!("{} {} {}", ions[i], achieved, status);
    }

    println!("");
    match maybe_hco3_index {
        Some(hco3_index) => {
            let result_alkalinity = alkalinity(&best_concentrations, hco3_index);
            println!("Alkalinity: {result_alkalinity:.3}");
        }
        None => {
            println!("No HCO3 ion present, cannot calculate alkalinity.")
        }
    }

    println!("");
    println!("Salt additions for {}l of water:", water_quantity);
    println!("");
    for s in 0..salts.len() {
        println!("{} {}", salts[s], best_quantities[s] * water_quantity);
    }
}
