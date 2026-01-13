use rand::Rng;
use rand_distr::{Distribution, Normal};
pub use std::{
    fs::{self, File},
    io::{self, BufRead, BufReader, Lines},
    process::exit,
};

#[derive(Clone, Debug)]
enum Constraint {
    Exact(usize, f32),                 // ion_index, target value
    Range(usize, f32, f32),            // ion_index, min, max (inclusive)
    Ratio(usize, usize, f32),          // ion1_index, ion2_index, target_ratio (ion1/ion2)
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

    for constraint in constraints {
        let error_contribution = match constraint {
            Constraint::Exact(ion_index, target) => {
                let diff = concentrations[*ion_index] - target;
                diff * diff
            }
            Constraint::Range(ion_index, min, max) => {
                let conc = concentrations[*ion_index];
                if conc < *min {
                    let diff = conc - min;
                    diff * diff
                } else if conc > *max {
                    let diff = conc - max;
                    diff * diff
                } else {
                    0.0 // Within range, no error
                }
            }
            Constraint::Ratio(ion1_index, ion2_index, target_ratio) => {
                let ion1_conc = concentrations[*ion1_index];
                let ion2_conc = concentrations[*ion2_index];

                // If denominator is very small, penalize heavily
                if ion2_conc < 0.01 {
                    10000.0
                } else {
                    let actual_ratio = ion1_conc / ion2_conc;
                    let diff = actual_ratio - target_ratio;
                    diff * diff
                }
            }
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

    let mut constraints: Vec<Constraint> = Vec::new();

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

            // Check if this is a ratio constraint: "ion1 : ion2 value"
            if words.len() == 4 && words[1] == ":" {
                let ion1: &str = words[0];
                let ion2: &str = words[2];

                let ion1_index;
                match ions.iter().position(|v| v == ion1) {
                    Some(value) => ion1_index = value,
                    None => return Some(Err(format!("unknown ion: {ion1}"))),
                }

                let ion2_index;
                match ions.iter().position(|v| v == ion2) {
                    Some(value) => ion2_index = value,
                    None => return Some(Err(format!("unknown ion: {ion2}"))),
                }

                let target_ratio: f32 = match words[3].parse() {
                    Ok(value) => {
                        if value < 0.0 {
                            return Some(Err("target ratio is negative".to_string()));
                        }
                        value
                    }
                    Err(parse_error) => {
                        return Some(Err(format!(
                            "error parsing target ratio: {parse_error}"
                        )));
                    }
                };

                constraints.push(Constraint::Ratio(ion1_index, ion2_index, target_ratio));
                return None;
            }

            let ion: &str = words[0];

            let ion_index;
            match ions.iter().position(|v| v == ion) {
                Some(value) => ion_index = value,
                None => return Some(Err(format!("unknown ion: {ion}"))),
            }

            if words[1] == "*" {
                // Unconstrained - don't add any constraint
                return None;
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

                constraints.push(Constraint::Range(ion_index, min_value, max_value));
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
                constraints.push(Constraint::Exact(ion_index, target_concentration));
            } else {
                return Some(Err("invalid constraint format (expected 'ion value', 'ion min - max', 'ion *', or 'ion1 : ion2 ratio')".to_string()));
            }

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

    // Display ion constraints
    for constraint in &constraints {
        match constraint {
            Constraint::Exact(ion_index, target) => {
                println!("{} {}", ions[*ion_index], target);
            }
            Constraint::Range(ion_index, min, max) => {
                println!("{} {} - {}", ions[*ion_index], min, max);
            }
            Constraint::Ratio(_, _, _) => {
                // Skip ratio constraints for now, display them separately
            }
        }
    }

    // Display ratio constraints
    for constraint in &constraints {
        if let Constraint::Ratio(ion1_index, ion2_index, target_ratio) = constraint {
            println!("{} : {} {}", ions[*ion1_index], ions[*ion2_index], target_ratio);
        }
    }

    println!("");
    println!("Achieved concentrations:");
    println!("");

    // Display achieved ion concentrations
    for constraint in &constraints {
        match constraint {
            Constraint::Exact(ion_index, target) => {
                let achieved = best_concentrations[*ion_index];
                let status = if (achieved - target).abs() < 1.0 {
                    "✓"
                } else {
                    "✗"
                };
                println!("{} {} {}", ions[*ion_index], achieved, status);
            }
            Constraint::Range(ion_index, min, max) => {
                let achieved = best_concentrations[*ion_index];
                let status = if achieved >= *min && achieved <= *max {
                    "✓"
                } else {
                    "✗"
                };
                println!("{} {} {}", ions[*ion_index], achieved, status);
            }
            Constraint::Ratio(_, _, _) => {
                // Skip ratio constraints, display them separately
            }
        }
    }

    // Display achieved ratios
    for constraint in &constraints {
        if let Constraint::Ratio(ion1_index, ion2_index, target_ratio) = constraint {
            let ion1_conc = best_concentrations[*ion1_index];
            let ion2_conc = best_concentrations[*ion2_index];
            let achieved_ratio = if ion2_conc > 0.01 {
                ion1_conc / ion2_conc
            } else {
                f32::INFINITY
            };
            let status = if (achieved_ratio - target_ratio).abs() < 0.1 {
                "✓"
            } else {
                "✗"
            };
            println!("{} : {} {} {}", ions[*ion1_index], ions[*ion2_index], achieved_ratio, status);
        }
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
