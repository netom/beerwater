# Beer brewing water salt addition calculator

A command-line tool that calculates how much of each brewing salt to add to
water in order to hit a set of target ion concentrations. Instead of solving
the system analytically, it uses a stochastic (random-nudge) optimizer that
searches for a combination of salt quantities that minimizes the error
against your constraints.

## How it works

1. `ion_contributions.txt` defines the ions you care about and, for each
   available salt, how many mg/L of each ion it contributes per g/L of salt
   dissolved.
2. `target.txt` defines the constraints you want to hit for each ion:
   an exact value, an acceptable range, an ion-to-ion ratio, or "don't care".
3. The optimizer starts from random salt quantities and, over 500,000
   iterations, repeatedly perturbs them and keeps any change that reduces
   the total constraint error, until it converges on a good combination.
4. The result is printed as target vs. achieved concentrations (with ✓/✗
   status), achieved ratios, alkalinity, and the actual salt additions in
   grams for your batch of water.

## Requirements

- [Rust](https://www.rust-lang.org/tools/install) (stable toolchain, 2021
  edition)

## Usage

Run from the repository root so the tool can find `ion_contributions.txt`
and `target.txt` in the current directory:

```bash
cargo run --release
```

(`--release` is strongly recommended - 500,000 iterations is noticeably
slower in a debug build.)

### `ion_contributions.txt` format

- First non-empty line: space-separated ion names, e.g.

  ```
  Ca+2   Mg+2   Na+    SO4-2  Cl-    HCO3-
  ```

- Each following line: a salt name followed by its contribution (mg/L of
  ion per g/L of salt) for every ion listed above, in the same order:

  ```
  CaCO3   198.7    0.0    0.0    0.0    0.0  607.6
  NaHCO3  0.0      0.0  283.8    0.0    0.0  723.0
  CaSO4   232.8    0.0    0.0  558.0    0.0    0.0
  ```

### `target.txt` format

One constraint per line. Ion names must match those from
`ion_contributions.txt`. Ions not mentioned default to unconstrained.

| Format                   | Meaning                                    | Example           |
|--------------------------|--------------------------------------------|-------------------|
| `Ion value`              | Exact target concentration (mg/L)          | `Ca+2 80`         |
| `Ion min - max`          | Acceptable range (mg/L, inclusive)         | `Mg+2 8 - 12`     |
| `Ion *`                  | Unconstrained                              | `Na+ *`           |
| `Ion1 : Ion2 ratio`      | Target ratio of Ion1 to Ion2 concentration | `SO4-2 : Cl- 2.5` |

Example `target.txt`:

```
Ca+2    80.0
Mg+2    15.0
Na+     15.0
HCO3-    0.0
SO4-2 : Cl- 2.5
```

### Output

The program prints progress every 10,000 iterations, then a summary of achieved
ion concentrations, indicating wether relevant constraints were met.

```
ERR @10000: 5250.417
ERR @20000: 3704.5386
...
ERR @490000: 0.0000048575102
ERR @500000: 0.0000048575102

Target constraints:

Ca+2 80
Mg+2 15
Na+ 15
HCO3- 0
SO4-2 : Cl- 2.5

Achieved concentrations:

Ca+2 80.00045 ✓
Mg+2 14.999464 ✓
Na+ 15.001669 ✓
HCO3- 0 ✓
SO4-2 : Cl- = 183.02328 : 73.24617 = 2.498742 ✓

Alkalinity: 0.000

Salt additions for 25l of water:

CaCO3 0
NaHCO3 0
CaSO4 5.5371614
CaCl 2.6090276
MgSO4 3.8108394
NaCl 0.95260787
```

Salt additions are printed in grams for the water volume configured by
`water_quantity` in `src/main.rs` (25 L by default).

## Development

```bash
cargo build            # debug build
cargo build --release  # optimized build
cargo check            # fast type/borrow check without building a binary
cargo test              # run tests (none currently defined)
```

Everything lives in `src/main.rs`:

- `ion_contributions.txt` and `target.txt` are parsed with a shared
  `process_data_file_or_exit()` helper that handles line-by-line reading,
  EOF, and formatted error reporting.
- `Constraint` (`Exact`, `Range`, `Ratio`) encodes each parsed target line.
- `nudge()`, `conc()`, and `err()` implement the perturb -> recompute
  concentrations -> score-against-constraints loop used by the optimizer.

To experiment with a different water profile or salt list, edit
`ion_contributions.txt` and `target.txt` and rerun `cargo run --release`.
