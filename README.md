# TKSA-MC (Python 3)

**Tanford-Kirkwood Surface Accessibility - Monte Carlo**

This software calculates protein charge-charge interactions via the Tanford-Kirkwood Surface Accessibility model, using either an Exact method (for small systems) or Monte Carlo sampling (for larger systems) to determine protonation states and electrostatic free energies.

## Description

This is a Python 3 port of the original TKSA-MC code. It uses `mdtraj` for structure parsing and SASA calculation, and `numpy`/`numba` for efficient computation of electrostatic energies.

The method calculates the electrostatic free energy contribution ($\Delta G_{qq}$) for each ionizable residue.

## Requirements

*   Python 3.x
*   `numpy`
*   `scipy`
*   `matplotlib`
*   `mdtraj`
*   `numba` (for performance optimization of MC solver)

Install dependencies:

```bash
pip install numpy scipy matplotlib mdtraj numba
```

## Installation

You can install the package using pip:

```bash
pip install .
```

This will make the `tksamc` command available in your terminal.

## Usage

### Method 1: Installed Package

Run the `tksamc` command:

```bash
tksamc -f sample_pdb_1ubq.pdb -ph 7.0 -T 300.0 -s MC
```

### Method 2: Running Locally (No Install)

If you prefer not to install the package (e.g., for development or testing), you can use the `run_local.py` script provided in the root directory:

```bash
python3 run_local.py -f sample_pdb_1ubq.pdb -ph 7.0 -T 300.0 -s MC
```

Alternatively, run as a module:

```bash
python3 -m tksamc.cli -f sample_pdb_1ubq.pdb -ph 7.0 -T 300.0 -s MC
```

### Arguments

*   `-f`: Input PDB file.
*   `-ph`: pH value (default: 7.0).
*   `-T`: Temperature in Kelvin (default: 300.0).
*   `-s`: Solver method. Choices: `EX` (Exact) or `MC` (Monte Carlo). Default: `MC`.
*   `-e`: Electrostatic method. Default: `TK`.
*   `-plot`: Generate plot (`yes` or `no`). Default: `yes`.
*   `-aref`: Reference Max SASA set (`header` or `mdtraj`). `header` uses legacy values. `mdtraj` uses theoretical values for Bondi radii (Tien et al. 2013). Default: `header`.

### Output

*   **CSV File**: A CSV file (e.g., `Output_MC_sample_pdb_1ubq_pH_7.0_T_300.0.csv`) containing detailed results for each charged residue, including coordinates, pKa, SASA, Charge, and $\Delta G_{qq}$.
*   **Plot**: A JPG plot (e.g., `Fig_MC_sample_pdb_1ubq_pH_7.0_T_300.0.jpg`) showing the electrostatic free energy per residue. Red bars indicate positive values (unfavorable), and blue bars indicate negative values (favorable).

## References

> "TKSA-MC: A Web Server for rational mutation through the optimization of protein charge interactions -
> Vinícius G Contessoto, Vinícius M de Oliveira, Bruno R Fernandes, Gabriel G Slade, Vitor B. P. Leite, bioRxiv 221556; doi: https://doi.org/10.1101/221556"
