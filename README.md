TKSA-MC
=======

|Citing-TKSA| |PyPI| |Python| |License| |GitHub-Stars|

.. |Citing-TKSA| image:: https://img.shields.io/badge/cite-TKSA--MC-informational
   :target: https://doi.org/10.1101/221556

.. |PyPI| image:: https://img.shields.io/pypi/v/tksamc.svg
   :target: https://pypi.org/project/tksamc/

.. |Python| image:: https://img.shields.io/pypi/pyversions/tksamc.svg
   :target: https://pypi.org/project/tksamc/

.. |License| image:: https://img.shields.io/github/license/contessoto/tksamc.svg
   :target: https://github.com/contessoto/tksamc/blob/main/LICENSE

.. |GitHub-Stars| image:: https://img.shields.io/github/stars/contessoto/tksamc.svg?style=social
   :target: https://github.com/contessoto/tksamc

`GitHub <https://github.com/contessoto/tksamc>`__
| `PyPI <https://pypi.org/project/tksamc/>`__
| `Issues <https://github.com/contessoto/tksamc/issues>`__

Overview
========

**TKSA-MC** (Tanford–Kirkwood Surface Accessibility – Monte Carlo) is a Python 3 implementation of the
Tanford–Kirkwood electrostatic model for proteins. The code computes charge–charge interaction free
energies by explicitly sampling protonation states, using either an **exact solver** (for small systems)
or **Monte Carlo sampling** (for large proteins).

The method estimates the electrostatic free energy contribution,  
:math:`\Delta G_{qq}`, for each ionizable residue, accounting for solvent accessibility via
solvent-accessible surface area (SASA).

TKSA-MC is designed to be **simple, fast, and reproducible**, with support for command-line execution
and scripting workflows.


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
