#!/usr/bin/env python3
"""
Wrapper script to run TKSA-MC locally without installing the package.
Usage:
    python3 run_local.py -f <pdb> ...
"""
import sys
import os

# Ensure the current directory is in the python path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from tksamc.cli import main
except ImportError as e:
    print(f"Error importing tksamc: {e}")
    print("Make sure you have installed the dependencies: numpy, scipy, matplotlib, mdtraj, numba")
    sys.exit(1)

if __name__ == "__main__":
    main()
