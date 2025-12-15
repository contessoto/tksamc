#!/usr/bin/env python3
"""
Plot script for TKSA-MC results.
Reads the CSV output file and generates a bar plot.
"""

import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description='Plot TKSA-MC results from CSV')
    parser.add_argument('csv_file', help='Input CSV file')
    args = parser.parse_args()

    csv_file = args.csv_file
    if not os.path.exists(csv_file):
        print(f"Error: File {csv_file} not found.")
        sys.exit(1)

    try:
        # Read CSV
        # Use line 0 as header if it starts with #
        # Pandas `comment` arg skips the line entirely, meaning we LOSE the header if it starts with #.
        # We should read the header separately or handle it.

        # Read the first line to check header
        with open(csv_file, 'r') as f:
            header_line = f.readline().strip()

        if header_line.startswith('#'):
            # It's our header but commented out.
            # We can use `names` argument or read with `header=0` but tell pandas to ignore the hash.
            # Easiest is to read, providing column names by cleaning the first line.

            col_names = header_line.lstrip('#').strip().split(',')
            col_names = [c.strip() for c in col_names]

            df = pd.read_csv(csv_file, comment='#', names=col_names)
        else:
            df = pd.read_csv(csv_file)

        # Check columns
        # Expected: 1-Name,2-Residue-index, ..., 12-dG_Energy
        # 12-dG_Energy is the value to plot.
        # We also want labels from 1-Name + 2-Residue-index

        # Clean column names (strip spaces if any)
        df.columns = [c.strip() for c in df.columns]

        energy_col = '12-dG_Energy'
        name_col = '1-Name'
        res_idx_col = '2-Residue-index'

        if energy_col not in df.columns:
            print(f"Error: Column {energy_col} not found in CSV.")
            print(f"Columns found: {df.columns}")
            sys.exit(1)

        plot_data = df[energy_col].values

        labels = []
        for index, row in df.iterrows():
            # Create label e.g., MET1
            # Assuming Name is 3-letter, convert to 1-letter if standard
            res = row[name_col]
            res_map = {
                'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
                'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
                'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
                'N_T':'NTR', 'C_T':'CTR', 'N_TER':'NTR', 'C_TER':'CTR'
            }
            res_code = res_map.get(res, res)
            labels.append(f"{res_code}{row[res_idx_col]}")

        # Plotting
        x_pos = range(len(plot_data))
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)

        colors = []
        for val in plot_data:
            if val > 0:
                colors.append('red')
            else:
                colors.append('blue')

        ax.bar(x_pos, plot_data, color=colors, width=1.0, edgecolor='black', linewidth=0.5)

        plt.xlabel('Residue')
        plt.ylabel(r'$\Delta G_{qq}$ (kJ/mol)')
        plt.title(f'Electrostatic Free Energy per Residue\n{os.path.basename(csv_file)}')

        # Ticks
        fontsize = 12
        if len(plot_data) > 35:
            fontsize = 8

        plt.xticks(x_pos, labels, rotation=90, fontsize=fontsize)
        plt.xlim([-0.5, len(x_pos)-0.5])

        # Save
        out_file = os.path.splitext(csv_file)[0] + '_plot.png'
        plt.tight_layout()
        plt.savefig(out_file, dpi=300)
        print(f"Plot saved to {out_file}")

    except Exception as e:
        print(f"Error plotting: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
