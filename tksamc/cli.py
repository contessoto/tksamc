#!/usr/bin/env python3
#coding: utf8

__description__ = \
"""
TKSA - Electrostatic Free Energy calculation for each ionizable residue
"""

__author__ = "Vinícius G. Contessoto"
__date__ = "21/12/2016"

################################################################
#
# Version 1.0 (Python 3 Port)
#
# python tksamc.py -h # for help
#
# The following programs are provided free of charge only for academic use.
# By downloading these programs you implicitly agree that they will be used exclusively in academic research.
#
################################################################

import sys
import numpy as np
import scipy as sc
import subprocess
import os
import argparse
import time
import matplotlib.pyplot as plt
from itertools import chain
from scipy.spatial import distance
from scipy.special import eval_legendre
from numpy import linalg as LA
from itertools import islice
from subprocess import call
import profile
import threading
import mdtraj as md
from . import solver
import math
from scipy.special import eval_legendre

# Constants
All_residues = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
# Reference values for Max SASA (Header / Legacy)
Area_residues = [113,140,151,183,218,85,194,182,211,180,204,158,143,189,241,122,146,160,259,229]

# Reference values for Max SASA (MDTraj / Bondi - Tien et al. 2013)
# ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR
Area_residues_Tien2013 = [121.0, 148.0, 187.0, 214.0, 228.0, 97.0, 216.0, 195.0, 230.0, 191.0, 203.0, 187.0, 154.0, 214.0, 265.0, 143.0, 163.0, 165.0, 264.0, 255.0]

Charged_residues = ['ARG','LYS','N_TER','HIS','GLU','ASP','C_TER']
Charge_values = [0,0,0,0,-1,-1,-1]
Charged_atoms = ['NH2','NZ','NE2','OE2','OD2']
PKA = [12.0,10.6,7.7,6.3,4.5,4.0,3.6]
e = (4.0,78.5,1.0,6.02,1.602,8.8541878176,1.3806488,8.314) #(ep,es,x,mol,e,eo,K,R) reduced units?

##################################################################################################
# Kirkwood Polynomial Function
##################################################################################################
def Kn(n,x):
   Kpf = np.sum([np.power(x,s)*np.divide(np.power(2.0,s)*math.factorial(n)*math.factorial(2*n-s),math.factorial(s)*math.factorial(2*n)*math.factorial(n-s)) for s in range(n)])
   return Kpf
##################################################################################################

def main():
   global Q,E,S,Pk,e,T,pH,total_charged_residues,G,G0,indiv_data,Gqq,Q0

   parser = argparse.ArgumentParser(description='Charge-charge energy calculation in python')
   parser.add_argument('-ph', action='store', default=7.0, dest='arg_pH', help='pH value')
   parser.add_argument('-T', action='store', default=300.0, dest='arg_T',  help='Temperature value')
   parser.add_argument('-f', metavar='input-file-PDB',help='insert a PDB file',type=argparse.FileType('rt'))
   parser.add_argument('-e', action='store',choices=['TK'], default="TK",dest='arg_e',type=str,help='Electrostatic energy calculation method')
   parser.add_argument('-s', action='store',choices=['EX','MC'], default="MC",dest='arg_s',type=str,help='Statistical method to protonation state amostration - EX = Exact; MC = Monte Carlo;')
   parser.add_argument('-plot', action='store',choices=['yes','no'], default="yes",dest='arg_plot',type=str,help='Save Plot figure file - EPS')
   parser.add_argument('-aref', action='store', choices=['header', 'mdtraj'], default='header', dest='arg_aref', type=str, help='Reference Max SASA set. header=Legacy (Richards), mdtraj=Bondi (Tien 2013). Default: header')

   try:
       arguments = parser.parse_args()
   except IOError as msg:
       parser.error(str(msg))

   print('################################################')
   print(u"\U0001F63A", "### TKSA started ###", u"\U0001F63A")
   print('Input file:', arguments.f.name)
   print('pH  =', arguments.arg_pH)
   print('T   =', arguments.arg_T)
   print('Elec. Energy Calc Method =', arguments.arg_e)
   print('Statistical Method =', arguments.arg_s)
   print('Plot =', arguments.arg_plot)

   file_pdb_name = arguments.f.name
   pH = np.float64(arguments.arg_pH)
   T = np.float64(arguments.arg_T)

   ##################################################################################################
   # SASA Calculation using MDTraj
   ##################################################################################################
   print('Running SASA - MDTraj (Shrake-Rupley)')

   try:
       traj = md.load(file_pdb_name)
       # Calculate SASA
       # mode='residue' returns area per residue
       sasa_by_residue = md.shrake_rupley(traj, mode='residue', probe_radius=0.01)[0] * 100.0 # Convert nm^2 to Angstrom^2

       # Select Reference Area List
       if arguments.arg_aref == 'mdtraj':
           print("Using MDTraj (Bondi/Tien2013) Reference Areas for Normalization.")
           current_area_refs = Area_residues_Tien2013
       else:
           print("Using Header (Legacy/Richards) Reference Areas for Normalization.")
           current_area_refs = Area_residues

       # Mapping MDTraj residues to the list
       SASA_data = []
       for i, residue in enumerate(traj.topology.residues):
           res_name = residue.name
           # MDTraj residue names might differ (e.g. HIE, HID -> HIS)
           # Standardizing to 3-letter code if possible
           if res_name == 'HIE' or res_name == 'HID' or res_name == 'HIP': res_name = 'HIS'
           if res_name == 'CYX': res_name = 'CYS'

           if res_name not in All_residues:
               continue

           sasa_val = sasa_by_residue[i]
           ref_area = current_area_refs[All_residues.index(res_name)]
           Area_norm = sasa_val / ref_area

           if Area_norm >= 1.0:
               print(f"Warning - ** SASA greater than 1.0 ** {res_name} {residue.index} {sasa_val} {ref_area} {Area_norm}")
               print("Automatically changed to 0.75")
               Area_norm = 0.750000000001

           SASA_data.append([res_name, sasa_val, Area_norm, residue.index])

   except Exception as e:
       print(f'Error calculating SASA with MDTraj: {e}')
       sys.exit(1)

   # Identify charged residues using MDTraj topology and coordinates
   S = []
   total_charged_residues = []

   topology = traj.topology
   # We need coordinates in Angstroms for TK calculation (compatible with original code)
   # MDTraj uses nanometers.
   xyz = traj.xyz[0] * 10.0

   # Create a map of residue index to SASA data
   # SASA_data list contains [res_name, sasa_val, Area_norm, res_index]
   # Map by res_index (topology index)
   sasa_map = {item[3]: item for item in SASA_data}

   for residue in topology.residues:
       if residue.chain.index != 0:
           continue # Only Chain A (index 0) as per original restriction "chain == 'A'"

       res_name = residue.name
       # Normalize name
       if res_name in ['HIE', 'HID', 'HIP']: res_name = 'HIS'
       if res_name == 'CYX': res_name = 'CYS'

       res_index = residue.index
       pdb_res_index = residue.resSeq # This corresponds to column 6 in PDB

       # Get SASA
       sasa_info = sasa_map.get(res_index)
       if not sasa_info:
           # Should not happen if filtered correctly
           sasa_val = 0.5
       else:
           sasa_val = sasa_info[2] # Area_norm

       # Check for N-Terminus
       # Logic: If it's the first residue of the chain?
       # Original code checked "atom_index == 1 and atom_type == 'N'".
       # This implies N-term is charged ONLY if it's the very first atom of the file?
       # Let's assume N-term charge applies to the first residue of the chain.
       # BUT, original code: `if ... and not residue_type in Charged_residues`.
       # This means if the first residue IS charged (e.g. Lys), it gets its SIDE CHAIN charge,
       # but does it also get N-term charge?
       # Original code:
       # Block 1: `if atom_index == 1 ... N ... not in Charged_residues`: Add N_T.
       # Block 2: `if residue_type in Charged_residues ...`: Add side chain charge.
       # This implies if residue 1 is LYS, it gets LYS charge but NOT N_T charge?
       # That seems chemically wrong (it should have both), but I must reproduce code behavior.
       # Actually, `sample_pdb_1ubq.pdb` starts with MET (residue 1). MET is not charged.
       # So it gets N_T.
       # If it started with LYS, the original code would SKIP N_T check.
       # This might be a bug in original code or intended simplified model.
       # I will reproduce the logic: "If first residue AND not a charged type -> Add N_T".

       is_first_residue = (residue.index == 0) # Assuming chain A starts at 0

       if is_first_residue and res_name not in Charged_residues:
           # Add N_T
           # We need coordinates of N atom
           n_atom = next((a for a in residue.atoms if a.name == 'N'), None)
           if n_atom:
               total_charged_residues.append(n_atom.index) # Using atom index
               # Structure of S:
               # ['N_T', residue_index, count, atom_serial, atom_type, X, Y, Z, PKA, SASA, Charge]
               # Column 4: atom_serial (was atom_name, user requested correction)

               pos = xyz[n_atom.index]
               # n_atom.serial is available in MDTraj atoms
               S.append(['N_T', pdb_res_index, len(total_charged_residues), n_atom.serial, 'N',
                         pos[0], pos[1], pos[2],
                         PKA[Charged_residues.index('N_TER')], sasa_val, Charge_values[Charged_residues.index('N_TER')]])

       # Charged Side Chains
       if res_name in Charged_residues:
           # Find the charged atom
           target_atom_names = Charged_atoms
           target_atom = next((a for a in residue.atoms if a.name in target_atom_names), None)

           if target_atom:
               total_charged_residues.append(target_atom.index)
               pos = xyz[target_atom.index]
               S.append([res_name, pdb_res_index, len(total_charged_residues), target_atom.serial, target_atom.name,
                         pos[0], pos[1], pos[2],
                         PKA[Charged_residues.index(res_name)], sasa_val, Charge_values[Charged_residues.index(res_name)]])

       # C-Terminus
       # Check for OXT atom
       oxt_atom = next((a for a in residue.atoms if a.name == 'OXT'), None)
       if oxt_atom and res_name not in Charged_residues:
           total_charged_residues.append(oxt_atom.index)
           pos = xyz[oxt_atom.index]
           S.append(['C_T', pdb_res_index, len(total_charged_residues), oxt_atom.serial, 'OXT',
                     pos[0], pos[1], pos[2],
                     PKA[Charged_residues.index('C_TER')], sasa_val, Charge_values[Charged_residues.index('C_TER')]])

   print("There are: %d Charged_residues" % len(total_charged_residues))

   if not S:
       print("No charged residues found. Check PDB format (Chain A, standard residues).")
       sys.exit(0)

   Restype=np.asarray([i[0] for i in S])
   X=np.asarray([i[5] for i in S])
   Y=np.asarray([i[6] for i in S])
   Z=np.asarray([i[7] for i in S])
   Pk=np.asarray([i[8] for i in S])
   SA=np.asarray([i[9] for i in S])
   Q=np.asarray([i[10] for i in S])

   # Replace residue names with single letters
   Restype = np.char.replace(Restype, 'HIS', 'H')
   Restype = np.char.replace(Restype, 'ASP', 'D')
   Restype = np.char.replace(Restype, 'ARG', 'R')
   Restype = np.char.replace(Restype, 'GLU', 'E')
   Restype = np.char.replace(Restype, 'LYS', 'K')

   # Center coordinates
   X = X - np.mean(X)
   Y = Y - np.mean(Y)
   Z = Z - np.mean(Z)
   XYZ = list(zip(X,Y,Z))
   Origin = np.zeros(np.shape(XYZ))

   dist = distance.cdist(XYZ, XYZ, 'euclidean')

   # TK Calculation
   if arguments.arg_e == 'TK':
        dist_origin = distance.cdist(XYZ, Origin, 'euclidean')
        angle = distance.cdist(XYZ, XYZ, 'cosine')

        # Radii (b, a)
        # raio = (b, a)
        raio = (np.max(dist)*0.5 + 3.4+2.0, np.max(dist)*0.5 + 2.0+2.0)

        np.seterr(invalid='ignore')
        np.seterr(divide='ignore')

        theta = np.arccos(1-angle)

        NormA = np.matrix([LA.norm(v) for v in np.array(XYZ)])
        rirj = np.array(np.dot(np.transpose(NormA),NormA))

        # Avoid division by zero in rirj/raio^2 if raio is huge? No, raio is large.

        A = np.divide(raio[1],e[0]*dist)

        # B calculation
        # Sum over n from 0 to 60
        # This can be slow.
        B_sum = np.zeros_like(dist)
        for n in range(0, 60):
            term = ((e[1]-e[0])/(e[1]-(n*e[0])/(n+1)))
            term *= (np.power((rirj/(raio[1]*raio[1])),n))
            term *= (eval_legendre(n, np.cos(theta)))
            B_sum += term

        B = B_sum / e[0]

        # C calculation
        C_sum = np.zeros_like(dist)
        for n in range(1, 60):
            term1_num = (2*n+1) * e[1]
            term1_den = (2*n-1) * ((n+1)*e[1]+n*e[0])
            term1 = (term1_num / term1_den)
            term1 *= (np.power((rirj/(raio[0]*raio[0])),n))
            term1 *= (eval_legendre(n, np.cos(theta)))

            term2_den1 = Kn(n+1,e[2]) / Kn(n-1,e[2])
            term2_den2 = (n*(e[1]-e[0])) / ((n+1)*e[1]+n*e[0])
            term2_den2 *= (e[2]**2) / (4.0*(n**2)-1)
            term2_den2 *= np.power(raio[1]/raio[0], 2*n+1)

            term2 = term2_den1 + term2_den2

            C_sum += term1 / term2

        term_pre_C = (e[2] / (1+e[2]))
        C = (term_pre_C + (e[2]**2)*C_sum) / e[1]

        Qe = np.divide(e[3]*e[4]*e[4]*np.power(10,7),4*np.pi*e[5])

        # SAij: average SASA of pair?
        # lambda u,v: (u+v)*0.5
        SA_zip = list(zip(SA))
        SAij_mat = distance.cdist(SA_zip, SA_zip, lambda u,v: (u[0]+v[0])*0.5)

        E = Qe*(np.divide(A-B,2*raio[1])-np.divide(C,2*raio[0]))*(1-SAij_mat)

        # Clean E
        E[np.isinf(E)]= 0
        E[np.isnan(E)]= 0
        if np.sum(np.where(E<0)) > 0:
            print('###############################################################')
            print("There are: %d negatives TK energy values - Please check the radius of TK method!" % int(np.sum(np.where(E<0))))
            print("Sugestion - Increase b radius")
            print("Current radius ratio b/a=", np.divide(raio[1],raio[0]))
            print('###############################################################')

        # E calculation complete
        # No longer saving E.dat as requested

   # Solving
   plot_data = []

   if arguments.arg_s == 'EX':
       print(u"\U0001F63A", "### TK - Exact ###", u"\U0001F63A")
       start = time.time()

       Gqq_result = solver.solve_exact(E, Q, Pk, pH, T)

       # Convert to kJ/mol? C code did: Gqq[b]/1000.0
       # My solve_exact returns Energy in Joules/mol (because terms are in RT).
       # Wait.
       # My `solve_exact` uses `RT`. `RT` is in J/mol.
       # The return value is `Gqq`.
       # Is Gqq in J/mol or RT units?
       # `Gqq += ... * Interaction`. Interaction is in Energy units (from Eij).
       # So Gqq is in J/mol.
       # To get kJ/mol, divide by 1000.

       plot_data = Gqq_result / 1000.0

       end = time.time()
       elapsed = end - start
       print("Ran in %f sec" % elapsed)

   elif arguments.arg_s == 'MC':
       print(u"\U0001F63A", "### TKSA - MC ###", u"\U0001F63A")
       start = time.time()

       G_result = solver.solve_mc(E, Q, Pk, pH, T)

       plot_data = G_result # Already formatted by solve_mc to match C output logic

       end = time.time()
       elapsed = end - start
       print("Ran in %f sec" % elapsed)

   # Plotting / Saving Results
   if arguments.arg_plot == 'yes':
       # Prepare data
       S_arr = np.array(S, dtype=object)

       # Add result to S
       # S has 11 columns currently.
       # We append plot_data (dG) to it.
       # plot_data is 1D array of size N.

       # Format residue names for plot
       # Restype is already single letter.
       # "Name+Index"

       # X-axis labels
       labels = []
       for i in range(len(S)):
           labels.append(f"{Restype[i]}{S[i][1]}")

       # Total dG
       total_dG = np.sum(plot_data)
       print("Total dG Energy: ", total_dG)

       # Plot
       x_pos = np.arange(len(plot_data))
       fig = plt.figure(figsize=(10, 6)) # Improved size
       ax = fig.add_subplot(111)
       width=1.0
       colors = []
       for position, value in enumerate(plot_data):
           if value > 0 and SA[position] > 0.5:
               colors.append('red')
           else:
               colors.append('blue')

       # Improved bar aesthetics
       ax.bar(x_pos, plot_data, width=width, color=colors, edgecolor='black', linewidth=0.5)

       ax.tick_params('both', length=5, width=2, which='major', labelsize=12)
       plt.setp(ax.spines.values(), linewidth=1.5)

       # Adjust x-ticks based on size
       fontsize = 12
       if len(plot_data) > 35:
           fontsize = 8

       # Labels and Title
       plt.xlabel('Residue', fontsize=14)
       plt.ylabel(r'$\Delta G_{qq}$ (kJ/mol)', fontsize=14)
       plt.title(f'Electrostatic Free Energy per Residue\n{os.path.basename(file_pdb_name)}', fontsize=14)

       plt.xticks(x_pos, labels, rotation=90, fontsize=fontsize)
       plt.xlim([-0.5, len(x_pos)-0.5])

       plt.tight_layout()

       fig_filename = 'Fig_'+arguments.arg_s+'_'+ os.path.splitext(os.path.basename(file_pdb_name))[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.jpg'
       fig.savefig(fig_filename, dpi=300)
       # plt.show() # Can't show in headless env

       # Save Output CSV
       # We need to construct S with added column

       # S is list of lists.
       for i in range(len(S)):
           S[i].append(plot_data[i])

       # Updated header with commas
       header='1-Name,2-Residue-index,3-Position,4-Atom-Index,5-Atom-type,6-X,7-Y,8-Z,9-PKA,10-SASA,11-Charge,12-dG_Energy,13-Total_dG='+str(total_dG)+''

       out_filename = 'Output_'+arguments.arg_s+'_'+os.path.splitext(os.path.basename(file_pdb_name))[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.csv'

       with open(out_filename, 'w') as f:
           f.write('# ' + header + '\n')
           for row in S:
               line = ",".join(map(str, row))
               f.write(line + '\n')

   # Cleanup
   # Move files to aux?
   # cmd2 = 'mv result.txt *.exe E.dat out.dat SASA* '+os.path.splitext(arguments.f.name)[0]+'*.txt ./aux'
   # We are not generating .exe or result.txt.
   # We generated E.dat.
   # We generated Output_... .dat and Fig_... .jpg.

   print(u"\U0001F63A", "### Finished ###", u"\U0001F63A")

if __name__ == "__main__":
    main()
