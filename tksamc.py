#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
TKSA - Electrostatic Free Energy calculation for each ionizable residue
"""

__author__ = "Vinícius G. Contessoto"
__date__ = "21/12/2016"

################################################################
#
# Version 1.0
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

parser = argparse.ArgumentParser(description='Charge-charge energy calculation in python')
parser.add_argument('-ph', action='store', default=7.0, dest='arg_pH', help='pH value')              
parser.add_argument('-T', action='store', default=300.0, dest='arg_T',  help='Temperature value')           
parser.add_argument('-f', metavar='input-file-PDB',help='insert a PDB file',type=argparse.FileType('rt'))
parser.add_argument('-e', action='store',choices=['TK'], default="TK",dest='arg_e',type=str,help='Electrostatic energy calculation method')
parser.add_argument('-s', action='store',choices=['EX','MC'], default="MC",dest='arg_s',type=str,help='Statistical method to protonation state amostration - EX = Exact; MC = Monte Carlo;')
parser.add_argument('-plot', action='store',choices=['yes','no'], default="yes",dest='arg_plot',type=str,help='Save Plot figure file - EPS')

try:
    arguments = parser.parse_args()
    print '################################################'
    print u"\U0001F63A", "### TKSA started ###", u"\U0001F63A"
    print 'Input file:', arguments.f.name
    print 'pH  =', arguments.arg_pH
    print 'T   =', arguments.arg_T
    print 'Elec. Energy Calc Method =', arguments.arg_e
    print 'Statistical Method =', arguments.arg_s
    print 'Plot =', arguments.arg_plot
    
except IOError, msg:
    parser.error(str(msg))                    
                  

All_residues = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
Area_residues = [113,140,151,183,218,85,194,182,211,180,204,158,143,189,241,122,146,160,259,229] # colocar a referencia
Charged_residues = ['ARG','LYS','N_TER','HIS','GLU','ASP','C_TER']
Charge_values = [0,0,0,0,-1,-1,-1]
Charged_atoms = ['NH2','NZ','NE2','OE2','OD2']
PKA = [12.0,10.6,7.7,6.3,4.5,4.0,3.6]
e = (4.0,78.5,1.0,6.02,1.602,8.8541878176,1.3806488,8.314) #(ep,es,x,mol,e,eo,K,R) reduzidos

##################################################################################################
# Kirkwood Polynomial Function
##################################################################################################
def Kn(n,x):
   
   Kpf = np.sum([np.power(x,s)*np.divide(np.power(2.0,s)*np.math.factorial(n)*np.math.factorial(2*n-s),np.math.factorial(s)*np.math.factorial(2*n)*np.math.factorial(n-s)) for s in range(n)])
   return Kpf
##################################################################################################

def main():
   global Q,E,S,Pk,e,T,pH,total_charged_residues,G,G0,indiv_data,Gqq,Q0

   file_pdb = arguments.f     
   pH = np.float(arguments.arg_pH)
   T = np.float(arguments.arg_T)
   
   ##################################################################################################
   # runs the standalone version of ©Surfrace
   ##################################################################################################

      
   print 'Running SASA - ©Surfrace'
   cmd1 = 'echo 1' + arguments.f.name + ' 1.4 1| ./surfrace5_0_linux_64bit > SASA_'+os.path.splitext(arguments.f.name)[0]+'_all.trash' ## Roda o programa para a SASA
   os.system(cmd1)
   try: 
     file_sasa = open(os.path.splitext(arguments.f.name)[0] + '_residue.txt', 'r') ## Abre o arquivo que vem do programa acima
   except (IOError) as errno:
     print ('I/O error - ** Check the files of SASA calculation - something went wrong **. %s' % errno)
     sys.exit()
   
   SASA_data=[]
   for line2 in file_sasa:
     list2 = line2.split()
     Area_norm = np.float(list2[2])/np.float(Area_residues[All_residues.index(list2[1])])
     if Area_norm >= 1.0:
       print "Warning - ** SASA greater than 1.0 **",list2[1],list2[0],list2[2],np.float(Area_residues[All_residues.index(list2[1])]),Area_norm
       print "Automatically changed to 0.75"
       Area_norm = 0.750000000001
     SASA_data.append([list2[1],list2[2],Area_norm])
   indiv_data=[]
   S=[]
   SAij=[]
   total_atoms=[]
   total_residues=[]
   total_charged_residues=[]
   for line in file_pdb: ## Reading file.pdb
     lista = line.split()
     id = lista[0]
     if id == 'ATOM':
       atom_index = np.int(lista[1]) 
       atom_type = lista[2]
       residue_type = lista[3]
       chain = lista[4]
       residue_index = np.int(lista[5])
       total_atoms.append([atom_index])
       if atom_type == 'CA' and chain == 'A':
	 total_residues.append([residue_index])
       if atom_index == 1 and atom_type == 'N' and chain == 'A' and residue_index == 1 and not residue_type in Charged_residues: ## Select the charged residues
	 total_charged_residues.append([atom_index])
	 S.append(['N_T',residue_index,np.size(total_charged_residues),lista[1],lista[2],np.float(lista[6]),np.float(lista[7]),np.float(lista[8]),PKA[Charged_residues.index('N_TER')],SASA_data[np.size(total_residues)-1][2],Charge_values[Charged_residues.index('N_TER')]])
       if residue_type in Charged_residues and atom_type in Charged_atoms: ## Seleciona os resíduos carregados
	 total_charged_residues.append([atom_index])
         S.append([lista[3],residue_index,np.size(total_charged_residues),lista[1],lista[2],np.float(lista[6]),np.float(lista[7]),np.float(lista[8]),PKA[Charged_residues.index(residue_type)],SASA_data[np.size(total_residues)-1][2],Charge_values[Charged_residues.index(residue_type)]])
       if atom_type == 'OXT' and chain == 'A'  and not residue_type in Charged_residues:
	 total_charged_residues.append([atom_index])
	 S.append(['C_T',residue_index,np.size(total_charged_residues),lista[1],lista[2],np.float(lista[6]),np.float(lista[7]),np.float(lista[8]),PKA[Charged_residues.index('C_TER')],SASA_data[np.size(total_residues)-1][2],Charge_values[Charged_residues.index('C_TER')]])

   print "There are: %d Charged_residues" % np.size(total_charged_residues)
   Restype=np.asarray([i[0] for i in S])
   X=np.asarray([i[5] for i in S])
   Y=np.asarray([i[6] for i in S])
   Z=np.asarray([i[7] for i in S])
   Pk=np.asarray([i[8] for i in S])
   SA=np.asarray([i[9] for i in S])
   Q=np.asarray([i[10] for i in S])
   Restype=np.char.replace(np.char.replace(np.char.replace(np.char.replace(np.char.replace(Restype, 'HIS','H'), 'ASP','D'), 'ARG','R'), 'GLU','E'), 'LYS','K')
   X = X - np.mean(X)
   Y = Y - np.mean(Y)
   Z = Z - np.mean(Z)
   XYZ = zip(X,Y,Z)
   Origin = np.zeros(np.shape(XYZ))
   dist = distance.cdist(XYZ, XYZ, 'euclidean')
   if arguments.arg_e == 'TK':
    dist_origin = distance.cdist(XYZ, Origin, 'euclidean')
    angle = distance.cdist(XYZ, XYZ, 'cosine')
    raio = (np.max(dist)*0.5 + 3.4+2.0, np.max(dist)*0.5 + 2.0+2.0)
    np.seterr(invalid='ignore')
    np.seterr(divide='ignore')
    theta = np.arccos(1-angle)
    NormA = np.matrix([LA.norm(v) for v in np.array(XYZ)])
    rirj = np.array(np.dot(np.transpose(NormA),NormA))
    A = np.divide(raio[1],e[0]*dist)
    B = (np.nansum(np.array([((e[1]-e[0])/(e[1]-(n*e[0])/(n+1)))*(np.power((rirj/(raio[1]*raio[1])),n))*(eval_legendre(n, np.cos(theta))) for n in range(0,50)]),axis=0))/(e[0]) 
    C = (np.divide(e[2],1+e[2]) + np.power(e[2],2)*np.sum(np.array([np.divide(np.divide(2*n+1,2*n-1)*np.divide(e[1],(n+1)*e[1]+n*e[0])*(np.power((rirj/(raio[0]*raio[0])),n))*(eval_legendre(n, np.cos(theta))),np.divide(Kn(n+1,e[2]),Kn(n-1,e[2])) + np.divide(n*(e[1]-e[0]),(n+1)*e[1]+n*e[0])*np.divide(np.power(e[2],2),4.0*np.power(n,2)-1)*np.power(np.divide(raio[1],raio[0]),2*n+1)) for n in range(1,50)]),axis=0))/(e[1])
    Qe = np.divide(e[3]*e[4]*e[4]*np.power(10,7),4*np.pi*e[5])
    SAij = distance.cdist(zip(SA), zip(SA), lambda u,v: (u+v)*0.5)
    E = Qe*(np.divide(A-B,2*raio[1])-np.divide(C,2*raio[0]))*(1-SAij)
 
    if np.sum(np.where(E<0)) > 0:
      print '###############################################################'
      print "There are: %d negatives TK energy values - Please check the radius of TK method!" % np.int(np.sum(np.where(E<0)))
      print "Sugestion - Increase b radius"
      print "Current radius ratio b/a=", np.divide(raio[1],raio[0])
      print '###############################################################'
    E[np.isinf(E)]= 0
    E[np.isnan(E)]= 0   
    E_out=np.vstack([np.vstack([Q,E]),Pk])
    np.savetxt('E.dat',E_out)
    
   if arguments.arg_s == 'EX': 
      
      print u"\U0001F63A", "### TK - Exact ###", u"\U0001F63A"
      start = time.time()        
      p = subprocess.Popen([r"c++","./src/tksaex.c",'-lm','-O3','-o','tksaex.exe'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
      p.communicate()
      p = subprocess.Popen(["./tksaex.exe",np.str(pH),np.str(T)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
      i=0
      j=1
      while p.poll() is None:
          sys.stdout.write('\r')
          sys.stdout.write("TKSA EX is running, please wait - [%20s%-20s]" % ('='*i,'='*i))
          sys.stdout.write(u"\U0001F63A")
          sys.stdout.flush()
          if i>19:
           j=j+1
          if j%2 == 0:
           i=i-1
          if j%2 == 1:
           i=i+1
          if i == 0:
            j=1
          sys.stdout.flush()
          time.sleep(0.1)
          
      output,err = p.communicate()
      print output
      print err
      end = time.time()
      elapsed = end - start
      print "Ran in %f sec" % elapsed
      
   if arguments.arg_s == 'MC': 
      
      print u"\U0001F63A", "### TKSA - MC ###", u"\U0001F63A"
      start = time.time()        
      p = subprocess.Popen([r"c++","./src/tksamc.c",'-lm','-O3','-o','tksamc.exe'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
      p.communicate()
      p = subprocess.Popen(["./tksamc.exe",np.str(pH),np.str(T)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
      i=0
      j=1
      while p.poll() is None:
          sys.stdout.write('\r')
          sys.stdout.write("TKSA MC is running, please wait - [%20s%-20s]" % ('='*i,'='*i))
          sys.stdout.write(u"\U0001F63A")
          sys.stdout.flush()
          if i>19:
           j=j+1
          if j%2 == 0:
           i=i-1
          if j%2 == 1:
           i=i+1
          if i == 0:
            j=1
          sys.stdout.flush()
          time.sleep(0.1)
          
      output,err = p.communicate()
      print output
      print err
      end = time.time()
      elapsed = end - start
      print "Ran in %f sec" % elapsed
      
   if arguments.arg_plot == 'yes' and arguments.arg_s =='EX': 
      try: 
        file_plot = open("out.dat", 'r')
      except (IOError) as errno:
        print ('I/O error - ** Output file with issues  - out.dat **. %s' % errno)
        sys.exit()
   
      plot_data=[]
      for line3 in file_plot: # Plotting
        list3 = line3.split()
        plot_data.append(list3)
      
      Restype=np.char.replace(np.char.replace(["%s%02d" % t for t in zip(Restype,np.asarray([i[1] for i in S]))],'C_T'+np.str(S[-1][1]),'CTR'),'N_T0'+np.str(S[0][1]),'NTR')
      S=np.hstack((S,plot_data))
      
      plot_data=list(map(float, np.asarray(plot_data).flatten()))
      print "Total dG Energy: ",np.sum(np.asarray(plot_data))
      x_pos = np.arange(len(total_charged_residues))
      fig = plt.figure()
      ax = fig.add_subplot(111)
      width=1.0
      colors = []
      for position, value in enumerate(plot_data):
        if value > 0 and SA[position] > 0.5:
           colors.append('r')
        else:
           colors.append('b')
      ax.bar(x_pos, plot_data,width=width,color=colors,linewidth=2)
      ax.tick_params('both', length=5, width=2, which='major',labelsize=13)   
      plt.setp(ax.spines.values(), linewidth=2)
      plt.xticks(x_pos+width/2.0,Restype,rotation=90,fontsize=15)
      plt.xlim([0,np.size(x_pos)])
      plt.ylabel(r'$\Delta G_{qq}$(kJ/mol)',fontsize=20)
      fig.savefig('Fig_EX_'+ os.path.splitext(arguments.f.name)[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.jpg', dpi = 300)
   
      header='1-Name	2-Residue-index	3-Position	4-Atom	5-Atom-type	6-X	7-Y	8-Z	9-PKA	10-SASA	11-Charge	12-dG_Energy 13-Total_dG= '+str(np.sum(np.asarray(plot_data)))+''
      np.savetxt('Output_EX_'+os.path.splitext(arguments.f.name)[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.dat',S,fmt='%s', delimiter="	",header=str(header))
   
   if arguments.arg_plot == 'yes' and(arguments.arg_s =='MC'): 
      try: 
        file_plot = open("out.dat", 'r')
      except (IOError) as errno:
        print ('I/O error - ** Output file with issues - out.dat **. %s' % errno)
        sys.exit()
   
      plot_data=[]
      for line3 in file_plot: ## Plotting
        list3 = line3.split()
        plot_data.append(list3)
      
      Restype=np.char.replace(np.char.replace(["%s%02d" % t for t in zip(Restype,np.asarray([i[1] for i in S]))],'C_T'+np.str(S[-1][1]),'CTR'),'N_T0'+np.str(S[0][1]),'NTR')
      S=np.hstack((S,plot_data))
      
      plot_data=list(map(float, np.asarray(plot_data).flatten()))
      print "Total dG Energy: ",np.sum(np.asarray(plot_data))
      x_pos = np.arange(len(total_charged_residues))
      fig = plt.figure()
      ax = fig.add_subplot(111)
      width=1.0
      colors = []
      for position, value in enumerate(plot_data):
        if value > 0 and SA[position] > 0.5:
           colors.append('r')
        else:
           colors.append('b')
      ax.bar(x_pos, plot_data,width=width,color=colors,linewidth=2)
      ax.tick_params('both', length=5, width=2, which='major',labelsize=13)   
      plt.setp(ax.spines.values(), linewidth=2)
      plt.xticks(x_pos+width/2.0,Restype,rotation=90,fontsize=15)
      plt.xlim([0,np.size(x_pos)])
      plt.ylabel(r'$\Delta G_{qq}$(kJ/mol)',fontsize=20)
      plt.show()
      fig.savefig('Fig_MC_'+ os.path.splitext(arguments.f.name)[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.jpg', dpi = 300)
   
      header='1-Name	2-Residue-index	3-Position	4-Atom	5-Atom-type	6-X	7-Y	8-Z	9-PKA	10-SASA	11-Charge	12-dG_Energy 13-Total_dG= '+str(np.sum(np.asarray(plot_data)))+''
      np.savetxt('Output_MC_'+os.path.splitext(arguments.f.name)[0]+'_pH_'+str(pH)+'_T_'+str(T)+'.dat',S,fmt='%s', delimiter="	",header=str(header))
   cmd2 = 'mv result.txt *.exe E.dat out.dat SASA* '+os.path.splitext(arguments.f.name)[0]+'*.txt ./aux'
   os.system(cmd2)
   print u"\U0001F63A", "### Finished ###", u"\U0001F63A"
      
if __name__ == "__main__": main()
