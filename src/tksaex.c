// ##########################################################
// #
// #   		Vinícius Contessoto - 12/2014
// #
// ##########################################################


#include <iostream>
#include <fstream>  //iofstream
#include <sstream>
#include <cmath>
#include <math.h>     
#include <stdlib.h>
#include <stdio.h>  //getchar
#include <vector>
#include "legendre.h"

using namespace std;
using namespace Legendre ;

#define PI 3.14159265359
#define ep 4.0
#define es 78.5
#define K  1.3806488
#define R 8.314
#define eo  8.8541878176
#define e 1.602
#define x 1.0


int main (int argc, char *argv[]) {
  
	int i,j,k,n,a,kk,jj,b,m;
	double Eij[501][501],Gqq[501],X[501],Wqq,Gn,Gu,vi,A,B,AB,Q2,termo1,termo2,termo3,G1,G2,GA,GB,GC,PH,T;
	long int states;
	long double Zn,Zu;
	ofstream out_data("out.dat");
k=0;
PH=atof(argv[1]);
T=atof(argv[2]);
Q2 = (6.02*e*e*pow(10,7))/(4.0*PI*eo);
ifstream f_in_ener("E.dat");
m = 0;
string line;
while (getline(f_in_ener,line)){ m++;}
f_in_ener.clear();
f_in_ener.seekg(0, f_in_ener.beg);

for (int i = 1; i <= m-2; i++){f_in_ener >> Eij[i][0];} // Charges

for (int i = 1; i <=m-2; i++){
  for (int j = 1; j <=m-2; j++){f_in_ener >> Eij[i][j];}}
  
 for (int i = 1; i <= m-2; i++){f_in_ener >> Eij[i][m];}  // PKA
  
n=m-2;

   states=pow(2,n)-1;
   Zn = Zu = 0.0;
   for(jj=1;jj<=states;jj++)
   {
      termo1 = Gn = Gu = vi = 0.0;

	    i=0;
	    kk=jj;
	   	 
	for(j=0;j<=n;j++)
                {
                        X[j]=0.0;
                }


	   while(kk>0)
	    { 
	    X[i]=kk%2; 
	    i++; 
	    kk=kk/2;
	   
	   }
      
   
     for(a=1;a<=n;a++) 
      {
	  termo2 = 0.0;
	  termo1 = (Eij[a][m])*(Eij[a][0] + X[a-1]);
   
   	for(k=1;k<=n;k++)
	  {
	       termo2 = termo2 + 0.5*Eij[a][k]*(Eij[k][0] + X[k-1])*(Eij[a][0] + X[a-1]);
	    
	  }
	
	  Gn = Gn -termo1*(log(10))*R*T + termo2;
	  Gu = Gu -termo1*(log(10))*R*T;
	
	  vi = vi + X[a-1];
	
      }
	
	  termo3 = vi*(log(10))*PH;
	  Zu = Zu + exp( -(Gu)/(R*T) - termo3);
	  Zn = Zn + exp( -(Gn)/(R*T) - termo3);

   }

   G1 = -R*T*log(Zu/Zn);
   G2 = GC = 0.0;   
   
   
   for(b=1;b<=n;b++)
      {
	Wqq = 0.0;  
	
   for(jj=1;jj<=states;jj++)
      {
	  termo1 = Gn = Gu = vi = 0.0;
		for(j=0;j<=n;j++)
      		{
			X[j]=0.0;
		}
	
	    i=0;
            kk=jj;

	   while(kk>0)
      	      {
            X[i]=kk%2;
            i++;
            kk=kk/2;

           }
           
for(a=1;a<=n;a++)
      {
          termo2 = 0.0;
          termo1 = (Eij[a][m])*(Eij[a][0] + X[a-1]);
        for(k=1;k<=n;k++)
          {
               termo2 = termo2 + 0.5*Eij[a][k]*(Eij[k][0] + X[k-1])*(Eij[a][0] + X[a-1]);
          }
                  Gn = Gn -termo1*(log(10))*R*T + termo2;
	  vi = vi + X[a-1];
	}	
	  GC = (exp( -(Gn)/(R*T)  -  vi*(log(10))*PH))/Zn ;
	  termo2 = 0.0;
	
	for(k=1;k<=n;k++)
	  {
	    termo2 = termo2 + 0.5*Eij[b][k]*(Eij[k][0] + X[k-1])*(Eij[b][0] + X[b-1]);
	  } 

	Wqq = Wqq + termo2*GC;
	
      }
   
     Gqq[b] = 0.5*Wqq;
     out_data << Gqq[b]/1000.0 << endl;
     G2 = G2 + Gqq[b];  
      }
}



