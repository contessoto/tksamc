// ##########################################################
// #
// #   		Vinícius Contessoto - 04/2014
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
#define K  1.3806488*pow(10,-23)
#define R 8.314
//#define T  300.0
//#define PH  5.5
#define eo  8.8541878176*pow(10,-12)
#define e 1.602*pow(10,-19)
#define x 0.1

float aleatorio (int a,int b)
{
float r;
 r =  (a - b)*(rand()/(RAND_MAX + 1.0)) + b;
return r;
}

int main (int argc, char *argv[]) {
	
	//char* Args[];
	int i,j,k,n,a,kk,jj,kkk,b,totalp,passos,AUX_MON,Mudadas,passos_des,Pequilib,descorrela,sinalProt,m;
	double Eij[501][501];
	double G1,G2,ionMedio[500],ionQuadrado[500],Q2,count,ionDesvio[500],ionM,ionQ,H,T,EtotalQuadrado[500],Edesvio[500],dg,PH;
        double Etotal[500],cargaN[500],parte2,cont,ET1,DeltaE,Enew,Eold,cargaMon[500], SomaEnergia, imonom[200],teste,dgquadrado,convert;
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

convert = 0.0083145*T;

srand(time(NULL));


for (a=1;a<=n;a++)
      {
       Etotal[a] = 0.0;
       EtotalQuadrado[a] = 0.0;
       Edesvio[a] = 0.0;
       ionDesvio[a] = 0.0;
       ionMedio[a] = 0.0;
       ionQuadrado[a] = 0.0;
       cargaMon[a] = 0;
      }

totalp = 100000;
descorrela = 200;
Pequilib =1000;

for (passos=1;passos<=totalp;passos++)

  {
for (passos_des=1;passos_des<=descorrela;passos_des++)

  {
  Enew = 0.0;
  Eold = 0.0;
  DeltaE = 0.0;
  parte2 = 0.0;
  Mudadas = 1;//int(0.05*n);

  for (a=1;a<=n;a++)
  { 
    imonom[a] = 0.0;         
    cargaN[a] = cargaMon[a]; 
  } 

  for (k=1;k<=Mudadas;k++)
  {
    REFAZ: imonom[k] = int(aleatorio(n+1,1));
    AUX_MON = imonom[k];

    for (j=1;j<k;j++)
      {
        if (imonom[j] == imonom[k])
        {
          goto REFAZ;
        }
      }    

    for (a=1;a<=n;a++)
      {
        if (a!=AUX_MON)
          {
            Eold = Eold + 0.0005*Eij[AUX_MON][a]*cargaMon[AUX_MON]*cargaMon[a];
          }
       }
 

    if (Eij[AUX_MON][0] == 0)
      { 
        if (cargaMon[AUX_MON] == 0)
          { 
            cargaN[AUX_MON] = cargaMon[AUX_MON] + 1;
            parte2 = parte2 + (PH - Eij[AUX_MON][m]);
          }
        else
          {
            cargaN[AUX_MON] = cargaMon[AUX_MON] - 1;
            parte2 = parte2 - (PH - Eij[AUX_MON][m]);
          }
      }
    else
      {
        if (cargaMon[AUX_MON] == 0) 
          {
            cargaN[AUX_MON] = cargaMon[AUX_MON] - 1;
            parte2 = parte2 - (PH - Eij[AUX_MON][m]);
          }
        else
          {
            cargaN[AUX_MON] = cargaMon[AUX_MON] + 1;
            parte2 = parte2 + (PH - Eij[AUX_MON][m]);
          }
      }
    }

 for (k=1;k<=Mudadas;k++)
 {  
    AUX_MON = int(imonom[k]);
    for (a=1;a<=n;a++)
      {    
        if (a!=AUX_MON)
        {
          Enew = Enew + 0.0005*Eij[AUX_MON][a]*cargaN[AUX_MON]*cargaN[a];
        }
      }
  }

   
     DeltaE = (Enew - Eold);
        
     DeltaE = DeltaE + parte2*log(10);

     if (DeltaE <= 0.0)
       {
         for (k=1;k<=Mudadas;k++)
         {
           AUX_MON = int(imonom[k]);
           cargaMon[AUX_MON] = cargaN[AUX_MON];
         }
         cont=cont+1;
       }
     else if (exp(-DeltaE) > aleatorio(1,0))
       {
         for (k=1;k<=Mudadas;k++)
         {
           AUX_MON = int(imonom[k]);
           cargaMon[AUX_MON] = cargaN[AUX_MON];
         }
         cont=cont+1;
       }

   }

     if (passos >= Pequilib){

     for (a=1;a<=n;a++)
     {
      ET1 = 0.0;
      for (k=1;k<=n;k++)
      {
      
       ET1 = ET1 + 0.0005*Eij[a][k]*cargaMon[a]*cargaMon[k];
       
      }
      ionMedio[a] = ionMedio[a] + cargaMon[a];
      ionQuadrado[a] = ionQuadrado[a] + pow(cargaMon[a],2.0);
      Etotal[a] = Etotal[a] + ET1;
      EtotalQuadrado[a] = EtotalQuadrado[a] + pow(ET1,2.0);
     }
     }
  
   } 
    
   
     SomaEnergia = 0.0;

      for (k=1;k<=n;k++)
        {
         dgquadrado = EtotalQuadrado[k]/(totalp-Pequilib);
         dg = Etotal[k]/(totalp-Pequilib);
         Edesvio[k] = sqrt(dgquadrado - pow(dg,2.0));

         ionQ = ionQuadrado[k]/(totalp-Pequilib);
         ionM = ionMedio[k]/(totalp-Pequilib);
         ionDesvio[k] = sqrt(ionQ - pow(ionM,2.0));

         out_data << convert*dg/5.0 << endl;
	 SomaEnergia = SomaEnergia + Etotal[k]/(totalp-Pequilib);
        }  
  
}


