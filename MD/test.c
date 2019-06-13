#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*Metropolis por Pedro y Enzo. Mayo 2019*/

#define M 2147483647
#define A 16807
#define Q 127773
#define R 2836  
#define S 2


float irandom(int *semilla);
float gaussiana (float mu, float sigma, int *semilla);
int main (int argc, char *argv[])
{ 
  FILE *fp, *fp2, *fp3;
  int i,N,k;
  float p,x,delta,xp,pa,xm2,xm,*c,*cx,sigma;
  int *semilla;
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;


   fp2=fopen("gaussiana.dat","w");


   N=100;
   if(argc==2)
    {
     sscanf(argv[1],"%d",&N);
     }
    c=(float*) malloc(N/2*sizeof(float));
    cx=(float*) malloc(N*sizeof(float));

     sigma=1.0;
     delta=2.9;
     pa=0.0;   
     x=0.0;
     xm2=0.0;
     xm=0.0;
     for(i=0;i<N;i++)
     {
      gaussiana(0.0,sigma,semilla);

      fprintf(fp2,"%f \n",  gaussiana(0.0,sigma,semilla));

     }

  free(cx);  
  free(c);
  free(semilla);
  fclose(fp2);
  return 0;
}

float irandom(int *semilla)
{ 
  int k;
  float x;
  k=(*semilla)/Q;       
  *semilla=A*(*semilla-k*Q)-k*R;
  if(*semilla<0) *semilla=(*semilla)+M;
  x=(*semilla)*(1.0/M);
  return (float)x;
}

float gaussiana(float mu, float sigma, int *semilla)
{ 
  int n=10;
  float z=0.0;
  float x;
  for(int i=0; i<n; i++)
  {
   z+=irandom(semilla);
  }
  z=sqrt(12.0*n)*(z/n - 0.5);
  return (float)z*sigma+mu;
}


