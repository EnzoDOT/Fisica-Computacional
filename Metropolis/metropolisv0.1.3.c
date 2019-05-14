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

int main (int argc, char *argv[])
{ 
  FILE *fp, *fp2, *fp3;
  int i,N,k;
  float p,x,delta,xp,pa,xm2,xm,*c,*cx;
  int *semilla;
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;

   fp=fopen("metropolis1a.dat","w");
   fp2=fopen("1amuestra.dat","w");
   fp3=fopen("1acoefc.dat","w");   

   N=100;
   if(argc==2)
    {
     sscanf(argv[1],"%d",&N);
     }
    c=(float*) malloc(N/2*sizeof(float));
    cx=(float*) malloc(N*sizeof(float));

   for(delta=0.1;delta<5.0;delta=delta+0.1)
    {
     pa=0.0;   
     x=0.0;
     xm2=0.0;
     xm=0.0;
     for(i=0;i<N;i++)
     {
      xp=x; 
      x=irandom(semilla);
      x=xp+(2.0*x-1.0)*delta;
      p=exp((xp*xp-x*x)/2.0);
      if(irandom(semilla)<=p)  
      {     
       pa=pa+1.0;
      }
      else
      {
       x=xp;
      }
      fprintf(fp2,"%f \n", x);
      *(cx+i)=x;
      xm2=xm2+x*x;
      xm=xm+x;
    }
    pa=pa/(float)N;
    xm=xm/(float)N;
    xm2=xm2/(float)N;
    fprintf(fp,"%f %f \n", delta, pa);

    for(k=0;k<N/2;k++)
    {
    *(c+k)=0;
    for(i=0;i<N-k;i++)
     {
      *(c+k)=*(c+k)+(*(cx+i+k))*(*(cx+i));
     }
     *(c+k)=*(c+k)/(N-k);
     *(c+k)=(*(c+k)-xm*xm)/(xm2-xm*xm);
     fprintf(fp3,"%d %f \n", k, *(c+k));     
    } 
     fprintf(fp3,"\n ");       
  }  
  free(cx);  
  free(c);
  free(semilla);
  fclose(fp);
  fclose(fp2);
  fclose(fp3);
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


