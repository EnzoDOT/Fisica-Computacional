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
  float p,x,delta,xp,fx,fxp,pa,xm2,xm,*c,*cx,cc;
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
    pa=0.0;   
    x=0.0;
    c=(float*) malloc(N/2*sizeof(float));
    cx=(float*) malloc(N*sizeof(float));
    xm2=0.0;
    xm=0.0;

   for(delta=0.1;delta<10.0;delta=delta+0.1)
    {
    for(i=0;i<N;i++)
    {
     xp=x;    
     x=irandom(semilla);   
     x=xp+(2.0*x-1.0)*delta;
     fx=exp((-x*x)/2.0);
     fxp=exp((-xp*xp)/2.0);
     if(fx<=fxp)
     {
      p=exp((xp*xp-x*x)/2.0);
      if(p>1.0) p=1.0;
     } 
     else
     {
      p=1.0;
     }
      if(irandom(semilla)<p)  
      {     
      pa=pa+1.0;//p;
//     fprintf(fp,"%f %f \n", x, pa);   
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
    pa=pa/N;
    fprintf(fp,"%f %f \n", delta, pa);
    
    xm=xm/N;
    xm2=xm2/N;
    
    for(k=0;k<N/2;k++)
    {
    *(c+k)=0;
    for(i=0;i<N-k;i++)
     {
      *(c+k)=(*(cx+i+k))*(*(cx+i));
     }
     *(c+k)=*(c+k)/(N-k);
     cc=(*(c+k)-xm2)/(xm2-xm*xm);
    } 
    
  }
  

  free(cx);  
  free(c);
  free(semilla);
  fclose(fp);
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


