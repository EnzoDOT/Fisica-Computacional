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
void set_box(float *posicion, int N, float L)
void set_v(float *v, int N, float T, int *semilla)
void imprimir (int *red, int dim);

int main (int argc, char *argv[])
{ 
  FILE *fp, *fp2;
  int i,N;
  float *posicion,dL,*Vlj,*Flj,*r,*r2;
  int *semilla;

  n=cbrt(N);
  float dL=L/n;  
  int gf=(int) 2.5/(100.0*dL)   //Grilla para interpolar potencial, fuerzas y distancias
  Vlj=(float*) malloc(gf*sizeof(float));
  Flj=(float*) malloc(gf*sizeof(float));
  r=(float*) malloc(gf*sizeof(float));
  r2=(float*) malloc(gf*sizeof(float));
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;
  posicion=(float*) malloc(3*N*sizeof(float)); 
  fp=fopen("md1energia.dat","w");
  fp2=fopen("md.dat","w");

   N=100;
   if(argc==2)
    {
     sscanf(argv[1],"%d",&N);
     }

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

float gaussiana(mu, sigma, *semilla)
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

void set_box(float *posicion; int N; float L)
{ int n=cbrt(N); i=0;
  float dL=L/n;
  for(int x=0; x<n; x++)
  {for (int y=0; y<n; y++)
    {for (int z=0; z<n; z++)
     {*(posicion+3*i)=dL*(x+0.5);
      *(posicion+3*i+1)=dL*(y+0.5);
      *(posicion+3*i+2)=dL*(z+0.5);
      i++;
     }
    }  
  }
}

void set_v(float *v; int N; float T,int *semilla)
{ float sigma=sqrt(N);
  for (i=0; i<3*N;i++)
  {
  v[i]=gaussiana(0.0,sigma,semilla);
  }
  float vcm[3]={0,0,0};
  for(i=0;i<N;i++)
  {
   for(k=0;k<3;k++)
   {
    vcm[k]+=v[3*i+k]/N;
   }
  }
  for(i=0;i<N;i++)
  {
   for(k=0;k<3;k++)
   {
    v[3*i+k]-=vcm[k];
   }
  }
}

void set_tablas(float *Vlj,float *Flj, float *r, float *r2, int N, float L)
{  
   int n=cbrt(N); 
   float dL=L/n,r; 
   for (i=1; i<gf+1;i=i++)
   {
    *(r+i)=i*dL;
    *(r2+i)=(i*dL)*(i*dL);
    *(Vlj+i)=4.0*(pow(*(r+i), -12.0)-pow(*(r+i), -6.0));
    *(Flj+i)=24.0*(2.0*pow(*(r+i), -13.0)-pow(*(r+i), -7.0));
   }
}

float norma(float *posicion, int i, int j)
{  
   int n=cbrt(N); 
   float dL=L/n,r; 
   Dx=*(posicion+3*i)-*(posicion+3*j);
   Dy=*(posicion+3*i+1)-*(posicion+3*j+1);
   Dz=*(posicion+3*i+2)-*(posicion+3*j+2);
   n2=sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
  return (float)n2;
}




