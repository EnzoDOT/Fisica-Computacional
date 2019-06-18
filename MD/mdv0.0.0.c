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
  float *posicion,dL,*Vlj,*Flj,*r,*r2,*fuerzas;
  int *semilla;

  n=cbrt(N);
  float dL=L/n;  
  int gf=(int) 2.5/(5000.0)   //Grilla para interpolar potencial, fuerzas y distancias
  Vlj=(float*) malloc(gf*sizeof(float));
  Flj=(float*) malloc(gf*sizeof(float));
  r=(float*) malloc(gf*sizeof(float));
  r2=(float*) malloc(gf*sizeof(float));
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;
  posicion=(float*) malloc(3*N*sizeof(float));
  fuerzas=(float*) malloc(3*N*sizeof(float)); 
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

void set_box(float *posicion; int n; float L)
{ int i=0;
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
    *(r+i)=i*dL+r0;
    *(r2+i)=(i*dL)*(i*dL)+r0*r0;
    *(Vlj+i)=4.0*(pow(*(r+i), -12.0)-pow(*(r+i), -6.0));
    *(Flj+i)=24.0*(2.0*pow(*(r+i), -13.0)-pow(*(r+i), -7.0));
   }
}

void interaccion(float *posicion, *float *fuerzas)
{  

   float n2; 

   for (i=1; i<N;i=i++)
   {
    for (j=0; j<i;j=j++)    
    {/* Acà hay que implementar el DeltaX para que evalùe con las 
       partìculas fantasmas o reales. Si la distancia entre partìculas reales
       es mayor a L/2, uso la fantasma*/
     Dx=deltax(posicion+3*i,posicion+3*j);
     Dy=deltax(posicion+3*i+1,posicion+3*j+1);
     Dz=deltax(posicion+3*i+2,posicion+3*j+2);
     n2=Dx*Dx*+Dy*Dy+Dz*Dz;
//     n2=norma(posicion,i,j);
     {//if(r^2<rc^2):...
      if(n2<rc*rc)
      {
       Fx=Dx*(*(Flj+indice))/(*(r+indice));
       Fy=Dy*(*(Flj+indice))/(*(r+indice));
       Fz=Dz*(*(Flj+indice))/(*(r+indice));
       *(fuerzas+3*i)+=Fx;
       *(fuerzas+3*i+1)+=Fy;
       *(fuerzas+3*i+2)+=Fz;
       *(fuerzas+3*j)=*(fuerzas+3*j)-Fx;
       *(fuerzas+3*j+1)=*(fuerzas+3*j+1)-Fy;
       *(fuerzas+3*j+2)=*(fuerzas+3*j+2)-Fz;
      }
     }
    }  
   }
}

float deltax(float *posicion1, float *posicion2)
{  
   float Dx; 
   Dx=*(posicion1)-*(posicion2);
   if(Dx>L/2.0) Dx=Dx-L;
   if(Dx<-L/2.0) Dx=Dx+L;

  return (float)Dx;
}




