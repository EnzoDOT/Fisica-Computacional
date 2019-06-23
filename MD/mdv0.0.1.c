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

double irandom(int *semilla);
double gaussiana (double mu, double sigma, int *semilla);
void set_box(double *posicion, int N, double L)
void set_v(double *v, int N, double T, int *semilla)
void set_tablas(double *Vlj,double *Flj, double *r, double *r2 , int * gf, double *dLt)
void interaccion(double *posicion, *double *fuerzas, double *Vlj, double *Flj, double *r, double *r2, double *dLt )
double deltax(double *posicion1, double *posicion2)

void imprimir (int *red, int dim);

int main (int argc, char *argv[])
{ 
  FILE *fp, *fp2;
  int i,N;
  double *posicion,dL,*Vlj,*Flj,*r,*r2,*fuerzas,*n;
  int *semilla;

//  n=(double*) cbrt(N); //Para generar ci en sc
//  double dL=L/n;  
  dLt=(double*) 2.5/(5000.0)   //Grilla para interpolar potencial, fuerzas y distancias
  int gf=(int*) 2.5/dLt
  Vlj=(double*) malloc(gf*sizeof(double));
  Flj=(double*) malloc(gf*sizeof(double));
  r=(double*) malloc(gf*sizeof(double));
  r2=(double*) malloc(gf*sizeof(double));
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;
  posicion=(double*) malloc(3*N*sizeof(double));
  fuerzas=(double*) malloc(3*N*sizeof(double)); 
  fp=fopen("md1energia.dat","w");
  fp2=fopen("md.dat","w");

   N=100;
   if(argc==2)
    {
     sscanf(argv[1],"%d",&N);
     }
  n=(double*) cbrt(N);

  free(semilla);
  fclose(fp);
  return 0;
}

double irandom(int *semilla)
{ 
  int k;
  double x;
  k=(*semilla)/Q;       
  *semilla=A*(*semilla-k*Q)-k*R;
  if(*semilla<0) *semilla=(*semilla)+M;
  x=(*semilla)*(1.0/M);
  return (double)x;
}

double gaussiana(mu, sigma, *semilla)
{ 
  int n=10;
  double z=0.0;
  double x;
  for(int i=0; i<n; i++)
  {
   z+=irandom(semilla);
  }
  z=sqrt(12.0*n)*(z/n - 0.5);
  return (double)z*sigma+mu;
}

void set_box(double *posicion; int *n; double L)
{ int i=0;
  double dL=L/n;
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

void set_v(double *v; int N; double T,int *semilla)
{ double sigma=sqrt(N);
  for (i=0; i<3*N;i++)
  {
  *(v+i)=gaussiana(0.0,sigma,semilla);
  }
  double vcm[3]={0,0,0};
  for(i=0;i<N;i++)
  {
   for(k=0;k<3;k++)
   {
    vcm[k]+=*(v+3*i+k)/N;
   }
  }
  for(i=0;i<N;i++)
  {
   for(k=0;k<3;k++)
   {
    *(v+3*i+k)-=vcm[k];
   }
  }
}

void set_tablas(double *Vlj,double *Flj, double *r, double *r2 , int * gf, double *dLt)
{  
   int i;
   for (i=1; i<gf+1;i=i++)
   {
    *(r+i)=i*dLt;
    *(r2+i)=(i*dLt)*(i*dLt);
    *(Vlj+i)=4.0*(pow(*(r+i), -12.0)-pow(*(r+i), -6.0))\
       - 4.0*(pow(2.5, -12.0)-pow(2.5, -6.0));
    *(Flj+i)=24.0*(2.0*pow(*(r+i), -13.0)-pow(*(r+i), -7.0));
   }
}

void interaccion(double *posicion, *double *fuerzas, double *Vlj, double *Flj, double *r, double *r2, double *dLt )
{  int i,j,indice;
   double n2,Fx,Fy,Fz,Dx,Dy,Dz; 
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
       indice=(int)((n2-*(r+1))/dLt);
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

double deltax(double *posicion1, double *posicion2)
{  
   double Dx; 
   Dx=*(posicion1)-*(posicion2);
   if(Dx>L/2.0) Dx=Dx-L;
   if(Dx<-L/2.0) Dx=Dx+L;

  return (double)Dx;
}

void verlet(double *posicion, double *v, double *fuerzas, double *dt, \
            double *Vlj,double *Flj, double *r, double *r2, double *dLt)
{ 
   double dt; 
   for (i=0; i<N;i=i++)//Primer paso de Verlet
   {   
    *(posicion+3*i)+=*(v+3*i)*dt+*(fuerzas+3*i)*dt*dt*0.5;
    *(posicion+3*i+1)+=*(v+3*i+1)*dt+0.5**(fuerzas+3*i+1)*dt*dt;
    *(posicion+3*i+2)+=*(v+3*i+2)*dt+0.5**(fuerzas+3*i+2)*dt*dt;
    *(v+3*i)+=*(fuerzas+3*i)*dt*0.5;
    *(v+3*i+1)+=*(fuerzas+3*i+1)*dt*0.5;
    *(v+3*i+2)+=*(fuerzas+3*i+2)*dt*0.5;
   }

   interaccion(double *posicion, double *fuerzas, double *Vlj, \
                     double *Flj, double *r, double *r2, double *dLt) //actualizo \
                                                                  las fuerzas
   for (i=0; i<N;i=i++)//Segundo paso de Verlet
   {   
    *(v+3*i)+=*(fuerzas+3*i)*dt*dt*0.5;
    *(v+3*i+1)+=*(fuerzas+3*i+1)*dt*dt*0.5;
    *(v+3*i+2)+=*(fuerzas+3*i+2)*dt*dt*0.5;
   }

}




