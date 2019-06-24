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

void imprimir(double *posicion, int N, double L,double t);
double irandom(int *semilla);
double gaussiana (double mu, double sigma, int *semilla);
void set_box(double *posicion, double n, double L);
void set_v(double *v, int N, double T, int *semilla);
void set_tablas(double *Vlj,double *Flj, double *r, double *r2 , int gf, double dLt);
void interaccion(double *posicion, double *fuerzas, double *Vlj, double *Flj, double *r, double *r2, double dLt, double rc2, int N, double L);
double deltax(double *posicion, int i ,int j ,int p, double L);
void verlet(double *posicion, double *v, double *fuerzas, double dt, double *Vlj,double *Flj, double *r, double *r2, double dLt, int N, double L, double rc2);
double hamiltoniano(double *posicion, double *v, double *fuerzas, double dt, double *Vlj,double *Flj, double *r, double *r2, double dLt, double rc2, int N, double L);

int main (int argc, char *argv[])
{ 
  FILE *fp, *fp2;
  int i,N,sfreq,tsteps;
  double *posicion,*v,*Vlj,*Flj,*r,*r2,*fuerzas,n,rc2, dt,L,tmax,T,dLt;
  int *semilla,gf;
  fp=fopen("md1energia.dat","w");
  fp2=fopen("md.dat","w");
   L=2.0;
   T= 1.0;
   tmax=10.0;
   N=27;
   gf=250;
   tsteps=100;

   if(argc==7)
    {sscanf(argv[1],"%lf",&L);
     sscanf(argv[2],"%lf",&T);
     sscanf(argv[3],"%d",&N);
     sscanf(argv[4],"%lf",&tmax);
     sscanf(argv[5],"%d",&tsteps);
     sscanf(argv[6],"%d",&gf);
     }
  printf("Hola");
  dt=tmax/tsteps;
  sfreq=tsteps/10;
  n=cbrt(N);
  dLt=2.5/((double) gf); /*Grilla para interpolar potencial,
                              fuerzas y distancias, Ntablas=5000.0*/
  rc2=(double) 2.5*2.5;
  gf=(int)2.5/dLt;
  Vlj=(double*) malloc(gf*sizeof(double));
  Flj=(double*) malloc(gf*sizeof(double));
  r=(double*) malloc(gf*sizeof(double));
  r2=(double*) malloc(gf*sizeof(double));
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;
  posicion=(double*) malloc(3*N*sizeof(double));
  fuerzas=(double*) malloc(3*N*sizeof(double)); 
  v=(double*) malloc(3*N*sizeof(double));
//  printf("Hola")  
  set_box(posicion,n,L);
  set_v(v, N, T, semilla);
  set_tablas(Vlj, Flj, r, r2, gf, dLt);
  interaccion(posicion, fuerzas, Vlj, Flj, r, r2, dLt, rc2, N, L);

  for(i=0;i<tsteps;i++)
  {
   verlet(posicion, v, fuerzas, dt, Vlj, Flj,  r, r2, dLt, N, L, rc2);
    // num % 2 computes the remainder when num is divided by 2
   if(i%sfreq==0) imprimir(posicion,N,L,i*dt);
  }
  
  free(posicion);
  free(fuerzas);
  free(v);
  free(semilla);
  free(Vlj);
  free(Flj);
  free(r);
  free(r2);
  fclose(fp);
  fclose(fp2);
  return 0;
}

void imprimir (double *posicion, int N, double L,double t)
{ int i;
  char filename[255]; 
   sprintf(filename,"Posicion_L=%lf_N=%4.2d_t=%lf.dat",L,N,t);
   FILE *fp=fopen(filename,"w");
  for(i=0;i<N;i++)
  {
   fprintf(fp,"%lf %lf %lf \n",*(posicion+3*i),*(posicion+3*i+1),*(posicion+3*i+2));       
  }
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

double gaussiana(double mu, double sigma, int *semilla)
{
  int i;
  int n=10;
  double z=0.0;
  for(i=0;i<n; i++)
  {
   z+=irandom(semilla);
  }
  z=sqrt(12.0*n)*(z/n - 0.5);
  return (double)z*sigma+mu;
}

void set_box(double *posicion, double n, double L)
{ int i,x,y,z;
  i=0;
  double dL=L/n;
  for(x=0; x<(int)n; x++)
  {for (y=0; y<(int)n; y++)
    {for(z=0; z<(int)n; z++)
     {*(posicion+3*i)=dL*(x+0.5);
      *(posicion+3*i+1)=dL*(y+0.5);
      *(posicion+3*i+2)=dL*(z+0.5);
      i++;
     }
    }  
  }
}

void set_v(double *v, int N, double T,int *semilla)
{ double sigma=sqrt(T);
  int i,k;
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

void set_tablas(double *Vlj,double *Flj, double *r, double *r2 , int gf, double dLt)
{  
   int i;
   for(i=1;i<gf+1;i++)
   {
    *(r+i-1)=i*dLt;
    *(r2+i-1)=(i*dLt)*(i*dLt);
    *(Vlj+i-1)=4.0*(pow(*(r+i), -12.0)-pow(*(r+i), -6.0)) - 4.0*(pow(2.5, -12.0)-pow(2.5, -6.0));
    *(Flj+i-1)=24.0*(2.0*pow(*(r+i), -13.0)-pow(*(r+i), -7.0));
   }
}

void interaccion(double *posicion, double *fuerzas, double *Vlj, double *Flj, double *r, double *r2, double dLt, double rc2, int N, double L)
{  int i,j,indice;
   double n2,Fx,Fy,Fz,Dx,Dy,Dz; 

   for (i=1; i<N;i++)
   {
    for (j=0; j<i;j++)    
    {
     Dx=deltax(posicion,i,j,0,L);
     Dy=deltax(posicion,i,j,1,L);
     Dz=deltax(posicion,i,j,2,L);
     n2=Dx*Dx*+Dy*Dy+Dz*Dz;

      if(n2<rc2)
      {       
       indice=(int)((n2-*(r+1))/dLt);
       Fx=Dx*(*(Flj+indice))/(*(r+indice));
       Fy=Dy*(*(Flj+indice))/(*(r+indice));
       Fz=Dz*(*(Flj+indice))/(*(r+indice));
//.... interpolación lineal
       Fx+=Dx*Dx*(*(Flj+indice+1)-*(Flj+indice+1))/dLt;                           // PREGUNTAR A GUILLE SI ESTO ES CORRECTO
       Fy+=Dy*Dy*(*(Flj+indice+1)-*(Flj+indice+1))/dLt;
       Fz+=Dz*Dz*(*(Flj+indice+1)-*(Flj+indice+1))/dLt;
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

double deltax(double *posicion, int i ,int j ,int p, double L)
{   
   double Dx; 
   Dx=*(posicion+3*i+p)-*(posicion+3*j+p);
   if(Dx>L/2.0) Dx=Dx-L;
   if(Dx<-L/2.0) Dx=Dx+L;

  return (double)Dx;
}

void verlet(double *posicion, double *v, double *fuerzas, double dt, double *Vlj,double *Flj, double *r, double *r2, double dLt, int N, double L, double rc2)
{ 
   int i;
   for (i=0; i<N;i++)//Primer paso de Verlet
   {   
    *(posicion+3*i)+=*(v+3*i)*dt+*(fuerzas+3*i)*dt*dt*0.5;
    *(posicion+3*i+1)+=*(v+3*i+1)*dt+*(fuerzas+3*i+1)*dt*dt*0.5;
    *(posicion+3*i+2)+=*(v+3*i+2)*dt+*(fuerzas+3*i+2)*dt*dt*0.5;
//  Aplico PBC
    if(*(posicion+3*i)>L) *(posicion+3*i)=*(posicion+3*i)-L;
    if(*(posicion+3*i+1)>L) *(posicion+3*i+1)=*(posicion+3*i+1)-L;
    if(*(posicion+3*i+2)>L) *(posicion+3*i+2)=*(posicion+3*i+2)-L;
    if(*(posicion+3*i)<0.0) *(posicion+3*i)=*(posicion+3*i)+L;
    if(*(posicion+3*i+1)<0.0) *(posicion+3*i+1)=*(posicion+3*i+1)+L;
    if(*(posicion+3*i+2)<0.0) *(posicion+3*i+2)=*(posicion+3*i+2)+L;

/*       ii=(i+dim-2)%(dim-2); ****PREGUNTAR A GUILLE SI PUEDO USAR %L*****
       jj=(j+dim-2)%(dim-2); */

    *(v+3*i)+=*(fuerzas+3*i)*dt*0.5;
    *(v+3*i+1)+=*(fuerzas+3*i+1)*dt*0.5;
    *(v+3*i+2)+=*(fuerzas+3*i+2)*dt*0.5;
   }
   interaccion(posicion,fuerzas,Vlj,Flj,r,r2,dLt,rc2,N,L); 
   for (i=0; i<N;i++)//Segundo paso de Verlet
   {   
    *(v+3*i)+=*(fuerzas+3*i)*dt*dt*0.5;
    *(v+3*i+1)+=*(fuerzas+3*i+1)*dt*dt*0.5;
    *(v+3*i+2)+=*(fuerzas+3*i+2)*dt*dt*0.5;
   }
}

double hamiltoniano(double *posicion, double *v, double *fuerzas, double dt, double *Vlj,double *Flj, double *r, double *r2, double dLt, double rc2, int N, double L)
{ 
   double p2,Vint,inter,Etot,Dx,Dy,Dz,n2;
   int i,j,indice;
  // término de energía cinética
   p2=0.0;
   inter=0.0;
   for (i=1; i<N;i++)
   {    
    p2+=*(v+3*i)*(*(v+3*i))+*(v+3*i+1)*(*(v+3*i+1))+*(v+3*i+2)*(*(v+3*i+2));
   }
   p2=0.5*p2;

  //término de interacción
   for (i=1; i<N;i++)
   {
    for (j=0; j<i;j++)    
    {
     Vint=0.0;
     Dx=deltax(posicion,i,j,0,L);
     Dy=deltax(posicion,i,j,1,L);
     Dz=deltax(posicion,i,j,2,L);
     n2=Dx*Dx*+Dy*Dy+Dz*Dz;
      if(n2<rc2)
     {
       indice=(int)((n2-*(r+1))/dLt);
       Vint=(*(Vlj+indice))/(*(r+indice));
//.... interpolación lineal
       Vint+=(*(Vlj+indice+1)-*(Vlj+indice+1))/dLt;  // PREGUNTAR A GUILLE SI ESTO ES CORRECTO
       inter+=Vint;
     }
   }
  }
   Etot=p2+inter;
   return Etot;
}




