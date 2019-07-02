#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*Metropolis por Pedro y Enzo. Julio 2019
Esta version incluye el cálculo de lambda para determinar
el tiempo de termalizacion*/

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
void interaccion(double *posicion, double *fuerzas, double *Vlj, double *Flj, \
double *r, double *r2, double dLt, double rc2, int N, double L);
double deltax(double *posicion1, double *posicion2, double L);
void verlet(double *posicion, double *v, double *fuerzas, double dt, \
double *Vlj,double *Flj, double *r, double *r2, double dLt, int N, double L, double rc2);
void hamiltoniano(double *posicion, double *v, double *fuerzas, double dt, \
double *Vlj,double *Flj, double *r, double *r2, double dLt, double rc2, \
int N, double L, double *Etot, double *Ecin, double *Epot);
void qlfta(double *fuerzas, int N);
int save_lammpstrj(char *filename, double* x, double* v, int N, double L, int frame);
float calc_T(double *v, int N);

void escaleo_v(double *v, int N, double *Tr, double Td);
double lambda(double *posicion, double a, int N, double Pi);

int main (int argc, char *argv[])
{ 
  int i,N,sfreq,tsteps;
  double *posicion,*v,*Vlj,*Flj,*r,*r2,*fuerzas,n,rc2, dt,L,tmax,T,Td,dLt;
  double *Etot,*Ecin,*Epot,*Tr,Pi,a,rho,Pid,V;
  int *semilla,gf;

   Pi=3.1415926535897932384626433832795;
   L=8.0;
   T= 2.0; //En 1.0 estoy cerca de la transicion de fase.
   tmax=10.0; //es razonable hasta 10.0
   N=27;
   gf=5000;
   tsteps=10000;

   if(argc==8)
    {sscanf(argv[1],"%lf",&L);
     sscanf(argv[2],"%lf",&T); //T inicial
     sscanf(argv[3],"%d",&N); 
     sscanf(argv[4],"%lf",&tmax);
     sscanf(argv[5],"%d",&tsteps);
     sscanf(argv[6],"%d",&gf);
     sscanf(argv[7],"%lf",&Td); //T deseada
     }

  rho=N/(L*L*L);
  dt=tmax/tsteps; //10^-3 como mìnimo
  sfreq=tsteps/1000;
  n=cbrt(N);
  a=L/n;
  V=1.0/(L*L*L);
  dLt=2.5/(gf); /*Grilla para interpolar potencial,
                              fuerzas y distancias, Ntablas=5000.0*/
  rc2=(double) 2.5*2.5;
  Vlj=(double*) malloc(gf*sizeof(double));
  Flj=(double*) malloc(gf*sizeof(double));
  r=(double*) malloc(gf*sizeof(double));
  r2=(double*) malloc(gf*sizeof(double));
  semilla=(int*) malloc(sizeof(int));
  Etot=(double*) malloc(sizeof(double));
  Ecin=(double*) malloc(sizeof(double));
  Epot=(double*) malloc(sizeof(double));
  Tr=(double*) malloc(sizeof(double));
  *semilla=S;
  posicion=(double*) malloc(3*N*sizeof(double));
  fuerzas=(double*) malloc(3*N*sizeof(double)); 
  v=(double*) malloc(3*N*sizeof(double));

 for(T=0.5;T<=2.0;T=T+0.15)
 {

  char filename[255]; 
  sprintf(filename,"Etot_L=%lf_N=%4.2d_T=%lf.dat",L,N,T);
  FILE *fp=fopen(filename,"w");

  sprintf(filename,"Ecin_L=%lf_N=%4.2d_T=%lf.dat",L,N,T);
  FILE *fp3=fopen(filename,"w");

  sprintf(filename,"Epot_L=%lf_N=%4.2d_T=%lf.dat",L,N,T);
  FILE *fp4=fopen(filename,"w");

  sprintf(filename,"Treal_L=%lf_N=%4.2d_T=%lf.dat",L,N,T);
  FILE *fp5=fopen(filename,"w");

  sprintf(filename,"Lambda_L=%lf_N=%4.2d_T=%lf.dat",L,N,T);
  FILE *fp6=fopen(filename,"w");
  
  char filename2[255]; 
  sprintf(filename2,"Visual_L=%lf_N=%4.2d.lammpstrj",L,N);
  FILE *fp2=fopen(filename2,"w");

  set_box(posicion,n,L);
  set_v(v, N, T, semilla);
  set_tablas(Vlj, Flj, r, r2, gf, dLt);

  qlfta(fuerzas, N);   
  interaccion(posicion,fuerzas,Vlj,Flj,r,r2,dLt,rc2,N,L); 

  for(i=0;i<tsteps;i++)
   {
   hamiltoniano(posicion,v,fuerzas,dt,Vlj,Flj,r,r2,dLt,rc2,N,L,\
       Etot,Ecin,Epot);
   fprintf(fp,"%lf %lf \n",i*dt,*Etot); 
   fprintf(fp3,"%lf %lf \n",i*dt,*Ecin); 
   fprintf(fp4,"%lf %lf \n",i*dt,*Epot); 
   *Tr=calc_T(v,N);
   fprintf(fp5,"%lf %lf \n",i*dt,*Tr); 
   fprintf(fp6,"%lf %lf \n",i*dt,lambda(posicion,a,N,Pi));
//   if(i%sfreq==0) imprimir(posicion,N,L,i*dt);
   if(i%sfreq==0) save_lammpstrj(filename2, posicion, v, N, L, i);

//   imprimir(posicion,N,L,i*dt);
   verlet(posicion, v, fuerzas, dt, Vlj, Flj,  r, r2, dLt, N, L, rc2);
 //  escaleo_v(v,N,Tr,Td);

   }
  }
  
  free(posicion);
  free(fuerzas);
  free(v);
  free(semilla);
  free(Vlj);
  free(Flj);
  free(r);
  free(r2);
//  fclose(fp);
 // fclose(fp2);
  return 0;
}

void imprimir (double *posicion, int N, double L,double t)
{ int i;
  char filename[255]; 
   sprintf(filename,"Posicion_L=%lf_N=%4.2d_t=%lf.dat",L,N,t);
   FILE *fp=fopen(filename,"w");
  for(i=0;i<N;i++)
  {
   fprintf(fp,"%f %f %f \n",*(posicion+3*i),*(posicion+3*i+1),*(posicion+3*i+2));       
  }
}

double lambda(double *posicion, double a, int N, double Pi)
{ int i;
  float lx,ly,lz;
  lx=0.0;
  ly=0.0;
  lz=0.0;
  for(i=0;i<N;i++)
  {
    lx+=cos((*(posicion+3*i)-a/2.0)*2.0*Pi/a);       
    ly+=cos((*(posicion+3*i+1)-a/2.0)*2.0*Pi/a);   
    lz+=cos((*(posicion+3*i+2)-a/2.0)*2.0*Pi/a);   
  }
  return (double)(lx+ly+lz)/(3.0*N);  
}

void qlfta(double *fuerzas, int N) //¡Que la fuerza te acompañe!
{ int i;
  for(i=0;i<3*N;i++)
  {
    *(fuerzas+i)=0.0;       
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
  double dL=L/n;
  i=0;
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

float calc_T(double *v, int N)
{ double vx2,vy2,vz2,Tr;
  int i;
   vx2=0.0;
   vy2=0.0;
   vz2=0.0;
   for (i=0; i<N;i++)
   {   
    vx2+=*(v+3*i)*(*(v+3*i));
    vy2+=*(v+3*i+1)*(*(v+3*i+1));
    vz2+=*(v+3*i+2)*(*(v+3*i+2));
   }
   vx2=vx2/N;
   vy2=vy2/N;
   vz2=vz2/N;
   Tr=(vx2+vy2+vz2)/3.0;
   return (double)Tr;
}

void escaleo_v(double *v, int N, double *Tr, double Td)
{ double esc;
  int i;
  esc=sqrt(Td/(*Tr));
   for (i=0; i<N;i++)
   {   
    *(v+3*i)=esc*(*(v+3*i));
    *(v+3*i+1)=esc*(*(v+3*i+1));
    *(v+3*i+2)=esc*(*(v+3*i+2));
   }
}
  

void set_tablas(double *Vlj,double *Flj, double *r, double *r2 , int gf, double dLt)
{  
   int i;
   for(i=1;i<gf+1;i++)
   {
    *(r+i-1)=i*dLt;
    *(r2+i-1)=(i*dLt)*(i*dLt);
    *(Vlj+i-1)=4.0*(pow(*(r2+i-1), -6.0)-pow(*(r2+i-1), -3.0))-\
    4.0*(pow(2.5, -12.0)-pow(2.5, -6.0));;
    *(Flj+i-1)=24.0*(2.0*pow(*(r2+i-1), -6.0)-pow(*(r2+i-1), -3.0))/(*(r2+i-1));
   }
}

void interaccion(double *posicion, double *fuerzas, double *Vlj, double *Flj, double *r, double *r2, double dLt, double rc2, int N, double L)
{  int i,j,indice;
   double n2,Fx,Fy,Fz,Dx,Dy,Dz,dr; 

   for (i=1; i<N;i++)
   {
    for (j=0; j<i;j++)    
    {
     Dx=deltax(posicion+3*i,posicion+3*j,L);
     Dy=deltax(posicion+3*i+1,posicion+3*j+1,L);
     Dz=deltax(posicion+3*i+2,posicion+3*j+2,L);
     n2=Dx*Dx+Dy*Dy+Dz*Dz;
      if(n2<rc2)
      {       
       n2=sqrt(n2);
       indice=(int)((n2-*(r))/dLt);
       
       dr=n2-*(r+indice);
       Fx=Dx*(*(Flj+indice));
       Fy=Dy*(*(Flj+indice));
       Fz=Dz*(*(Flj+indice));
//.... interpolación lineal
       Fx+=Dx*dr*(*(Flj+indice+1)-*(Flj+indice))/dLt; // PREGUNTAR A GUILLE SI ESTO ES CORRECTO
       Fy+=Dy*dr*(*(Flj+indice+1)-*(Flj+indice))/dLt; //Puedo no hacer la interpolacion para 5000
       Fz+=Dz*dr*(*(Flj+indice+1)-*(Flj+indice))/dLt; //valores en la tabla.-
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

double deltax(double *posicion1, double *posicion2, double L)
{   
   double Dx; 
   Dx=*(posicion1)-*(posicion2);
   if(Dx>L/2.0) Dx=Dx-L;
   if(Dx<-L/2.0) Dx=Dx+L;

  return (double)Dx;
}

void verlet(double *posicion, double *v, double *fuerzas, double dt, double *Vlj, \
double *Flj, double *r, double *r2, double dLt, int N, double L, double rc2)
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
/*    ****PREGUNTAR A GUILLE SI PUEDO USAR %L***** (NO SE PUEDE) */

    *(v+3*i)+=*(fuerzas+3*i)*dt*0.5;
    *(v+3*i+1)+=*(fuerzas+3*i+1)*dt*0.5;
    *(v+3*i+2)+=*(fuerzas+3*i+2)*dt*0.5;
   }

   qlfta(fuerzas, N);   
   interaccion(posicion,fuerzas,Vlj,Flj,r,r2,dLt,rc2,N,L); 
   
   for (i=0; i<N;i++)//Segundo paso de Verlet
   {   
    *(v+3*i)+=*(fuerzas+3*i)*dt*0.5;
    *(v+3*i+1)+=*(fuerzas+3*i+1)*dt*0.5;
    *(v+3*i+2)+=*(fuerzas+3*i+2)*dt*0.5;
   }
}

void hamiltoniano(double *posicion, double *v, double *fuerzas, double dt, \
 double *Vlj,double *Flj, double *r, double *r2, double dLt, double rc2, int N, \ 
 double L, double *Etot, double *Ecin, double *Epot)
{ 
   double p2,Vint,inter,Dx,Dy,Dz,n2,dr;
   int i,j,indice;
  // término de energía cinética
   p2=0.0;
   inter=0.0;
   for (i=0; i<N;i++)
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
     Dx=deltax(posicion+3*i,posicion+3*j, L);
     Dy=deltax(posicion+3*i+1,posicion+3*j+1, L);
     Dz=deltax(posicion+3*i+2,posicion+3*j+2, L);
     n2=Dx*Dx+Dy*Dy+Dz*Dz;
      if(n2<rc2)
     { 
       n2=sqrt(n2);
       indice=(int)((n2-*(r))/dLt);
       dr=n2-*(r+indice);
       Vint=(*(Vlj+indice));
//.... interpolación lineal
       Vint+=dr*(*(Vlj+indice+1)-*(Vlj+indice))/dLt;  // PREGUNTAR A GUILLE SI ESTO ES CORRECTO
       inter+=Vint;
     }
   }
  }
   *Ecin=p2/N;
   *Epot=inter/N;
   *Etot=(p2+inter)/N;
}


int save_lammpstrj(char *filename, double* x, double* v, int N, double L, int frame){
  FILE *fp;
  if (frame) fp = fopen(filename, "a"); // Si frame==0, es el primero y por lo tanto
  else fp = fopen(filename, "w");       // tiene que crear un nuevo archivo
  // Header que usa lammps
	fprintf(fp, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n", frame, N);
	for(int l = 0; l < 3; l++){
		fprintf(fp, "0 %f\n", L); // Limites de la caja en x-y-z
	}
	fprintf(fp, "ITEM: ATOMS id x y z vx vy vz \n"); // "Nombre de las columnas"
	for(int i = 0; i < N; i++){
		fprintf(fp, "%d %f %f %f %f %f %f\n", i, x[3*i], x[3*i+1], x[3*i+2], v[3*i], v[3*i+1], v[3*i+2]);
	}
  fclose(fp);
  return 0;
}
