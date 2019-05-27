#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*Percolar por Pedro y Enzo. Mayo 2019*/

#define M 2147483647
#define A 16807
#define Q 127773
#define R 2836  
#define S 2


float irandom(int *semilla);
void poblar (int *red,float p, int dim, int *semilla);
void imprimir (int *red, int dim);
void cborde (int *red, int dim);
float hamiltoniano (int *red, int dim, float Jint, float B);
void flipeo (int *red, int dim, float Jint, float B, float *energia, float *mag, int *semilla, float *pacep,FILE *fp, FILE *fp2);
void aceptacion (int *red, int dim, float Delta, float *energia, float *mag, int i, int j, int *semilla, float *pacep, FILE *fp, FILE *fp2);
float magnetizacion (int *red, int dim);

int main (int argc, char *argv[])
{ int *red;
//  FILE *fp,*fp2;
  char filename[255]; 
  int dim;
  int i,N;
  float p,B,Jint,*pacep,*energia,*mag;
  int *semilla;
  semilla=(int*) malloc(sizeof(int));
  pacep=(float*) malloc(sizeof(float));
  energia=(float*) malloc(sizeof(float));
  mag=(float*) malloc(sizeof(float));
//  fp=fopen("MMCEn.dat","w");
//  fp2=fopen("MMCMag.dat","w");
   *semilla=S;
   p=0.5;
   dim=8;
   N=100;
   B=1.0;
   Jint=0.0;
   if(argc==6)
    {sscanf(argv[1],"%d",&dim);
     sscanf(argv[2],"%f",&p);
     sscanf(argv[3],"%d",&N);
     sscanf(argv[4],"%f",&Jint);
     sscanf(argv[5],"%f",&B);
     }
   sprintf(filename,"MMCEn_L=%d_J=%f_B=%f.dat",dim,Jint,B);
   FILE *fp=fopen(filename,"w");
   sprintf(filename,"MMCMag_L=%d_J=%f_B=%f.dat",dim,Jint,B);
   FILE *fp2=fopen(filename,"w");
   *mag=0.0;
   *energia=0.0;

   red=(int*) malloc(dim*dim*sizeof(int));  
   poblar(red,p,dim,semilla);
   cborde(red,dim);
   *energia=hamiltoniano(red,dim,Jint,B);
   *mag=magnetizacion(red,dim);
   imprimir(red,dim);
   *pacep=0.0;
   
  for(i=0; i<N;i++)
   {
//   cborde(red,dim);
   flipeo (red,dim,Jint,B,energia,mag,semilla,pacep,fp,fp2); 
//   cborde(red,dim); 
   }
//   cborde(red,dim);
  imprimir(red,dim);
  free(red);
  free(semilla);
  fclose(fp);
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

void poblar (int *red,float p, int dim, int *semilla)
{ int i;
  for(i=0;i<dim*dim;i++)
      {*(red+i)=-1;
        {if(irandom(semilla)<p)
         *(red+i)=1;
        }
       }      
}

void imprimir (int *red, int dim)
{ int i,j;
  for(i=0;i<dim;i++)
        {for(j=0;j<dim;j++)
         {
          if(*(red+dim*i+j)==1 || *(red+dim*i+j)==0)
          {
           printf(" %d ",*(red+dim*i+j));
          }
          else
          {
           printf("%d ",*(red+dim*i+j));                   
          }
         }
         printf("\n");
        }
 printf("\n");
}

void cborde (int *red, int dim)
{ int i;          
  for(i=0;i<dim;i++)
    {  
    *(red+i*dim)=*(red+(i+1)*dim-2); //cborde para primer columna
    *(red+i)=*(red+dim*dim-2*dim+i); //cborde para primer fila
    *(red+(i+1)*dim-1)=*(red+i*dim+1); //cborde para ultima columna
    *(red+dim*dim-dim+i)=*(red+dim+i); //cborde para ultima fila
    }
}


void flipeo (int *red, int dim, float Jint, float B, float *energia, float *mag, int *semilla, float *pacep,FILE *fp, FILE *fp2)
{ int i,j,suma;
  float Delta;  
  for(j=1;j<dim-1;j++)
   {
   for(i=1;i<dim-1;i++)
    {
     suma=*(red+i-1+j*dim)+*(red+i+1+j*dim)+*(red+i+(j-1)*dim)+*(red+i+(j+1)*dim);
     Delta=*(red+i+j*dim)*2*suma*Jint+*(red+i+j*dim)*2*B;
     aceptacion(red,dim,Delta,energia,mag,i,j,semilla,pacep,fp,fp2);
    }
   }
}   

void aceptacion (int *red, int dim, float Delta, float *energia, float *mag, int i, int j, int *semilla, float *pacep, FILE *fp, FILE *fp2)
{     float p;
//      int ii,jj;
      p=exp(-Delta);
      *pacep=*pacep+1;
      if(irandom(semilla)<=p)  
      {  
//      ii=(i+dim-1)%(dim-1);
//      jj=(j+dim-1)%(dim-1);
       *(red+i+j*dim)=*(red+i+j*dim)*(-1);
 //      *(red+ii+jj*dim)=*(red+i+j*dim);
      if(j==1 || j==dim-2 || i==1 || i==dim-2) cborde(red,dim);     
      *energia=*energia-Delta;
      *mag=*mag+*(red+i+j*dim);   
      fprintf(fp2," %f %f \n", *pacep, *mag);
      fprintf(fp,"%f %f \n", *pacep, *energia); 
      }
}

float hamiltoniano (int *red, int dim, float Jint, float B)
{ int i,j;          
  float energiaB,energiaJ;
  energiaB=0.0;
  energiaJ=0.0;
  for(j=1;j<dim-1;j++)
   {for(i=1;i<dim-1;i++)
    {      
    energiaB=energiaB+*(red+i+j*dim); 
    energiaJ=energiaJ+*(red+i+1+j*dim)+*(red+i+(j+1)*dim);
    }
   }
    energiaB=B*energiaB+Jint*energiaJ;
    return (float)energiaB;    
}

float magnetizacion (int *red, int dim)
{ int i,j;          
  float stot;
  stot=0;
  for(j=1;j<dim-1;j++)
   {for(i=1;i<dim-1;i++)
    {      
    stot=stot+*(red+i+j*dim);
    }
   }
  return (float)stot/((float)(dim-2)*(dim-2));
}
