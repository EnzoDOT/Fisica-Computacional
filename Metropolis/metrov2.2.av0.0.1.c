#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*Percolar por Pedro y Enzo. Marzo 2019*/

#define M 2147483647
#define A 16807
#define Q 127773
#define R 2836  
#define S 2


float irandom(int *semilla);
void poblar (int *red,float p, int dim, int *semilla);
void imprimir (int *red, int dim);
void cborde (int *red, int dim);
void hamiltoniano (int *red, int dim, float Jint, float B);
void flipeo (int *red, int dim, float *tabla);


int main (int argc, char *argv[])
{ int *red;
  FILE *fp;
  int dim;
  int i,N;
  float p,B,Jint,*tabla;
  int *semilla;
  semilla=(int*) malloc(sizeof(int));
  fp=fopen("MMC.dat","a");
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
  dim=dim+2;
  red=(int*) malloc(dim*dim*sizeof(int));  
  tabla=(float*) malloc(5*sizeof(float));
  *(tabla)= -8*Jint;
  *(tabla+1)= -4*Jint;
  *(tabla+2)= 0.0;
  *(tabla+3)= 4*Jint;
  *(tabla+4)= 8*Jint;
    
//  imprimir(red,dim);
  
 // for(i=0; i<N;i++)
 //  {
   *semilla=S;
   poblar(red,p,dim,semilla);
   cborde(red,dim);
   imprimir(red,dim);
   hamiltoniano(red,dim,Jint,B);
   flipeo(red,dim,tabla);
//   }
//  printf("%f",prom);  
  free(red);
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
          if(*(red+dim*i+j)==1)
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
    *(red+(i+1)*dim-1)=*(red+i*dim+1); //cborde para ultima columna
    
    *(red+i)=*(red+dim*dim-2*dim+i); //cborde para primer fila
    *(red+dim*dim-dim+i)=*(red+dim+i); //cborde para ultima fila
    }
}

void hamiltoniano (int *red, int dim, float Jint, float B)
{ int i,j;          
  float energiaB;
  energiaB=0.0;
  for(j=1;j<dim-1;j++)
   {for(i=1;i<dim-1;i++)
    {      
    energiaB=energiaB+*(red+i+j*dim); //cborde para primer fila
    }
   }
    energiaB=B*energiaB;
//    printf("%f ",energiaB);
}

void flipeo (int *red, int dim, float *tabla)
{ int i,j,suma;
  int Delta;  
        
  for(j=1;j<dim-1;j++)
  {for(i=1;i<dim-1;i++)
   {      
    flip1=*(red+i+j*dim)*(-1);
    suma=*(red+i-1+j*dim)+*(red+i+1+j*dim)+*(red+i+(j-1)*dim)+*(red+i+(j+1)*dim);
    Delta=*(red+i+j*dim)*(-1)*2*suma
   }
  }
  
    
    //printf("%d %d",suma,*(red+i+j*dim));

}


