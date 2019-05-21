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

int main (int argc, char *argv[])
{ int *red;
  FILE *fp;
  int dim;
  int i,N;
  float p;
  int *semilla;
  semilla=(int*) malloc(sizeof(int));
  fp=fopen("MMC.dat","a");
   *semilla=S;
   p=0.5;
   dim=8;
   N=100;
   if(argc==4)
    {sscanf(argv[1],"%d",&dim);
     sscanf(argv[2],"%f",&p);
     sscanf(argv[3],"%d",&N);
     }
  red=(int*) malloc(dim*dim*sizeof(int));  
  poblar(red,p,dim,semilla);
  cborde(red,dim);
  imprimir(red,dim);
  for(i=0; i<N;i++)
   {
   *semilla=S+i;
   poblar(red,p,dim,semilla);
   }
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
