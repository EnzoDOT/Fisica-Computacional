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
void poblar (int *red,float p, int dim, int *semilla);
void imprimir (int *red, int dim);

int main (int argc, char *argv[])
{ int *red,*ns,*ns2,*historial;
  FILE *fp;
  int dim;
  int per,perc;
  int i,j,N;
  float p,x,delta,xp;
  int *semilla;
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;

   fp=fopen("miarchivo.dat","w");
   p=0.5;
   dim=8;
   N=100;
   red=(int*) malloc(dim*dim*sizeof(int));  
   ns=(int*) malloc(dim*dim*sizeof(int));         
   ns2=(int*) malloc(dim*dim*sizeof(int));    
   historial=(int*) malloc(dim*dim*sizeof(int));     
    j=0;
    xp=0.0;
    delta=1.0;
    for(i=0;i<27000;i++)
    {
    x=irandom(semilla);

    
     printf("%f \n", x);   
    
    }

  free(red);
  free(ns);
  free(ns2);
  free(historial);
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
  x=-1.0+(*semilla)*(2.0/M);
  return (float)x;
}

void poblar (int *red,float p, int dim, int *semilla)
{ int i;
  for(i=0;i<dim*dim;i++)
      {*(red+i)=0;
        {if(irandom(semilla)<p)
         *(red+i)=1;
        }
       }      
}

void imprimir (int *red, int dim)
{ int i,j;
  int width =1;
  char letter = ' ';
  for(i=0;i<dim;i++)
        {for(j=0;j<dim;j++)
         {printf("%d",*(red+dim*i+j));
          printf( "%*c", width, letter );
         }
         printf("\n");
        }
 printf("\n");
}

