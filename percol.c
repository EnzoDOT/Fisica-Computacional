#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*Percolar por Pedro y Enzo. Marzo 2019*/

#define M 2147483647
#define A 16807
#define Q 127773
#define R 2836  
#define S 591987


float irandom(int *semilla);
void poblar (int *red,float p, int dim, int *semilla);
void imprimir (int *red, int dim);
void clasificar(int *red,int i);
void etiqueta_falsa(int *red,int *historial,int s1, int s2, int i);

int main ()
{ int *red;
  int dim;
  float p;
  float *x;
  int *semilla;
  int i;
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;
  dim=8;
//  *x=irandom(semilla);
//  sscanf(argv[1],"%d",&dim);
  *semilla=S;
  red=(int*) malloc(dim*dim*sizeof(int));
  p=0.5;
  poblar(red,p,dim,semilla);
//  clasificar(red,i);
  //percola();
  imprimir(red,dim);
  free(red);
  free(semilla);
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
//  printf("%f",x);
  return (float)x;
}

void poblar (int *red,float p, int dim, int *semilla)
{ int i;
  float *x;
  for(i=0;i<dim*dim;i++)
      {*(red+i)=0;
        *x=irandom((int*) semilla);
        {if(*x<p)
         *(red+i)=1;
        }
       }
}

void imprimir (int *red, int dim)
{ int i,j;
  for(i=0;i<dim;i=i+1)
        {for(j=0;j<dim;j=j+1)
         {printf("%d",*(red+dim*i+j));
         }
         printf("\n");
        }
 printf("\n");
}

void clasificar(int *red,int i)
{ int frag;
  int s1;
  int dim;
  int *historial;
  historial=(int*) malloc(dim*dim*sizeof(int));

  *(red+i)=0;
  frag=2;
  if(*(red+i)==1)
   {*(red+i)=frag;
    frag=frag+1;
   }
  s1=*(red+i-1); 
  if(*(red+i))
   {*(red+i)=s1;
   }
  while(*(historial+s1)<0)
   {s1=-*(historial+s1);
   }
}

void etiqueta_falsa(int *red,int *historial,int s1, int s2, int i)
{int minimo, maximo;
 while(*(historial+s1)<0)
  {s1=-(*historial+s1);
  }
  while(*(historial+s2)<0)
  {s2=-(*(historial+s2));
  }
  if(s1<s2){minimo=s1;
            maximo=s2;
            }
  else{minimo=s2;
       maximo=s1;
       }
  *(red+i)=minimo;
  *(historial+maximo)=-minimo;
  *(historial+minimo)=minimo;
  if(minimo==maximo){*(historial+maximo)=minimo;
                     }                
}


