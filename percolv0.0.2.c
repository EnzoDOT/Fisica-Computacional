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
void clasificar(int *red,int dim);
void etiqueta_falsa(int *red,int *historial,int s1, int s2, int i);
void actualizar(int *red, int *historial,int s,int frag);
void corregir_etiqueta(int *red,int *historial,int dim);

int main (int argc,char *argv[])
//int main ()
{ int *red;
  int dim;
  float p;
  int *semilla;
  semilla=(int*) malloc(sizeof(int));
  *semilla=S;
  p=0.5;
  dim=8;
// *x=irandom(semilla);
  if(argc==3)
   {sscanf(argv[1],"%d",&dim);
    sscanf(argv[2],"%f",&p);
    }
  *semilla=S;
  red=(int*) malloc(dim*dim*sizeof(int));
  poblar(red,p,dim,semilla);
  imprimir(red,dim);
  clasificar(red,dim);
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

void clasificar(int *red, int dim)
{ int frag;
  int i,j;
  int s1,s2;
  int *historial;
  historial=(int*) malloc(dim*dim*sizeof(int));

  for(j=0;j<dim*dim;j++) //Inicializo historial
    {*(historial+j)=j;
    }
    
   s1=0;
   s2=0;

  frag=2;
  if(*(red)==1) //El primer elemento de la red.
   {*(red)=frag;
    frag++;
   }
   
   
  for(i=1;i<dim;i++) //barro la primera fila
   {
    s1=*(red+i-1);
    if((*(red+i)!=0) & (s1!=0)) //red de i está ocupado y s1 distinto de 0.
     {
      actualizar(red+i,historial,s1,frag);
     }
    if((*(red+i)!=0) & (s1==0))   //red de i está ocupado y s1=0.
     { 
      *(red+i)=frag; //asigno una nueva etiqueta
       frag++; //actualizo etiqueta
     }
    }

  for(j=1;j<dim;j++) //ahora recorro la primer fila
    {
     s2=*(red+(j-1)*dim);
     if((*(red+j*dim)!=0) & (s2!=0)) //red está ocupado y s2 distinto de 0
     {//valor igual al vecino de la fila anterior<
      actualizar(red+j*dim,historial,s2,frag);     
     }
     if((*(red+dim*j)!=0) & (s2==0))
     {*(red+dim*j)=frag;
      frag++;
     }
     
     for(i=1;i<dim;i++)
     {s1=*(red+dim*j+i-1);
      s2=*(red+(j-1)*dim+i);
     if(*(red+dim*j+i))
     {if(s1*s2)
       {
      etiqueta_falsa(red,historial,s1,s2,dim*j+i);
       }
      else
       {if(s1!=0)
        {
         actualizar(red+dim*j+i,historial,s1,frag);
        }
       else 
        {if(s2!=0)
          {
            actualizar(red+dim*j+i,historial,s2,frag);
          }
         else
          {
            *(red+dim*j+i)=frag;
            frag++;
          }
        }
       }   
     }
     }
    }
   
   corregir_etiqueta(red,historial,dim);

   imprimir(historial,dim);
}

void actualizar(int *local,int *historial,int s,int frag) // modificado por Guillermo Frank
{ 
 while(*(historial+s)<0)
 {s=-(*(historial+s));
 }
 *local=s;
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
  if(minimo==maximo)
    {*(historial+maximo)=minimo;
    }                
}

void corregir_etiqueta(int *red,int *historial,int dim) // Guillermo Frank
{
  int i,s;

 for(i=0;i<dim*dim;i++)
   {
     s=(*(red+i));
     if(s)
      {
        while(*(historial+s)<0) s=-(*(historial+s));
      }
     *(red+i)=s;
   }
}

