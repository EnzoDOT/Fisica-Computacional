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
void clasificar(int *red,int dim, int *historial);
void etiqueta_falsa(int *red,int *historial,int s1, int s2, int i);
void actualizar(int *red, int *historial,int s,int frag);
void corregir_etiqueta(int *red,int *historial,int dim);
void percola(int *red,int dim,int *per, int *cper);
void distfrag (int *red, int dim,int *ns, int *ns2, int *cper);

int main (int argc, char *argv[])
{ int *red,*ns,*ns2,*historial;
  int *nsp;
  FILE *fp2,*fp;
  int dim,perc;
  int per;
  int i,k,N;
  float p,m2,pc,gm;
  int *semilla,*cper;
  semilla=(int*) malloc(sizeof(int));
  cper=(int*) malloc(sizeof(int));
   gm=-2.64;
   p=0.5;
   dim=8;
   N=100;
   if(argc==4)
    {sscanf(argv[1],"%d",&dim);
     sscanf(argv[2],"%f",&p);
     sscanf(argv[3],"%d",&N);
     }
   fp2=fopen("Problema6L_128.dat","w"); 
//   fp=fopen("Problema6gamma.dat","w"); 
   red=(int*) malloc(dim*dim*sizeof(int));  
   ns=(int*) malloc(dim*dim*sizeof(int));         
   ns2=(int*) malloc(dim*dim*sizeof(int));    
   nsp=(int*) malloc(dim*dim*sizeof(int));   
   historial=(int*) malloc(dim*dim*sizeof(int));     

   perc=0;

// pc(L=6)=0.574717  dispersionA
// pc(L=128)=0.592488


   pc=0.574717 ;

   for(p=0.0;p<=1.0;p=p+0.01)
   {
    for(k=0;k<dim*dim;k++)
     {
     *(nsp+k)=0;
     }
    perc=0;
    for(i=0;i<N;i++)
    {
    *semilla=S+i;
    *cper=0;
    per=0;      
    poblar(red,p,dim,semilla);
    clasificar(red,dim,historial);
    percola(red,dim,&per,cper);
    perc=perc+per;
    distfrag(red, dim, ns, ns2, cper); 
    for(k=0;k<dim*dim;k++)
      {
       *(nsp+k)=*(nsp+k)+*(ns2+k);
      }
    }
    m2=0.0;
//    fprintf(fp2,"%f ", p);
    for(k=0;k<dim*dim;k++)
      {
      m2=m2+*(nsp+k)*(1.0)*k*k/N;
      }
      fprintf(fp2,"%f %f \n", p, m2);
   //  fprintf(fp,"%f %f \n", p-pc,6.0*powf(fabs((p-pc)/pc),gm));
   }
  free(red);
  free(ns);
  free(ns2);
  free(historial);
  free(semilla);
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

void clasificar(int *red, int dim, int *historial)
{ int frag;
  int i,j;
  int s1,s2;

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


void percola (int *red, int dim, int *per, int *cper)
{ int i,j;    
  for(i=0;i<dim;i++)
      {
       for(j=0;j<dim;j++)
       {  
         if(*(red+i)!=0 && *(red+i)==*(red+dim*dim-dim+j))
         {
          *per=1;
           *cper=*(red+i);
         }
       }
      }             
}

void distfrag (int *red, int dim, int *ns, int *ns2, int *cper)
{ int i,j;    
  for(j=0;j<dim*dim;j++)
    {  
     *(ns+j)=0;
     *(ns2+j)=0;
    }

  for(i=0;i<dim*dim;i++) //Determina el tamaño de cada cluster
       {       
       if(*(red+i)!=0 && *(red+i)!=*cper)
         {
          j=*(red+i);          
          *(ns+j)=*(ns+j)+1;        
         }        
        }     

  for(i=0;i<dim*dim;i++) //Determinamos la cantidad de clusters de tamaño s                    
       {                 //y construimos el vector ns2 de distribuciòn de tamaños
       if(*(ns+i)!=0)
         {
          j=*(ns+i);          
          *(ns2+j)=*(ns2+j)+1;
         }        
        }  


}
