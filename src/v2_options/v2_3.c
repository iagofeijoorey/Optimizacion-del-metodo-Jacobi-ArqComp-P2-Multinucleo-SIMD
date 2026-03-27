#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


#include "counter.h"


#define tol 1e-5
#define max_iter 15000
#define TAM_LINEA_CACHE 64
#define tamanho_desenrollo 5 // Desenrollamos los lazos de 5 en 5


int n;
int n_real;


double resolver_jacobi(double **a, double *b, double *x);
void generar_matriz(double **matriz);
void generar_vector(double *vector);




int main(int argc, char *argv[])
{
   int c = 1;
   n_real = n - (n % tamanho_desenrollo); // Para gestionar el caso en el que n no es múltiplo de tamanho_desenrollo


   if (argc != 2 && argc != 3)
   {
       printf("Uso incorrecto: %s n [c]\n", argv[0]);
       exit(1);
   }


   n = atoi(argv[1]);
   if (argc == 3)
       c = atoi(argv[2]);


   printf("%s: n=%d c=%d\n", argv[0], n, c);


   /* Reserva de matriz alineada */
   double **a = aligned_alloc(TAM_LINEA_CACHE, sizeof(double *) * n);
   double *data = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n * n);


   // Desenrollamos el lazo (5 en 5)
   for (int i = 0; i < n_real; i += tamanho_desenrollo)
   {
       a[i] = &data[i * n];
       a[i+1] = &data[(i+1) * n];
       a[i+2] = &data[(i+2) * n];
       a[i+3] = &data[(i+3) * n];
       a[i+4] = &data[(i+4) * n];
   }


   // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
   for (int i = n_real; i < n; i++) a[i] = &data[i * n];


   /* Reserva de vectores */
   double *b = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n);
   double *x = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n);


   /* Inicializar x */
   memset(x, 0, n * sizeof(double));


   srand(time(NULL));


   generar_matriz(a);
   generar_vector(b);


   double ck;


   start_counter();
   double solucion = resolver_jacobi(a, b, x);
   ck = get_counter();


   printf("%s: solucion = %e\n", argv[0], solucion);
   printf("%s: counter = %lf\n", argv[0], ck);


   free(a);
   free(data);
   free(b);
   free(x);


   return 0;
}




double resolver_jacobi(double **a, double *b, double *x)
{
   double norm2, sigma;
   double x_new[n];


   for (int iter = 0; iter < max_iter; iter++)
   {
       norm2 = 0;


       for (int i = 0; i < n; i++)
       {
           sigma = 0.0;


           // Desenrollamos el lazo (5 en 5)
           for (int j = 0; j < n_real; j += tamanho_desenrollo)
           {
               if (i != j) sigma += a[i][j] * x[j];
               if (i != (j+1)) sigma += a[i][j+1] * x[j+1];
               if (i != (j+2)) sigma += a[i][j+2] * x[j+2];
               if (i != (j+3)) sigma += a[i][j+3] * x[j+3];
               if (i != (j+4)) sigma += a[i][j+4] * x[j+4];
           }


           // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
           for (int j = n_real; j < n; j++) if (i != j) sigma += a[i][j] * x[j];


           x_new[i] = (b[i] - sigma) / a[i][i];
           norm2 += pow(x_new[i] - x[i], 2);
       }
      
       /* copiar nuevo vector */
       // Desenrollamos el lazo (5 en 5)
       for (int i = 0; i < n_real; i += tamanho_desenrollo)
       {
           x[i] = x_new[i];
           x[i+1] = x_new[i+1];
           x[i+2] = x_new[i+2];
           x[i+3] = x_new[i+3];
           x[i+4] = x_new[i+4];
       }
      
       // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
       for (int i = n_real; i < n; i++) x[i] = x_new[i];
          
       if (sqrt(norm2) < tol)
           return norm2;
   }


   return -1;
}




/* Genera matriz diagonalmente dominante */
void generar_matriz(double **matriz)
{
   for (int i = 0; i < n; i++)
   {
       double suma = 0.0;
       // Desenrollamos el lazo para realizar las operaciones de 5 en 5
       for (int j = 0; j < n_real; j += tamanho_desenrollo)
       {
           matriz[i][j] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j];


           matriz[i][j+1] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j+1];


           matriz[i][j+2] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j+2];


           matriz[i][j+3] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j+3];


           matriz[i][j+4] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j+4];
       }


       // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
       for (int j = n_real; j < n; j++)
       {
           matriz[i][j] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j];
       }
       // La diagonal es la suma de toda la fila (garantiza dominancia estricta)
       matriz[i][i] += suma;
   }
}




/* Genera vector aleatorio */
void generar_vector(double *vector)
{
   // Modificamos para desenrollar el lazo de 5 en 5
   for (int i = 0; i < n_real; i += tamanho_desenrollo)
   {
       vector[i] = ((double)rand() / RAND_MAX) * 100.0;
       vector[i+1] = ((double)rand() / RAND_MAX) * 100.0;
       vector[i+2] = ((double)rand() / RAND_MAX) * 100.0;
       vector[i+3] = ((double)rand() / RAND_MAX) * 100.0;
       vector[i+4] = ((double)rand() / RAND_MAX) * 100.0;
   }
  
   // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
   for (int i = n_real; i < n; i++) vector[i] = ((double)rand() / RAND_MAX) * 100.0;
  
}







