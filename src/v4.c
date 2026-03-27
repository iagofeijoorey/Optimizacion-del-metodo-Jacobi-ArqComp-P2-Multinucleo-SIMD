#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <immintrin.h>

#include "counter.h"


#define tol 1e-5
#define max_iter 15000
#define TAM_LINEA_CACHE 64
#define tamanho_desenrollo 4 // Desenrollamos los lazos de 5 en 5


int n;
int tamanho_bloque;


double resolver_jacobi(double **a, double *b, double *x);
void generar_matriz(double **matriz);
void generar_vector(double *vector);




int main(int argc, char *argv[])
{
   int c = 1;
   if (argc != 2 && argc != 3)
   {
       printf("Uso incorrecto: %s n [c]\n", argv[0]);
       exit(1);
   }


   n = atoi(argv[1]);
   tamanho_bloque = n - (n % tamanho_desenrollo); // Para gestionar el caso en el que n no es múltiplo de tamanho_desenrollo
   if (argc == 3)
       c = atoi(argv[2]);


   printf("%s: n=%d c=%d\n", argv[0], n, c);


   /* Reserva de matriz alineada */
   double **a = aligned_alloc(TAM_LINEA_CACHE, sizeof(double *) * n);
   double *data = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n * n);


   // Desenrollamos el lazo (5 en 5)
   for (int i = 0; i < tamanho_bloque; i += tamanho_desenrollo)
   {
       a[i] = &data[i * n];
       a[i+1] = &data[(i+1) * n];
       a[i+2] = &data[(i+2) * n];
       a[i+3] = &data[(i+3) * n];
   }


   // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
   for (int i = tamanho_bloque; i < n; i++) a[i] = &data[i * n];


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


           // Acumulamos en SIMD y reducimos una sola vez al final del bucle.
           __m256d vsum = _mm256_setzero_pd();
           for (int j = 0; j < tamanho_bloque; j += tamanho_desenrollo)
           {
               if (i >= j && i < j + tamanho_desenrollo)
               {
                   for (int k = j; k < j + tamanho_desenrollo; k++)
                       if (k != i) sigma += a[i][k] * x[k];
               }
               else
               {
                   __m256d v1 = _mm256_loadu_pd(&a[i][j]);
                   __m256d v2 = _mm256_loadu_pd(&x[j]);
                   __m256d v3 = _mm256_mul_pd(v1, v2);
                   vsum = _mm256_add_pd(vsum, v3);
               }

           }

           __m128d low = _mm256_castpd256_pd128(vsum);          //Pilla los dos primeos doubles
           __m128d high = _mm256_extractf128_pd(vsum, 1);       //PIlla los dos ultimos doubles
           __m128d sum2 = _mm_add_pd(low, high);                //Suma la parte baja con la alta y lo fuarda en una nueva variable
           sum2 = _mm_hadd_pd(sum2, sum2);                      //Ahora suma ambas partes y las guarda en los dos slots de sum2
           sigma += _mm_cvtsd_f64(sum2);                        //Como sum 2 tiene el resultado duplicado, basta con guardar la priemra parte en sigma


           // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
           for (int j = tamanho_bloque; j < n; j++) if (i != j) sigma += a[i][j] * x[j];


           x_new[i] = (b[i] - sigma) / a[i][i];
           // Cambiamos el pow de la version 1 por esta multiplicacion, evitando así llamar a pow, lo que puede
           //  reducir el número de instrucciones y mejorar el rendimiento
           double diff = x_new[i] - x[i];
           norm2 += diff * diff;
       }
      
       /* copiar nuevo vector */
       // Desenrollamos el lazo (5 en 5)
       for (int i = 0; i < tamanho_bloque; i += tamanho_desenrollo)
       {
           __m256d v1 =  _mm256_loadu_pd(&x_new[i]);
           _mm256_storeu_pd(&x[i], v1);

           //x[i] = x_new[i];

       }
      
       // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
       for (int i = tamanho_bloque; i < n; i++) x[i] = x_new[i];
          
       if (norm2 < tol * tol) // Comparamos con tol^2 para evitar la llamada a sqrt, lo que puede reducir el número de instrucciones
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
       for (int j = 0; j < tamanho_bloque; j += tamanho_desenrollo)
       {
           matriz[i][j] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j];


           matriz[i][j+1] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j+1];


           matriz[i][j+2] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j+2];


           matriz[i][j+3] = ((double)rand() / RAND_MAX) * 100.0;
           suma += matriz[i][j+3];
       }


       // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
       for (int j = tamanho_bloque; j < n; j++)
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
   for (int i = 0; i < tamanho_bloque; i += tamanho_desenrollo)
   {
       vector[i] = ((double)rand() / RAND_MAX) * 100.0;
       vector[i+1] = ((double)rand() / RAND_MAX) * 100.0;
       vector[i+2] = ((double)rand() / RAND_MAX) * 100.0;
       vector[i+3] = ((double)rand() / RAND_MAX) * 100.0;
   }
  
   // Gestionamos el caso en el que n no es múltiplo de tamanho_desenrollo
   for (int i = tamanho_bloque; i < n; i++) vector[i] = ((double)rand() / RAND_MAX) * 100.0;
  
}





