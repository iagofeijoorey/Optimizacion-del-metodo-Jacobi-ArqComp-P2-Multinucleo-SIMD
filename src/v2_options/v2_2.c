#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


#include "counter.h"


#define tol 1e-5
#define max_iter 15000
#define TAM_LINEA_CACHE 64


int n;


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
   if (argc == 3)
       c = atoi(argv[2]);


   printf("%s: n=%d c=%d\n", argv[0], n, c);


   /* Reserva de matriz alineada */
   double **a = aligned_alloc(TAM_LINEA_CACHE, sizeof(double *) * n);
   double *data = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n * n);


   for (int i = 0; i < n; i++)
       a[i] = &data[i * n];


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
       // Dividimos el bucle original en dos partes:


       // 1: Calculamos el nuevo vector solucion
       for (int i = 0; i < n; i++)
       {
           sigma = 0.0;


           for (int j = 0; j < n; j++)
           {
               if (i != j)
                   sigma += a[i][j] * x[j];
           }


           x_new[i] = (b[i] - sigma) / a[i][i];
       }


       // 2: Calculamos la norma del vector al cuadrado
       norm2 = 0;
       for (int i = 0; i < n; i++)
       {
           norm2 += pow(x_new[i] - x[i], 2);
       }


       /* copiar nuevo vector */
       for (int i = 0; i < n; i++)
           x[i] = x_new[i];


       if (sqrt(norm2) < tol)
           return norm2;
   }


   return -1;
}




/* Genera matriz diagonalmente dominante */
void generar_matriz(double **matriz)
{
   // Dividimos el lazo en dos partes:


   // 1: generamos los valores aleatorios para la matriz
   for (int i = 0; i < n; i++)
   {
       for (int j = 0; j < n; j++)
       {
           matriz[i][j] = ((double)rand() / RAND_MAX) * 100.0;
       }
   }


   // 2: Calculamos la suma para garantizar que la diagonal es dominante
   for (int i = 0; i < n; i++)
   {
       double suma = 0.0;
       for (int j = 0; j < n; j++)
       {
           suma += matriz[i][j];
       }


       // La diagonal es la suma de toda la fila (garantiza dominancia estricta)
       matriz[i][i] += suma;
   }
}




/* Genera vector aleatorio */
void generar_vector(double *vector)
{
   for (int i = 0; i < n; i++)
       vector[i] = ((double)rand() / RAND_MAX) * 100.0;
}





