#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "counter.h"

#define tol 1e-5
#define max_iter 15000
#define TAM_LINEA_CACHE 64

int n;

double resolver_jacobi(double** a, double *b, double* x);
void generar_matriz(double** matriz);
void generar_vector(double* vector);

int main(int argc, char* argv[])
{
    int c = 1;

    if (argc != 2 && argc != 3)
    {
        printf("Uso incorrecto: ./v1 n [c]\n");
        exit(1);
    }

    n = atoi(argv[1]);
    if (argc == 3) c = atoi(argv[2]);

    printf("v1: n=%d c=%d\n", n, c);

    /* Reserva de matriz alineada */
    double **a = aligned_alloc(TAM_LINEA_CACHE, sizeof(double*) * n);
    double *data = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n * n);

    for (int i = 0; i < n; i++)
        a[i] = &data[i * n];

    /* Reserva de vectores */
    double *b = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n);
    double *x = aligned_alloc(TAM_LINEA_CACHE, sizeof(double) * n);

    /* Inicializar x */
    for(int i = 0; i < n; i++)
        x[i] = 0.0;

    srand(time(NULL));

    generar_matriz(a);
    generar_vector(b);

    double ck;

    start_counter();
    double solucion = resolver_jacobi(a, b, x);
    ck = get_counter();

    printf("v1: solucion = %e\n", solucion);
    printf("v1: counter = %lf\n", ck);

    free(a);
    free(data);
    free(b);
    free(x);

    return 0;
}

double resolver_jacobi(double** a, double *b, double *x)
{
    double norm2, sigma;
    double x_new[n];

    for (int iter = 0; iter < max_iter; iter++)
    {
        norm2 = 0;

        for (int i = 0; i < n; i++)
        {
            sigma = 0.0;

            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    sigma += a[i][j] * x[j];
            }

            x_new[i] = (b[i] - sigma) / a[i][i];
            norm2 += pow(x_new[i] - x[i], 2);
        }

        /* copiar nuevo vector */
        for(int i = 0; i < n; i++)
            x[i] = x_new[i];

        if (sqrt(norm2) < tol)
            return norm2;
    }

    return -1;
}

/* Genera matriz diagonalmente dominante */
void generar_matriz(double** matriz)
{
    for (int i = 0; i < n; i++)
    {
        double suma = 0.0;
        for (int j = 0; j < n; j++)
        {
            matriz[i][j] = ((double) rand() / RAND_MAX) * 100.0;
            suma += matriz[i][j];
        }
        // La diagonal es la suma de toda la fila (garantiza dominancia estricta)
        matriz[i][i] += suma;
    }
}

/* Genera vector aleatorio */
void generar_vector(double* vector)
{
    for (int i = 0; i < n; i++)
        vector[i] = ((double) rand() / RAND_MAX) * 100.0;
}