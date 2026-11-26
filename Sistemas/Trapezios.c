#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

double Trapezios(double* Xis, double* FXis, int n)
{
    double area = 0.0;
    for(int i = 0 ; i < n-1 ; i++)
    {
        area += (FXis[i] + FXis[i+1]) / 2 * (Xis[i+1] - Xis[i]);
    }
    return area;
}

int main(void)
{
    int n;
    printf("Quantidade de pontos: ");
    if(scanf("%d", &n) != 1 || n < 2 || n > 10000)
    {
        fprintf(stderr, "Quantidade inválida (use de 2 a 10000 pontos).\n");
        return EXIT_FAILURE;
    }

    double* Xis = malloc(sizeof(double) * (size_t)n);
    double* FXis = malloc(sizeof(double) * (size_t)n);
    if(Xis == NULL || FXis == NULL)
    {
        fprintf(stderr, "Não foi possível alocar memória.\n");
        free(Xis);
        free(FXis);
        return EXIT_FAILURE;
    }

    for(int i = 0; i < n; i++)
    {
        printf("x[%d] e f(x[%d]): ", i, i);
        if(scanf("%lf %lf", &Xis[i], &FXis[i]) != 2)
        {
            fprintf(stderr, "Entrada inválida.\n");
            free(Xis);
            free(FXis);
            return EXIT_FAILURE;
        }
    }

    printf("Integral aproximada: %.12g\n", Trapezios(Xis, FXis, n));
    free(Xis);
    free(FXis);
    return EXIT_SUCCESS;
}
