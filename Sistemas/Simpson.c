#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

double Simpson(double* Xis, double* FXis, int n)
{
    if(n < 3)
    {
        printf("Número de pontos Insuficiente! Precisa-se de pelo menos 3 pontos\n");
        return -1;
    }

    double step = Xis[1] - Xis[0];
    double sum = 0.0;
    int i = 0;
    while(i+2 < n)
    {
        if(modulo((Xis[i+1] - Xis[i]) - step) > 1e-12 || modulo((Xis[i+2] - Xis[i+1]) - step) > 1e-12)
        {
            printf("Pontos inseridos não estao igualmente espaçados!\nOperação cancelada");
            return 0.0;
        }
        sum += (FXis[i] + FXis[i+1] * 4 + FXis[i+2]) * step / 3;
        i += 2;
    }
    if(i != n-1)
    {
        sum += (FXis[i] + FXis[i+1]) / 2 * (Xis[i+1] - Xis[i]);
    }
    return sum;
}

int main(void)
{
    int n;
    printf("Quantidade de pontos: ");
    if(scanf("%d", &n) != 1 || n < 3 || n > 10000)
    {
        fprintf(stderr, "Quantidade inválida (use de 3 a 10000 pontos).\n");
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

    printf("Integral aproximada: %.12g\n", Simpson(Xis, FXis, n));
    free(Xis);
    free(FXis);
    return EXIT_SUCCESS;
}
