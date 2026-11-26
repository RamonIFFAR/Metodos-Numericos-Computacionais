#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

/***********************************************************************************
************************************* Newton ***************************************
************************************************************************************/

void Newton_defineCoeficientes(double* Xis, double* FXis, int n, double coeficientes[][n])
{
    for(int i = 0 ; i < n ; i++)
    {
        coeficientes[i][0] = FXis[i];
    }

    // Sendo n a quantidade de pontos, Xis o vetor contendo os pontos X e coeficientes a matriz contendo as diferenças divididas
    for(int i = 1 ; i < n ; i++) // Percorre as colunas
    {
        for(int j = 0 ; j < n - i ; j++) // Percorre as linhas
        {
            coeficientes[j][i] = (coeficientes[j+1][i-1] - coeficientes[j][i-1]) / (Xis[j+i] - Xis[j]);
        }
    }
}

void Newton_printaFuncao(double* Xis, int n, double coeficientes[][n])
{
    for(int i = 0 ; i < n ; i++)
    {
        printf("%f", coeficientes[0][i]);
        for(int j = 0 ; j < i ; j++)
        {
            printf("(x - %f)", Xis[j]);
        }
        if(i != n-1)
        {
            printf(" + ");
        }
    }
}

void Newton(double* Xis, double* FXis, int n)
{
    double Diferencas[n][n]; // Referente aos coeficientes das diferenças divididas
    Newton_defineCoeficientes(Xis, FXis, n, Diferencas);
    Newton_printaFuncao(Xis, n, Diferencas);
    //printa_matriz_matricial(n, Diferencas);
}

int main(void)
{
    int n;
    printf("Quantidade de pontos: ");
    if(scanf("%d", &n) != 1 || n < 1 || n > 500)
    {
        fprintf(stderr, "Quantidade inválida (use de 1 a 500 pontos).\n");
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

    printf("Polinômio interpolador na forma de Newton:\n");
    Newton(Xis, FXis, n);
    printf("\n");
    free(Xis);
    free(FXis);
    return EXIT_SUCCESS;
}
