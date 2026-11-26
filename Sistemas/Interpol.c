#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"

/***********************************************************************************
***************************** Interpolação Polinomial ******************************
************************************************************************************/

void InterpolPol_GeraMatX(double* MatX, double* Xis, int tamLinha)
{
    for(int i = 0 ; i < tamLinha ; i++)
    {
        for(int j = 1 ; j < tamLinha ; j++)
        {
            *(MatX + i*tamLinha + j) = pow(Xis[i], j);
        }
        *(MatX + i*tamLinha) = 1;
    }
}

int InterpolPol(double* Xis, double* FXis, int n, double* respostas)
{
    if(Xis == NULL || FXis == NULL || respostas == NULL || n <= 0)
    {
        return -1;
    }

    double* CopiaFXis = malloc(sizeof(double)*n);
    if(CopiaFXis == NULL)
    {
        return -1;
    }
    memcpy(CopiaFXis, FXis, sizeof(double)*n);
    double* MatX = malloc(sizeof(double)*n*n);
    if(MatX == NULL)
    {
        free(CopiaFXis);
        return -1;
    }

    int flag = busca_repetidos(Xis, n);

    if(flag == -1)
    {
        printf("Valores de X repetidos!\nAbortando operacao\n");
        free(MatX);
        free(CopiaFXis);
        return -1;
    }

    InterpolPol_GeraMatX(MatX, Xis, n);
    //printa_matriz(MatX, n);
    //printf("\n\n");

    flag = transforma_matriz(MatX, n, CopiaFXis);
    if(flag == -1)
    {
        printf("Operacao Abortada!\nMatriz sob transformacao é singular. Insira outra matriz.");
        free(MatX);
        free(CopiaFXis);
        return -1;
    }
    //printa_matriz(MatX, n);
    //printf("\n");
    //printa_lista(CopiaFXis, n);
    //printf("\n\n");

    flag = resolve_matriz(MatX, n, CopiaFXis, respostas);
    if(flag == -1)
    {
        printf("Operacao Abortada!\nMatriz sob resolucao é singular. Insira outra matriz.");
        free(MatX);
        free(CopiaFXis);
        return -1;
    }
    free(MatX);
    free(CopiaFXis);
    return 0;
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
    double* respostas = malloc(sizeof(double) * (size_t)n);
    if(Xis == NULL || FXis == NULL || respostas == NULL)
    {
        fprintf(stderr, "Não foi possível alocar memória.\n");
        free(Xis);
        free(FXis);
        free(respostas);
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
            free(respostas);
            return EXIT_FAILURE;
        }
    }

    if(InterpolPol(Xis, FXis, n, respostas) == 0)
    {
        printf("Coeficientes do polinômio (termo constante primeiro):\n");
        printa_lista(respostas, n);
        printf("\n");
    }
    else
    {
        fprintf(stderr, "Não foi possível calcular a interpolação.\n");
    }

    free(Xis);
    free(FXis);
    free(respostas);
    return EXIT_SUCCESS;
}