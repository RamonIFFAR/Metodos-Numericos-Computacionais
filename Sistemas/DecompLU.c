#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

/***********************************************************************************
******************************** Método Decomp. LU *********************************
************************************************************************************/

void DecompLU_resolveL(double* MatL, double* MatY, double* MatB, int tamLinha)
{
    double CopiaB[tamLinha];
    memcpy(CopiaB, MatB, sizeof(double)*tamLinha);
    for(int i = 0 ; i < tamLinha ; i++)
    {
        int j;
        for(j = i-1; j >= 0 ; j--)
        {
            CopiaB[i] -= *(MatL + tamLinha * i + j) * MatY[j];
        }
        MatY[i] = CopiaB[i] / *(MatL + tamLinha * i + i);
    }
}

void DecompLU_resolveU(double* MatU, double* MatX, double* MatY, int tamLinha)
{
    double CopiaY[tamLinha];
    memcpy(CopiaY, MatY, sizeof(double)*tamLinha);
    for(int i = tamLinha-1 ; i >= 0 ; i--)
    {
        int j;
        for(j = i+1; j < tamLinha ; j++)
        {
            CopiaY[i] -= *(MatU + tamLinha * i + j) * MatX[j];
        }
        MatX[i] = CopiaY[i] / *(MatU + tamLinha * i + i);
    }
}

int DecompLU(double* MatrizA, double* VetX, double* VetB, double* MatrizL, double* MatrizU, int tamLinha)
{
    /*  Ax = B
    *   A = LU
    *   LUx = B
    *   Ux = Y
    *   LY = B
    */
    double CopiaA[tamLinha*tamLinha];
    double VetY[tamLinha];

    memset(VetY, 0, sizeof(double)*tamLinha);
    memcpy(CopiaA, MatrizA, sizeof(double)*tamLinha*tamLinha);

    memset(MatrizL, 0, sizeof(double)*tamLinha*tamLinha);
    diagonaliza(MatrizL, tamLinha);

    int flag = transforma_matriz_simple(CopiaA, tamLinha, MatrizL);
    if(flag == -1)
    {
        if(flag == -1)
        {
            printf("Operacao Abortada!\nMatriz sob resolucao é singular. Insira outra matriz.");
            return -1;
        }
    }
    memcpy(MatrizU, CopiaA, sizeof(double)*tamLinha*tamLinha);
    DecompLU_resolveL(MatrizL, VetY, VetB, tamLinha);
    DecompLU_resolveU(MatrizU, VetX, VetY, tamLinha);
    return 0;
    /*
    printa_matriz(CopiaA, tamLinha);
    printf("\n");
    printa_matriz(MatrizL, tamLinha);
    printf("\n");
    printa_lista(VetY, tamLinha);
    printf("\n");
    printa_lista(VetX, tamLinha);
    */
}

int main(void)
{
    int n;
    printf("Ordem da matriz A: ");
    if(scanf("%d", &n) != 1 || n < 1 || n > 500)
    {
        fprintf(stderr, "Ordem inválida (use de 1 a 500).\n");
        return EXIT_FAILURE;
    }

    size_t matrix_size = (size_t)n * (size_t)n;
    double* A = malloc(sizeof(double) * matrix_size);
    double* L = malloc(sizeof(double) * matrix_size);
    double* U = malloc(sizeof(double) * matrix_size);
    double* b = malloc(sizeof(double) * (size_t)n);
    double* x = malloc(sizeof(double) * (size_t)n);
    if(A == NULL || L == NULL || U == NULL || b == NULL || x == NULL)
    {
        fprintf(stderr, "Não foi possível alocar memória.\n");
        free(A);
        free(L);
        free(U);
        free(b);
        free(x);
        return EXIT_FAILURE;
    }

    printf("Informe a matriz A por linhas (%d valores):\n", n * n);
    for(size_t i = 0; i < matrix_size; i++)
    {
        if(scanf("%lf", &A[i]) != 1)
        {
            fprintf(stderr, "Entrada inválida para a matriz.\n");
            free(A); free(L); free(U); free(b); free(x);
            return EXIT_FAILURE;
        }
    }
    printf("Informe os %d valores do vetor b:\n", n);
    for(int i = 0; i < n; i++)
    {
        if(scanf("%lf", &b[i]) != 1)
        {
            fprintf(stderr, "Entrada inválida para o vetor b.\n");
            free(A); free(L); free(U); free(b); free(x);
            return EXIT_FAILURE;
        }
    }

    if(DecompLU(A, x, b, L, U, n) == 0)
    {
        printf("Matriz L:\n");
        printa_matriz(L, n);
        printf("Matriz U:\n");
        printa_matriz(U, n);
        printf("Solução x:\n");
        printa_lista(x, n);
        printf("\n");
    }
    else
    {
        fprintf(stderr, "Não foi possível decompor/resolver o sistema.\n");
    }
    free(A); free(L); free(U); free(b); free(x);
    return EXIT_SUCCESS;
}
