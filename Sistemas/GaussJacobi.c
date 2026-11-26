#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"

/***********************************************************************************
******************************** Método Gauss-Jacob ********************************
************************************************************************************/

int GaussJacob_CalculaVetor(double* MatrizCoeficientes, int numLinhas, double* vetAtual, double* vetImagem, int tam)
{

    int i, j;
    double vetTemp[tam];
    memcpy(vetTemp, vetAtual, sizeof(double)*tam);
    for(i = 0 ; i < numLinhas ; i++)
    {
        double temp = vetImagem[i];
        for(j = 0 ; j < tam;  j++)
        {
            if(j == i) continue;
            temp -= vetTemp[j] * (*(MatrizCoeficientes + i * tam + j));
        }
        if(fabs(*(MatrizCoeficientes + i * tam + i)) < 1e-12)
        {
            return -1;
        }
        temp /= (*(MatrizCoeficientes + i * tam + i));
        //printf("Temp[%d]: %f\n", i, temp);
        vetAtual[i] = temp;
    }
    return 0;
}

int GaussJacob(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia, int max_iteracoes, double* respostas)
{
    double anterior[tam];
    memset(anterior, 0, sizeof(double)*tam);
    for(int j = 0 ; j < tam ; j++)
    {
        if(pivotiza(MatrizCoeficientes, FXis, tam, j, NULL) != 0)
        {
            return -1;
        }
    }
    int i = 0;
    do
    { 
        memcpy(anterior, respostas, sizeof(double)*tam);
        int flag = GaussJacob_CalculaVetor(MatrizCoeficientes, numLinhas, respostas, FXis, tam);
        if(flag == -1)
        {
            printf("Erro!\n Programa detectou matriz indeterminada! Encerrando método...\n");
            return -1;
        }
        i++;
        //printa_lista(anterior, tam);
        //printf("\n");
        //printa_lista(respostas, tam);
        //printf("\n");
        //printf("deu em: %d\n", GaussJacob_CalculaErro(respostas, anterior, tam) <= tolerancia);
    } while(GaussX_CalculaErro(respostas, anterior, tam) > tolerancia && max_iteracoes > i);
    return 0;
}

int main(void)
{
    int n, max_iteracoes;
    double tolerancia;
    printf("Ordem do sistema: ");
    if(scanf("%d", &n) != 1 || n < 1 || n > 500)
    {
        fprintf(stderr, "Ordem inválida (use de 1 a 500).\n");
        return EXIT_FAILURE;
    }

    double* matriz = malloc(sizeof(double) * (size_t)n * (size_t)n);
    double* b = malloc(sizeof(double) * (size_t)n);
    double* respostas = calloc((size_t)n, sizeof(double));
    if(matriz == NULL || b == NULL || respostas == NULL)
    {
        fprintf(stderr, "Não foi possível alocar memória.\n");
        free(matriz);
        free(b);
        free(respostas);
        return EXIT_FAILURE;
    }

    printf("Informe a matriz A por linhas (%d valores):\n", n * n);
    for(int i = 0; i < n * n; i++)
    {
        if(scanf("%lf", &matriz[i]) != 1)
        {
            fprintf(stderr, "Entrada inválida para a matriz.\n");
            free(matriz);
            free(b);
            free(respostas);
            return EXIT_FAILURE;
        }
    }
    printf("Informe os %d valores do vetor b:\n", n);
    for(int i = 0; i < n; i++)
    {
        if(scanf("%lf", &b[i]) != 1)
        {
            fprintf(stderr, "Entrada inválida para o vetor b.\n");
            free(matriz);
            free(b);
            free(respostas);
            return EXIT_FAILURE;
        }
    }
    printf("Tolerância e número máximo de iterações: ");
    if(scanf("%lf %d", &tolerancia, &max_iteracoes) != 2 || tolerancia <= 0.0 || max_iteracoes < 1)
    {
        fprintf(stderr, "Tolerância ou limite de iterações inválido.\n");
        free(matriz);
        free(b);
        free(respostas);
        return EXIT_FAILURE;
    }

    if(GaussJacob(b, n, matriz, n, tolerancia, max_iteracoes, respostas) == 0)
    {
        printf("Aproximação da solução:\n");
        printa_lista(respostas, n);
        printf("\n");
    }
    else
    {
        fprintf(stderr, "Falha ao calcular a solução.\n");
    }
    free(matriz);
    free(b);
    free(respostas);
    return EXIT_SUCCESS;
}
