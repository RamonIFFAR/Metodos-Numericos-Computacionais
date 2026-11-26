#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

void GaussSeidel_CalculaVetor(double* MatrizCoeficientes, int numLinhas, double* vetAtual, double* vetImagem, int tam)
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
        temp /= (*(MatrizCoeficientes + i * tam + i));
        vetTemp[i] = temp;
        vetAtual[i] = temp;
        printa_lista(vetTemp, numLinhas);
        printf("\n");
    }
}

double* GaussSeidel(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia, int max_iteracoes)
{
    double* respostas = calloc(tam, sizeof(double));
    double anterior[tam];
    memset(anterior, 0, sizeof(double)*tam);
    int i = 0;
    do
    {
        memcpy(anterior, respostas, sizeof(double)*tam);
        GaussSeidel_CalculaVetor(MatrizCoeficientes, numLinhas, respostas, FXis, tam);
        i++;
        //printa_lista(anterior, tam);
        //printf("\n");
        //printa_lista(respostas, tam);
        //printf("\n");
        //printf("deu em: %d\n", GaussJacob_CalculaErro(respostas, anterior, tam) <= tolerancia);
    } while(GaussX_CalculaErro(respostas, anterior, tam) > tolerancia && max_iteracoes > i);
    return respostas;
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
    if(matriz == NULL || b == NULL)
    {
        fprintf(stderr, "Não foi possível alocar memória.\n");
        free(matriz);
        free(b);
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
            return EXIT_FAILURE;
        }
    }
    printf("Informe os %d valores do vetor b:\n", n);
    for(int i = 0; i < n; i++)
    {
        if(scanf("%lf", &b[i]) != 1 || modulo(matriz[i * n + i]) < 1e-12)
        {
            fprintf(stderr, "Vetor inválido ou elemento diagonal nulo na linha %d.\n", i + 1);
            free(matriz);
            free(b);
            return EXIT_FAILURE;
        }
    }
    printf("Tolerância e número máximo de iterações: ");
    if(scanf("%lf %d", &tolerancia, &max_iteracoes) != 2 || tolerancia <= 0.0 || max_iteracoes < 1)
    {
        fprintf(stderr, "Tolerância ou limite de iterações inválido.\n");
        free(matriz);
        free(b);
        return EXIT_FAILURE;
    }

    double* respostas = GaussSeidel(b, n, matriz, n, tolerancia, max_iteracoes);
    if(respostas == NULL)
    {
        fprintf(stderr, "Falha ao calcular a solução.\n");
        free(matriz);
        free(b);
        return EXIT_FAILURE;
    }
    printf("Aproximação da solução:\n");
    printa_lista(respostas, n);
    printf("\n");
    free(respostas);
    free(matriz);
    free(b);
    return EXIT_SUCCESS;
}
