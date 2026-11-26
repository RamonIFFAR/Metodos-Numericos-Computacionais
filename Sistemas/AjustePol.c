#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

/***********************************************************************************
******************************** Ajuste Polinomial *********************************
************************************************************************************/

void AP_constroi_vandermonde_simplificado(double* Xis, double* Y, int n, int m, double* S, double* T)
{
    for(int i = 0 ; i < n ; i++)
    {
        double potencia = 1.0;
        for(int p = 0 ; p <= 2*m ; p++)
        {
            S[p] += potencia;
            if(p <= m)
            {
                T[p] += Y[i]*potencia;
            }
            potencia *= Xis[i];
        }
    }
}

void AP_constroi_matrizA(int m, double* matriz, double* S)
{
    for(int i = 0 ; i < m ; i++)
    {
        for(int j = 0 ; j < m ; j++)
        {
            matriz[m*i + j] = S[i+j];
        }
    }
}

int AjustePolinomial(double* Xis, double* Y, int n, int m, double* resposta)
{
    if(Xis == NULL || Y == NULL || resposta == NULL || n < m+1 || m < 0)
    {
        return -1;
    }
    // Possivelmente retornar um int pra sinalizar erros
    double S[2*m+1], T[m+1];
    memset(S, 0, sizeof(double)*(2*m+1));
    memset(T, 0, sizeof(double)*(m+1));
    AP_constroi_vandermonde_simplificado(Xis, Y, n, m, S, T);
    double A[(m+1)*(m+1)];
    double b[m+1];
    AP_constroi_matrizA(m+1, A, S);
    memcpy(b, T, sizeof(double)*(m+1));

    printf("Matriz Construida:\n");
    printa_matriz(A, m+1);
    printf("\n\n");

    printf("Vetor Resposta:\n");
    printa_lista(b, m+1);
    printf("\n\n");

    printf("Matriz Triangulada:\n");
    int flag = transforma_matriz(A, m+1, b);
    if(flag == -1)
    {
        printf("Operacao Abortada!\nMatriz sob resolucao é singular. Insira outra matriz.");
        return -1;
    }
    printa_matriz(A, m+1);
    printf("\n\n");

    printf("Vetor Resposta Matriz:\n");
    printa_lista(b, m+1);
    printf("\n\n");

    flag = resolve_matriz(A, m+1, b, resposta);
    if(flag == -1)
    {
        printf("Operacao Abortada!\nMatriz sob resolucao é singular. Insira outra matriz.");
        return -1;
    }
    //printa_lista(resposta, m+1);

    return 1;
}

int main(void)
{
    int n, grau;
    printf("Quantidade de pontos: ");
    if(scanf("%d", &n) != 1 || n < 1 || n > 500)
    {
        fprintf(stderr, "Quantidade inválida (use de 1 a 500 pontos).\n");
        return EXIT_FAILURE;
    }
    printf("Grau do polinômio de ajuste: ");
    if(scanf("%d", &grau) != 1 || grau < 0 || grau >= n)
    {
        fprintf(stderr, "Grau inválido (deve ser não negativo e menor que o número de pontos).\n");
        return EXIT_FAILURE;
    }

    double* Xis = malloc(sizeof(double) * (size_t)n);
    double* Y = malloc(sizeof(double) * (size_t)n);
    double* coeficientes = malloc(sizeof(double) * (size_t)(grau + 1));
    if(Xis == NULL || Y == NULL || coeficientes == NULL)
    {
        fprintf(stderr, "Não foi possível alocar memória.\n");
        free(Xis);
        free(Y);
        free(coeficientes);
        return EXIT_FAILURE;
    }
    for(int i = 0; i < n; i++)
    {
        printf("x[%d] e y[%d]: ", i, i);
        if(scanf("%lf %lf", &Xis[i], &Y[i]) != 2)
        {
            fprintf(stderr, "Entrada inválida.\n");
            free(Xis);
            free(Y);
            free(coeficientes);
            return EXIT_FAILURE;
        }
    }

    if(AjustePolinomial(Xis, Y, n, grau, coeficientes) >= 0)
    {
        printf("Coeficientes do polinômio (termo constante primeiro):\n");
        printa_lista(coeficientes, grau + 1);
        printf("\n");
    }
    else
    {
        fprintf(stderr, "Não foi possível calcular o ajuste.\n");
    }
    free(Xis);
    free(Y);
    free(coeficientes);
    return EXIT_SUCCESS;
}
