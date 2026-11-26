#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"

double somatorio_expoente(double* num, int tam, int potencia)
{
    double resultado = 0;
    for(int i = 0 ; i < tam; i++)
    {
        resultado += pow(num[i], potencia);
    }
    return resultado;
}

void constroi_matriz_resultado(double* FXis, double* Xis, double* decoy, int n, int grau)
{
    int m = grau+1;
    
    for(int i = 0; i < m; i++)
    {
        int grau_decresce = grau - i;
        for(int j = 0 ; j < n ; j++)
        {
            decoy[i] += FXis[j] * pow(Xis[j],  grau_decresce);
        }
    }
}

void constroi_matriz(double* matriz, double* Xis, int n, int grau)
{
    int m = grau+1;

    for(int i = 0 ; i < m; i++)
    {
        
        for(int j = 0 ; j < m ; j++)
        {
            int grau_decresce = 2 * grau - i - j;
            matriz[i*m+j] = somatorio_expoente(Xis, n, grau_decresce);
        }
    
    }
}

/*
* Xis e FXis. int tam se refere ao valor de "n".
* Retorna um ponteiro de int com os valores correspondetes à a (idx = 0) até m (idx = grau)
*/
double* Gauss(double* Xis, double* FXis, int n, int grau) 
{
    int m = grau+1;
    double matriz[m*m];
    double decoy[m];

    int flag = busca_repetidos(Xis, n);

    if(flag == -1)
    {
        printf("Valores de X repetidos!\nAbortando operacao\n");
        return NULL;
    }

    if(n < grau+1)
    {
        printf("Quantidade Insuficiente de Pontos!");
        return NULL;
    }

    memset(decoy, 0, sizeof(double) * m);

    constroi_matriz_resultado(FXis, Xis, decoy, n, grau);
    printf("Saída Decoy:\n");
    printa_lista(decoy, m);
    printf("\n\n");

    constroi_matriz(matriz, Xis, n, grau);
    printf("Saída Matriz Construída:\n");
    printa_matriz(matriz, m);
    printf("\n\n");

    printf("Saída Matriz Transformada:\n");
    flag = transforma_matriz(matriz, m, decoy);

    if(flag == -1)
    {
        printf("Operacao Abortada!\nMatriz sob transformação é singular. Insira outra matriz.");
        return NULL;
    }

    printa_matriz(matriz, m);
    printf("\n\n");

    printf("Saída Decoy Final:\n");
    printa_lista(decoy, m);
    printf("\n\n");

    
    double* resposta = malloc(sizeof(double) * m);

    flag = resolve_matriz(matriz, m, decoy, resposta);
    if(flag == -1)
    {
        printf("Operacao Abortada!\nMatriz sob resolucao é singular. Insira outra matriz.");
        return NULL;
    }

    printf("\nSaída Resposta:\n");
    printa_lista(resposta, m);
    return resposta;
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

    double* coeficientes = Gauss(Xis, FXis, n, grau);
    if(coeficientes != NULL)
    {
        printf("Coeficientes retornados:\n");
        printa_lista(coeficientes, grau + 1);
        printf("\n");
        free(coeficientes);
    }
    free(Xis);
    free(FXis);
    return EXIT_SUCCESS;
}
