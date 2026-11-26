#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tools.h"

#define VARIACAO_TOLERADA 1e-12

double modulo(double x)
{
    if(x < 0)
    {
        return x * -1;
    }
    return x;
}

int busca_repetidos(double* Xis, int tam)
{
    for(int i = 0 ; i < tam-1 ; i++)
    {
        for(int j = i+1 ; j < tam ; j++)
        {
            if(modulo(Xis[i] - Xis[j]) < VARIACAO_TOLERADA)
            {
                return -1;
            }
        }
    }
    return 0;
}

double GaussX_CalculaErro(double* atual, double* anterior, int tam)
{
    double modulo_temp = 0;
    double modulo_atual = 0;

    for(int i = 0 ; i < tam ; i++)
    {
        double temp = atual[i] - anterior[i];
        modulo_temp += temp * temp;
        modulo_atual += atual[i] * atual[i];
    }
    modulo_temp = sqrt(modulo_temp);
    modulo_atual = sqrt(modulo_atual);

    if(modulo_atual < VARIACAO_TOLERADA)
    {
        return modulo_temp;
    }
    return modulo_temp / modulo_atual;
}

static void troca_linhas(double* A, double* b, int n, int l1, int l2)
{
    if (l1 == l2) return;

    for (int j = 0; j < n; j++)
    {
        double tmp = A[l1*n + j];
        A[l1*n + j] = A[l2*n + j];
        A[l2*n + j] = tmp;
    }
    if (b != NULL)
    {
        double tmp = b[l1];
        b[l1] = b[l2];
        b[l2] = tmp;
    }
}

int pivotiza(double* A, double* b, int n, int coluna, int* perm)
{
    if (A == NULL || n <= 0 || coluna < 0 || coluna >= n)
        return -2;

    int linha_maior = coluna;
    double maior = modulo(A[coluna*n + coluna]);

    for (int i = coluna + 1; i < n; i++)
    {
        double atual = modulo(A[i*n + coluna]);
        if (atual > maior)
        {
            maior = atual;
            linha_maior = i;
        }
    }

    if (maior < VARIACAO_TOLERADA)
        return -1;

    troca_linhas(A, b, n, coluna, linha_maior);

    if (perm != NULL && linha_maior != coluna)
    {
        int tmp = perm[coluna];
        perm[coluna] = perm[linha_maior];
        perm[linha_maior] = tmp;
    }
    return 0;
}

void diagonaliza(double* mat, int tam)
{
    for(int i = 0 ; i < tam ; i++)
    {
        *(mat + i*tam + i) = 1;
    }
}

int transforma_matriz_simple(double* matriz, int tam, double* MatrizL)
{
    double multiplicador = 0;
    for( int linha_ref = 0 ; linha_ref < tam ; linha_ref++)
    {
        if(modulo(*(matriz + linha_ref * tam + linha_ref)) < VARIACAO_TOLERADA)
        {
            return -1;
        }
        for( int linha_atual = linha_ref+1 ; linha_atual < tam ; linha_atual++) // Passa pelas linhas pra serem alteradas
        {
            multiplicador = *(matriz + linha_atual * tam + linha_ref) / *(matriz + linha_ref * tam + linha_ref);
            multiplicador *= -1;
            for(int i = linha_ref ; i < tam ; i++) // Transforma a linha
            {
                *(matriz + linha_atual * tam + i) += *(matriz + linha_ref * tam + i) * multiplicador;
            }
            *(MatrizL + tam * linha_atual + linha_ref) = multiplicador * -1;
        }
    }
    return 0;
}

int transforma_matriz(double* matriz, int m, double* decoy)
{
    double multiplicador = 0;
    for( int linha_ref = 0 ; linha_ref < m - 1; linha_ref++)
    {
        int flag = pivotiza(matriz, decoy, m, linha_ref, NULL);
        if(flag != 0)
        {
            return -1;
        }
        for( int linha_atual = linha_ref+1 ; linha_atual < m ; linha_atual++) // Passa pelas linhas pra serem alteradas
        {
            multiplicador = *(matriz + linha_atual * m + linha_ref) / *(matriz + linha_ref * m + linha_ref);
            multiplicador *= -1;
            for(int i = linha_ref ; i < m ; i++) // Transforma a linha
            {
                *(matriz + linha_atual * m + i) += *(matriz + linha_ref * m + i) * multiplicador;
            }
            decoy[linha_atual] += decoy[linha_ref] * multiplicador;
            //printa_lista(decoy, m);
            //printf("\n");
        }
    }
    return 0;
}

int resolve_matriz(double* matriz, int tam, double* decoy, double* resultado)
{
    // Incluir checagem da matriz ser resolvível ou não *insira emoji de caveira agui
    for(int i = tam-1 ; i >= 0 ; i--)
    {
        int j = tam-1;
        double temp = decoy[i];

        for(; j > i ; j--)
        {
            temp -= *(matriz + i * tam + j) * resultado[j];
        }
        double flag = *(matriz + i * tam + j);
        
        if(modulo(flag) < VARIACAO_TOLERADA)
        {
            return -1;
        }

        resultado[j] = temp / *(matriz + i * tam + j);
    }
    return 0;
}

void printa_matriz_matricial(int n, double matriz[][n])
{
    for(int i = 0 ; i < n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            printf("%f ", matriz[i][j]);
        }
        printf("\n");
    }
}

void printa_matriz(double* matriz, int tam)
{
    for(int i = 0 ; i < tam ; i++)
    {
        for(int j = 0 ; j < tam ; j++)
        {
            printf("%f ", *(matriz + i * tam + j));
        }
        printf("\n");
    }
}

void printa_lista(double* lista, int tam)
{
    for(int i = 0; i < tam ; i++ )
    {
        printf("%f ", lista[i]);
    }
}