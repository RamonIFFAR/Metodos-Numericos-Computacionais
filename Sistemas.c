#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
*       FERRAMENTAS
*/

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

void diagonaliza(double* mat, int tam)
{
    for(int i = 0 ; i < tam ; i++)
    {
        *(mat + i*tam + i) = 1;
    }
}

/***********************************************************************************
**************************** Método Eliminação de Gauss ****************************
************************************************************************************/

// Recebe os valores a serem somados em num, a quantidade de números em tam e a potencia à qual serão elevados em potencia 
double Pow(double base, int expoente)
{
    double resultado = 1;
    for(int i = 0 ; i < expoente ; i++)
    {
        resultado *= base;
    }
    return resultado;
}

double somatorio_expoente(double* num, double tam, double potencia)
{
    double resultado = 0;
    for(int i = 0 ; i < tam; i++)
    {
        resultado += Pow(num[i], potencia);
    }
    return resultado;
}

void constroi_matriz_resultado(double* FXis, double* Xis, double* decoy, int tam, int grau)
{
    int grau_decresce = grau;
    for(int i = 0; i < tam; i++)
    {
        for(int j = 0 ; j < tam ; j++)
        {
            decoy[i] += FXis[j] * Pow(Xis[j],  grau_decresce);
        }
        grau_decresce--;
    }
}

void constroi_matriz(double* matriz, double* Xis, int tam, int grau, double* decoy)
{
    int grau_max = 2*grau, grau_decresce;
    for(int i = 0 ; i < tam; i++)
    {
        grau_decresce = grau_max;
        for(int j = 0 ; j < tam ; j++)
        {
            matriz[i*tam+j] = somatorio_expoente(Xis, tam, grau_decresce);
            grau_decresce--;
        }
        grau_max--;
    }
}

void transforma_matriz_simple(double* matriz, int tam, double* MatrizL)
{
    double multiplicador = 0;
    for( int linha_ref = 0 ; linha_ref < tam ; linha_ref++)
    {
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
}

void transforma_matriz(double* matriz, int tam, double* decoy)
{
    double multiplicador = 0;
    for( int linha_ref = 0 ; linha_ref < tam ; linha_ref++)
    {
        for( int linha_atual = linha_ref+1 ; linha_atual < tam ; linha_atual++) // Passa pelas linhas pra serem alteradas
        {
            multiplicador = *(matriz + linha_atual * tam + linha_ref) / *(matriz + linha_ref * tam + linha_ref);
            multiplicador *= -1;
            for(int i = linha_ref ; i < tam ; i++) // Transforma a linha
            {
                *(matriz + linha_atual * tam + i) += *(matriz + linha_ref * tam + i) * multiplicador;
            }
            decoy[linha_atual] += decoy[linha_ref] * multiplicador;
            printa_lista(decoy, tam);
            printf("\n");
        }
    }
    
}

double* resolve_matriz(double* matriz, int tam, double* decoy)
{
    // Incluir checagem da matriz ser resolvível ou não *insira emoji de caveira agui
    double* variaveis = malloc(sizeof(double)*tam);
    for(int i = tam-1 ; i >= 0 ; i--)
    {
        int j = tam-1;
        for(j; j > i ; j--)
        {
            decoy[i] -= *(matriz + i * tam + j) * variaveis[j];
        }
        variaveis[j] = decoy[i] / *(matriz + i * tam + j);
    }
    return variaveis;
}

/*
* Xis e FXis. int tam se refere ao valor de "n".
* Retorna um ponteiro de int com os valores correspondetes à a (idx = 0) até m (idx = grau)
*/
double* Gauss(double* Xis, double* FXis, int tam, int grau) 
{
    double matriz[tam*tam];
    double decoy[tam];

    memset(decoy, 0, sizeof(double) * tam);

    constroi_matriz_resultado(FXis, Xis, decoy, tam, grau);
    printa_lista(decoy, tam);

    constroi_matriz(matriz, Xis, tam, grau, decoy);
    //printa_matriz(matriz, tam);
    //printf("\n");

    transforma_matriz(matriz, tam, decoy);
    //printa_matriz(matriz, tam);
    //printf("\n");

    //printa_lista(decoy, tam);
    //printf("\n");

    double* resposta = resolve_matriz(matriz, tam, decoy);
    //printf("\n");
    //printa_lista(resposta, tam);
    return resposta;
}

/***********************************************************************************
******************************** Método Gauss-Jacob ********************************
************************************************************************************/

double GaussX_CalculaErro(double* atual, double* anterior, int tam) // Recebe os dois vetores na iteracao k+1 e k, respectivamente
{
    double temp[tam];

    for(int i = 0 ; i < tam ; i++)
    {
        temp[i] = atual[i] - anterior[i];
    }

    double modulo_temp = 0;
    for(int i = 0 ; i < tam ; i++)
    {
        modulo_temp += temp[i]*temp[i];
    }
    modulo_temp = sqrt(modulo_temp);
    
    double modulo_atual = 0;
    for(int i = 0 ; i < tam ; i++)
    {
        modulo_atual += atual[i]*atual[i];
    }
    modulo_atual = sqrt(modulo_atual);

    return modulo_temp / modulo_atual;
}

void GaussJacob_CalculaVetor(double* MatrizCoeficientes, int numLinhas, double* vetAtual, double* vetAnterior, double* vetImagem, int tam)
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
            temp -= vetTemp[i] * (*(MatrizCoeficientes + i * tam + j));
        }
        temp /= (*(MatrizCoeficientes + i * tam + i));
        //printf("Temp[%d]: %f\n", i, temp);
        vetAtual[i] = temp;
    }
}

double* GaussJacob(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia)
{
    double* respostas = calloc(tam, sizeof(double));
    double anterior[tam] = {};
    do
    { 
        memcpy(anterior, respostas, sizeof(double)*tam);
        GaussJacob_CalculaVetor(MatrizCoeficientes, numLinhas, respostas, anterior, FXis, tam);
        //printa_lista(anterior, tam);
        //printf("\n");
        //printa_lista(respostas, tam);
        //printf("\n");
        //printf("deu em: %d\n", GaussJacob_CalculaErro(respostas, anterior, tam) <= tolerancia);
    } while(GaussX_CalculaErro(respostas, anterior, tam) > tolerancia);
    return respostas;
}

/***********************************************************************************
******************************** Método Gauss-Seidel ********************************
************************************************************************************/

void GaussSeidel_CalculaVetor(double* MatrizCoeficientes, int numLinhas, double* vetAtual, double* vetAnterior, double* vetImagem, int tam)
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

double* GaussSeidel(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia)
{
    double* respostas = calloc(tam, sizeof(double));
    double anterior[tam] = {};
    do
    { 
        memcpy(anterior, respostas, sizeof(double)*tam);
        GaussSeidel_CalculaVetor(MatrizCoeficientes, numLinhas, respostas, anterior, FXis, tam);
        //printa_lista(anterior, tam);
        //printf("\n");
        //printa_lista(respostas, tam);
        //printf("\n");
        //printf("deu em: %d\n", GaussJacob_CalculaErro(respostas, anterior, tam) <= tolerancia);
    } while(GaussX_CalculaErro(respostas, anterior, tam) > tolerancia);
    return respostas;
}

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

void DecompLU(double* MatrizA, double* VetX, double* VetB, double* MatrizL, double* MatrizU, int tamLinha)
{
    /*  Ax = B
    *   A = LU
    *   LUx = B
    *   Ux = Y
    *   LY = B
    */
    double CopiaA[tamLinha*tamLinha];
    double VetY[tamLinha] = {};
    memcpy(CopiaA, MatrizA, sizeof(double)*tamLinha*tamLinha);

    diagonaliza(MatrizL, tamLinha);

    transforma_matriz_simple(CopiaA, tamLinha, MatrizL);
    memcpy(MatrizU, CopiaA, sizeof(double)*tamLinha*tamLinha);
    DecompLU_resolveL(MatrizL, VetY, VetB, tamLinha);
    DecompLU_resolveU(MatrizU, VetX, VetY, tamLinha);
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

/***********************************************************************************
******************************** Ajuste Polinomial *********************************
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

double* InterpolPol(double* Xis, double* FXis, int tam)
{
    int n = tam;
    double* respostas;
    double* CopiaFXis = malloc(sizeof(double)*n);
    memcpy(CopiaFXis, FXis, sizeof(double)*n);
    double* MatX = malloc(sizeof(double)*n*n);

    InterpolPol_GeraMatX(MatX, Xis, n);
    //printa_matriz(MatX, n);
    //printf("\n\n");

    transforma_matriz(MatX, n, CopiaFXis);
    //printa_matriz(MatX, n);
    //printf("\n");
    //printa_lista(CopiaFXis, n);
    //printf("\n\n");

    respostas = resolve_matriz(MatX, n, CopiaFXis);
    free(MatX);
    free(CopiaFXis);
    return respostas;
}

/***********************************************************************************
************************************ Lagrange **************************************
************************************************************************************/

void Lagrange_multiplicacaoPermutada(double* coeficientesPos, double* Xis, int n)
{

}

void Lagrange_constroi(double* Xis, double* FXis, int n)
{
    
}

void Lagrange(double* Xis, double* FXis, int n)
{

}

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

/***********************************************************************************
************************************ Trapézios *************************************
************************************************************************************/

double Trapezios(double* Xis, double* FXis, int n)
{
    double area = 0.0;
    for(int i = 0 ; i < n-1 ; i++)
    {
        area += (FXis[i] + FXis[i+1]) / 2 * (Xis[i+1] - Xis[i]);
    }
    return area;
}

/***********************************************************************************
************************************* Simpson **************************************
************************************************************************************/

double Simpson(double* Xis, double* FXis, int n)
{
    double step = Xis[1] - Xis[0];
    double sum = 0.0;
    int i = 0;
    while(i+2 < n)
    {
        sum += (FXis[i] + FXis[i+1] * 4 + FXis[i+2]) * step / 3;
        i += 2;
    }
    if(i != n-1)
    {
        sum += Trapezios(Xis+i, FXis+i, 2);
    }
    return sum;
}

int main(int argc, char** argv)
{
    int grau = 2;
    double* resposta;
    double Xis[3] = {2, 4, 6};
    double FXis[3] = {1, 2, 3};
    double MatrizCoeficientes[9] = {10, 1, 1, 1, 10, 1, 1, 1, 10};

    /* Gauss */
    /**
    resposta = Gauss(Xis, FXis, 3, 2);
    printf("\n\n");
    /**/

    /* Gauss Jacob */
    /**
    resposta = GaussJacob(FXis, 3, MatrizCoeficientes, 3, 0.01);
    printa_lista(resposta, 3);
    printf("\n\n");
    /**/

    /* Gauss Seidel*/
    /**
    resposta = GaussSeidel(FXis, 3, MatrizCoeficientes, 3, 0.01);
    printf("\n\n");
    /**/

    /* Decomp LU */
    /**
    double MatTeste [9] = {3,2,4,1,1,2,4,3,2};
    double VetX[3] = {};
    double MatB[3] = {1,2,3};
    double MatL[9] = {};
    double MatU[9] = {};
    int tam = 3;
    DecompLU(MatTeste, VetX, MatB, MatL, MatU, tam);
    printa_lista(VetX, tam);
    /**/

    //resposta = InterpolPol(Xis, FXis, 3);
    //printa_lista(resposta, 3);

    /* Newton */
    /**
    double pontosX[3] = {1,2,3};
    double pontosF[3] = {2,5,10};
    Newton(pontosX, pontosF, 3);
    /**/

    /* Trapézios */
    /**
    int tam = 5;
    double TPontos[5] = {0,1,2,3,4};
    double TRes[5] = {3,5,4,8,6};
    printf("Resultado: %f", Trapezios(TPontos, TRes, tam));
    /**/

    /* Simpson */
    /**/
    double lista1[3] = {0,1,2};
    double resposta1[3] = {2, 4, 8};
    printf("Teste 1: %f\n\n", Simpson(lista1, resposta1, 3));

    double lista2[6] = {0, 2, 4, 6, 8, 10};
    double resposta2[6] = {3,5,9,10,14,15};
    printf("Teste 2: %f\n\n", Simpson(lista2, resposta2, 6));

    double lista3[8] = {1,2,3,4,5,6,7,8};
    double resposta3[8] = {4,7,9,12,10,14,16,19};
    printf("Teste 3: %f\n\n", Simpson(lista3, resposta3, 8));
    /**/
}