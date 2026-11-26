#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define VARIACAO_TOLERADA 0.0000001

// Compilar com -ffloat-store!
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
        if(modulo(Xis[i] - Xis[i+1]) < VARIACAO_TOLERADA)
        {
            return -1;
        }
    }
    return 0;
}

/***********************************************************************************
**************************** Método Eliminação de Gauss ****************************
************************************************************************************/

// Recebe os valores a serem somados em num, a quantidade de números em tam e a potencia à qual serão elevados em potencia 
void troca_linhas(double* matriz, double* Y, int n, int dest, int origin)
{
    for(int i = 0 ; i < n ; i++)
    {
        double temp = matriz[dest*n + i];
        matriz[dest*n + i] = matriz[origin*n + i];
        matriz[origin*n + i] = temp;
    }
    double temp = Y[dest];
    Y[dest] = Y[origin];
    Y[origin] = temp;
}

void pivotiza(double* matriz, int n, double* Y, int coluna)
{
        int linha_maior = coluna;
        for(int j = coluna ; j < n ; j++)
        {
            if(modulo(matriz[n*j + coluna]) > modulo(matriz[n*linha_maior + coluna]))
            {
                linha_maior = j;
            }
        }
        troca_linhas(matriz, Y, n, coluna, linha_maior);
}

double Pow(double base, int expoente)
{
    double resultado = 1;
    for(int i = 0 ; i < expoente ; i++)
    {
        resultado *= base;
    }
    return resultado;
}

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

int transforma_matriz_simple(double* matriz, int tam, double* MatrizL)
{
    double multiplicador = 0;
    for( int linha_ref = 0 ; linha_ref < tam ; linha_ref++)
    {
        if(*(matriz + linha_ref * tam + linha_ref) < VARIACAO_TOLERADA)
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
}

int transforma_matriz(double* matriz, int m, double* decoy)
{
    double multiplicador = 0;
    for( int linha_ref = 0 ; linha_ref < m - 1; linha_ref++)
    {
        pivotiza(matriz, m, decoy, linha_ref);
        if(*(matriz + linha_ref*m + linha_ref) < VARIACAO_TOLERADA)
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
    
}

int resolve_matriz(double* matriz, int tam, double* decoy, double* resultado)
{
    // Incluir checagem da matriz ser resolvível ou não *insira emoji de caveira agui
    for(int i = tam-1 ; i >= 0 ; i--)
    {
        int j = tam-1;
        double temp = decoy[i];

        for(j; j > i ; j--)
        {
            temp -= *(matriz + i * tam + j) * resultado[j];
        }
        double flag = *(matriz + i * tam + j);
        
        if(flag < VARIACAO_TOLERADA)
        {
            return -1;
        }

        resultado[j] = temp / *(matriz + i * tam + j);
    }
    return 0;
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

    if(modulo_atual < VARIACAO_TOLERADA)
    {
        return modulo_temp;
    }
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
            temp -= vetTemp[j] * (*(MatrizCoeficientes + i * tam + j));
        }
        temp /= (*(MatrizCoeficientes + i * tam + i));
        //printf("Temp[%d]: %f\n", i, temp);
        vetAtual[i] = temp;
    }
}

double* GaussJacob(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia, int max_iteracoes)
{
    double* respostas = calloc(tam, sizeof(double));
    double anterior[tam];
    memset(anterior, 0, sizeof(double)*tam);
    int i = 0;
    do
    { 
        memcpy(anterior, respostas, sizeof(double)*tam);
        GaussJacob_CalculaVetor(MatrizCoeficientes, numLinhas, respostas, anterior, FXis, tam);
        i++;
        //printa_lista(anterior, tam);
        //printf("\n");
        //printa_lista(respostas, tam);
        //printf("\n");
        //printf("deu em: %d\n", GaussJacob_CalculaErro(respostas, anterior, tam) <= tolerancia);
    } while(GaussX_CalculaErro(respostas, anterior, tam) > tolerancia && max_iteracoes > i);
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

double* GaussSeidel(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia, int max_iteracoes)
{
    double* respostas = calloc(tam, sizeof(double));
    double anterior[tam];
    memset(anterior, 0, sizeof(double)*tam);
    int i = 0;
    do
    {
        memcpy(anterior, respostas, sizeof(double)*tam);
        GaussSeidel_CalculaVetor(MatrizCoeficientes, numLinhas, respostas, anterior, FXis, tam);
        i++;
        //printa_lista(anterior, tam);
        //printf("\n");
        //printa_lista(respostas, tam);
        //printf("\n");
        //printf("deu em: %d\n", GaussJacob_CalculaErro(respostas, anterior, tam) <= tolerancia);
    } while(GaussX_CalculaErro(respostas, anterior, tam) > tolerancia && max_iteracoes > i);
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
    if(n < 2)
    {

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
    /**
    printa_lista(S, 2*m+1);
    printf("\n\n");
    printa_lista(T, m+1);
    /**/
}

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
    double* CopiaFXis = malloc(sizeof(double)*n);
    memcpy(CopiaFXis, FXis, sizeof(double)*n);
    double* MatX = malloc(sizeof(double)*n*n);

    int flag = busca_repetidos(Xis, n);

    if(flag == -1)
    {
        printf("Valores de X repetidos!\nAbortando operacao\n");
        return NULL;
    }

    InterpolPol_GeraMatX(MatX, Xis, n);
    //printa_matriz(MatX, n);
    //printf("\n\n");

    flag = transforma_matriz(MatX, n, CopiaFXis);
    if(flag == -1)
    {
        printf("Operacao Abortada!\nMatriz sob transformacao é singular. Insira outra matriz.");
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
        return -1;
    }
    free(MatX);
    free(CopiaFXis);
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
        if(modulo((Xis[i+1] - Xis[i]) - step) > VARIACAO_TOLERADA || modulo((Xis[i+2] - Xis[i+1]) - step) > VARIACAO_TOLERADA)
        {
            printf("Pontos inseridos não estao igualmente espaçados!\nOperação cancelada");
            return 0.0;
        }
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
    double GX[4] = {1,2,3,4};
    double GY[4] = {6,11,18,27};
    resposta = Gauss(GX, GY, 4, 2);
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

    /* Ajuste Polinomial */
    /**/
    double PX[3] = {1,2,3};
    double PY[3] = {2,4,5};
    double Resp[3] = {};
    AjustePolinomial(PX, PY, 3, 1, Resp);

    double PY2[3] = {4, 4, 2};
    double Resp2[3] = {};
    AjustePolinomial(PX, PY2, 3, 2, Resp2);
    /**/

    /* Interpol */
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
    /**
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