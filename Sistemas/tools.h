#ifndef TOOLS_H
#define TOOLS_H

void diagonaliza(double* mat, int tam);
int transforma_matriz_simple(double* matriz, int tam, double* MatrizL);
int transforma_matriz(double* matriz, int m, double* decoy);
int resolve_matriz(double* matriz, int tam, double* decoy, double* resultado);
void printa_matriz_matricial(int n, double matriz[][n]);
void printa_matriz(double* matriz, int tam);
void printa_lista(double* lista, int tam);
int pivotiza(double* A, double* b, int n, int coluna, int* perm);
double modulo(double x);
int busca_repetidos(double* Xis, int tam);
double GaussX_CalculaErro(double* atual, double* anterior, int tam);
double* Gauss(double* Xis, double* FXis, int n, int grau);
int GaussJacob_CalculaVetor(double* MatrizCoeficientes, int numLinhas, double* vetAtual, double* vetImagem, int tam);
int GaussJacob(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia, int max_iteracoes, double* respostas);
void GaussSeidel_CalculaVetor(double* MatrizCoeficientes, int numLinhas, double* vetAtual, double* vetImagem, int tam);
double* GaussSeidel(double* FXis, int tam, double* MatrizCoeficientes, int numLinhas, double tolerancia, int max_iteracoes);
int InterpolPol(double* Xis, double* FXis, int n, double* respostas);
int AjustePolinomial(double* Xis, double* Y, int n, int m, double* resposta);
int DecompLU(double* MatrizA, double* VetX, double* VetB, double* MatrizL, double* MatrizU, int tamLinha);
void Newton(double* Xis, double* FXis, int n);
double Simpson(double* Xis, double* FXis, int n);
double Trapezios(double* Xis, double* FXis, int n);

#endif