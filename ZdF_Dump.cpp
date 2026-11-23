#include <iostream>
#include <functional>
#include <cmath>
#include <iomanip>
using namespace std;



int teste(double a)
{
    if (a) return 1;
    else return 0;
}

double ajustaDominioLog(double *calculoDerivada, double a1, int tentativas)
{
    int i = 0;
    double a2 = a1 - (*calculoDerivada);
    while(a2 <= 0 && i < tentativas)
    {
        (*calculoDerivada) /= 2.0;
        a2 = a1 - (*calculoDerivada);
        i++;
    }
    return a2;

}

double derivadaNumerica(std::function<double(double)> f, double x, double h = 1e-5) {
    /**
    cout << std::setprecision(15);
    cout << "Eu aqui: " << h << "\n\n";
    cout << "f(xs): " << (f(x+h)) << " e " << f(x) << "\n";
    long double i = f(x+h);
    long double j = f(x);
    cout << "i - j" << i - j << "\n";
    cout << "Diferença: " << f(x+h) - f(x) << "\n\n";
    /**/
    return (f(x + h) - f(x)) / h;
}

double lineari(double a)
{
    return (20*a + 30);
}

double quinta(double x)
{
    return (5*pow(x,5) + 8*pow(x,4) + 3*pow(x,3) + pow(x,2) + 7*x + 20);
}

double lograritmo(double x)
{
    /**
    cout << "x ---->" << x << "\n\n";
    cout << " 5 - 3 * log10(x) ---->" << 5 - 3 * log10(x) << "\n\n";
    /**/
    return (5 - 3*log10(x));
}

double ComSeno(double x)
{
    return(40*cos(x) - 8);
}


double modulo(double a)
{
    if( a < 0 )
    {
        return a*-1;
    }
    else
    {
        return a;
    }

}

void bissecao(std::function<double(double)> f, double a, double b, double tsi)
{
    double c, fc;

    cout << "Função: 20x + 30" << '\n';

    if(f(a) * f(b) >= 0)
    {
        cout << "Pontos inválidos. Escolha novamente.";
        return;
    }
    // int i = 0;
    do{
        c = (a+b)/2;
        fc = f(c);

        // cout << c << '\n';
        // cout << fc << '\n';

        if(fc * f(a) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }while(modulo(fc) > tsi && modulo(a-b) > tsi);

    cout << "f(x): " << fc << '\n';
    cout << "x:" << c << '\n';

    /**/
    if(fc <= tsi)
    {
        cout << "f(x) < tsi\n";
    }

    if(c <= tsi)
    {
        cout << "b - a < tsi\n";
    }
    /**/
}

void falsa_posicao(std::function<double(double)> f, double a, double b, double tsi)
{
    double c, fc;

    cout << "Função: 5*x^5 + 8*x^4 + 3*x^3 + x^2 + 7x + 20" << '\n';

    if(f(a) * f(b) >= 0)
    {
        cout << "Pontos inválidos. Escolha novamente.";
        return;
    }

    do{
        c = b - (f(b) * (b-a))/(f(b) - f(a));
        // cout << c << '\n';
        // cout << fc << '\n';
        fc = f(c);

        if(fc * f(a) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }while(modulo(fc) > tsi && modulo(a-b) > tsi);

    cout << "f(x): " << fc << '\n';
    cout << "x:" << c << '\n';
    cout << "confirmando: " << f(c) << '\n';

    /**/
    if(fc <= tsi)
    {
        cout << "f(x) < tsi\n";
    }

    if(c <= tsi)
    {
        cout << "b - a < tsi\n";
    }
    /**/

}

void newton(std::function<double(double)> f, double a1, double tsi, int Nmax)
{
    double temp = (f(a1)/derivadaNumerica(f, a1));
    double a2 = ajustaDominioLog(&temp, a1, 50), fa2 = f(a2);
    if (a2 < 0)
    {
        cout << "Falha ao ajustar domínio. Escolha outro número!!\n\n";
        return;
    }
    cout << "Função: 5 - 3*log10(x)" << '\n';
    int i;
    for( i = 0; modulo(fa2) >= tsi && modulo(a2) != 0 && modulo(a2-a1)/modulo(a2) >= tsi && i < Nmax; i++)
    {
        a1 = a2;
        temp = (f(a1)/derivadaNumerica(f, a1));
        a2 = ajustaDominioLog(&temp, a1, 50);
        if (a2 < 0)
        {
            cout << "Falha ao ajustar domínio. Escolha outro número!!\n\n";
            return;
        }
        fa2 = f(a2);
    }

    // cout << "i: " << i << '\n';
    cout << "f(x): " << fa2 << '\n';
    cout << "x:" << a2 << '\n';

    /**/
    if(modulo(fa2) < tsi)
    {
        cout << "f(x) < tsi\n";
    }

    if(modulo(a2-a1)/modulo(a2) < tsi)
    {
        cout << "variacao < tsi\n";
    }
    /**/
}

void secante(std::function<double(double)> f, double x1, double x0, double tsi, int Nmax)
{
    double x2 = x1 - (f(x1) * (x1-x0))/(f(x1) - f(x0));
    
    cout << "Função: 50*cos(x) - 8" << '\n';

    for( int i = 0; modulo(f(x2)) >= tsi && modulo(x1-x0) >= tsi && i < Nmax; i++)
    {
        x0 = x1;
        x1 = x2;
        if(f(x1) - f(x0) == 0)
        {
            break;
        }
        x2 = x1 - (f(x1) * (x1-x0))/(f(x1) - f(x0));
    }

    cout << "f(x): " << f(x2) << '\n';
    cout << "x:" << x2 << '\n';

    /**/
    if(modulo(f(x2)) < tsi)
    {
        cout << "f(x) < tsi\n";
    }

    if(modulo(x1-x0) < tsi)
    {
        cout << "variacão< tsi\n";
    }
    /**/
}

void chamarTodasPreDefinido()
{
    cout << "\n-------------------";
    cout << "\nBissecao Linear\n";
    bissecao(lineari, -5, 5, 0.001);

    cout << "\n-------------------";
    cout << "\nFalsa Posicao\n";
    falsa_posicao(quinta, -5, 5, 0.001);

    cout << "\n-------------------";
    cout << "\nNewton\n";
    newton(lograritmo, 200, 0.001, 100);

    cout << "\n-------------------";
    cout << "\nSecante\n";
    secante(ComSeno, 4, 1, 0.001, 100);
}

void chamarTodasCustom(double b1, double b2, double f1, double f2, double n1, double s1, double s2)
{
    cout << "\n-------------------";
    cout << "\nBissecao Linear\n";
    bissecao(lineari, b1, b2, 0.001);

    cout << "\n-------------------";
    cout << "\nFalsa Posicao\n";
    falsa_posicao(quinta, f1, f2, 0.001);

    cout << "\n-------------------";
    cout << "\nNewton\n";
    newton(lograritmo, n1, 0.001, 100);

    cout << "\n-------------------";
    cout << "\nSecante\n";
    secante(ComSeno, s1, s2, 0.001, 100);    
}

int main()
{
    double b1,b2,f1,f2,n1,s1,s2;
    int flag = 0;
    int comando = 0;
    do{
        cout << "Insira o comando a ser executado\n";
        cout << "0 - Encerrar\n1 - Executar com valores pre-definidos\n2 - Executar com valores customizados\n3 - Executar novamente\nEscolha: ";
        cin >> comando;

        if(comando == 1)
        {
            chamarTodasPreDefinido();
        }
        else if (comando == 2)
        {
            flag = 1;
            // Bisseção
            cout << "Insira o primeiro valor para a bissecao: ";
            cin >> b1;
            cout << "Insira o segundo valor para a bissecao: ";
            cin >> b2;

            // falsa posicao
            cout << "Insira o primeiro valor para a falsa posicao: ";
            cin >> f1;
            cout << "Insira o segundo valor para a falsa posicao: ";
            cin >> f2;

            // newton
            cout << "Insira o valor para o met. de newton: ";
            cin >> n1;

            // secante
            cout << "Insira o primeiro valor para a secante: ";
            cin >> s1;
            cout << "Insira o segundo valor para a secante: ";
            cin >> s2;

            chamarTodasCustom(b1,b2,f1,f2,n1,s1,s2);
        }
        else if (comando == 3)
        {
            if (flag == 0)
            {
                cout << "Sem valores customizados. Utilizando valores pre-definidos\n";
                chamarTodasPreDefinido();
            }
        }
        else if (comando != 0) 
        {
            cout << "\n\nBah. Escolhe outro\n\n";
        }
        cout << "\n-------------------\n";
    }while(comando != 0);
    cout << "Programa Encerrado\n";
    
}