#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"

double f1(double x);
double f2(double x);
double f3(double x);

void rombergCalculate(double** matrix, int m_size, double a, double b, double(*f)(double));

int main(int argc, char const *argv[])
{
    double** mD1 = dmatrix(0, 7, 0, 7);
    double** mD2 = dmatrix(0, 15, 0, 15);
    double** mD3 = dmatrix(0, 7, 0, 7);

    FILE* f_fun1 = fopen("fun1.dat", "w");
    FILE* f_fun2 = fopen("fun2.dat", "w");
    FILE* f_fun3 = fopen("fun3.dat", "w");

    rombergCalculate(mD1, 8, 0, 1, f1);
    rombergCalculate(mD2, 16, -1, 1, f2);
    rombergCalculate(mD3, 8, 0, 1, f3);


    const char* header = "D[w][0]\t\tD[w][w]\n";
    fprintf(f_fun1, header);
    fprintf(f_fun2, header);
    fprintf(f_fun3, header);

    for (int i = 0; i < 8; i++)
    {
        fprintf(f_fun1, "%.10f\t%.10f\n", mD1[i][0], mD1[i][i]);
        fprintf(f_fun3, "%.10f\t%.10f\n", mD3[i][0], mD3[i][i]);
    }

    for (int i = 0; i < 16; i++)
    {
        fprintf(f_fun2, "%.10f\t%.10f\n", mD2[i][0], mD2[i][i]);
    }
    

    fclose(f_fun1);
    fclose(f_fun2);
    fclose(f_fun3);

    free_dmatrix(mD1, 0, 7, 0, 7);
    free_dmatrix(mD2, 0, 15, 0, 15);
    free_dmatrix(mD3, 0, 7, 0, 7);
    return 0;
}

double f1(double x)
{
    if (fabs(x) < 1e-8) {x+= 1e-7;}
    return sin(x)/x;
}

double f2(double x)
{
    if (fabs(x) < 1e-8) {x+= 1e-7;}
    return (cos(x) - exp(x)) / sin(x);
}

double f3(double x)
{
    if (fabs(x) < 1e-8) {x+= 1e-7;}
    return exp(-1/x)/x;
}

void rombergCalculate(double** matrix, int m_size, double a, double b, double(*f)(double))
{
    matrix[0][0] = 0.5 * (b-a) * (f(a)+f(b));
    for (int w = 1; w < m_size; w++)
    {
        double hw = (b - a) / pow(2, w);

        double sum = 0.0;
        for (int i = 1; i <= pow(2, w-1); i++)
        {
            sum += f(a + (2*i - 1) * hw);
        }
        matrix[w][0] = 0.5 * matrix[w-1][0] + hw * sum;

        for (int k = 1; k < m_size; k++)
        {
            for (int w = k; w < m_size; w++)
            {
                matrix[w][k] = (pow(4, k) * matrix[w][k-1] - matrix[w-1][k-1]) / (pow(4, k) - 1);
            }
        }
    }
}