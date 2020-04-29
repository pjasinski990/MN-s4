#include <stdio.h>
#include <stdlib.h>

#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c"

float f1(float x);
float f2(float x);
float s(float h, int n, float* xw, float* c, float x);
float phi(float* xw, float x, int i, float h);

int main(int argc, char const *argv[])
{
    const float xmin = -5.0f;
    const float xmax = 5.0f;
    const float dx = 0.01f;
    FILE* fout1 = fopen("f1.dat", "w");
    FILE* fout2 = fopen("f2.dat", "w");

    const int sizes_for_f1[] = {5, 6, 10, 20};
    const int sizes_for_f2[] = {6, 7, 14};

    /*
     * Looping with function f1
     */
    for (int i = 0; i < sizeof(sizes_for_f1)/sizeof(*sizes_for_f1); i++)
    {
        int n = sizes_for_f1[i];
        float h = (xmax-xmin) / (n-1.0f);
        float* xw = vector(-2, n+3);
        float* yw = vector(1, n);
        float** mA = matrix(1, n, 1, n);
        float** b = matrix(1, n, 1, 1);
        float alpha = (f1(xmin+dx) - f1(xmin-dx)) / (2.0f*dx);
        float beta = (f1(xmax+dx) - f1(xmax-dx)) / (2.0f*dx);
        float* c = vector(0, n+1);
    
        // Filling stuff
        for (int i = -2; i <= n+3; i++) 
        {
            xw[i] = xmin + (i-1)*h;
        }
        for (size_t i = 1; i <= n; i++)
        {
            for (size_t j = 1; j <= n; j++)
            {
                mA[i][j] = 0.0f;
            }
            mA[i][i] = 4.0f;
            yw[i] = b[i][1] = f1(xw[i]);
        }
        b[1][1] += h/3.0f * alpha;
        b[n][1] -= h/3.0f * beta;

        for (size_t i = 1; i < n; i++)
        {
            mA[i+1][i] = 1.0f;
            mA[i][i+1] = 1.0f;
        }
        mA[1][2] = mA[n][n-1] = 2.0f;

        // Calculating c vector
        gaussj(mA, n, b, 1);
        for (size_t i = 1; i <= n; i++)
        {
            c[i] = b[i][1];
        }
        c[0] = c[2] - h/3.0f * alpha;
        c[n+1] = c[n-1] + h/3.0f * beta;
        
        for (double d = xmin; d <= xmax; d+=dx)
        {
            fprintf(fout1, "%g %g\n", d, s(h, n, xw, c, d));
        }
        fprintf(fout1, "\n\n");

        free_vector(xw, -2, n+3);
        free_vector(yw, 1, n);
        free_matrix(mA, 1, n, 1, n);
        free_matrix(b, 1, n, 1, 1);
        free_vector(c, 0, n+1);
    }
    
    /*
     * Looping with function f2
     */
    for (int i = 0; i < sizeof(sizes_for_f2)/sizeof(*sizes_for_f2); i++)
    {
        int n = sizes_for_f2[i];
        float h = (xmax-xmin) / (n-1.0f);
        float* xw = vector(-2, n+3);
        float* yw = vector(1, n);
        float** mA = matrix(1, n, 1, n);
        float** b = matrix(1, n, 1, 1);
        float alpha = (f2(xmin+dx) - f2(xmin-dx)) / (2.0f*dx);
        float beta = (f2(xmax+dx) - f2(xmax-dx)) / (2.0f*dx);
        float* c = vector(0, n+1);
    
        // Filling stuff
        for (int i = -2; i <= n+3; i++) 
        {
            xw[i] = xmin + (i-1)*h;
        }
        for (size_t i = 1; i <= n; i++)
        {
            for (size_t j = 1; j <= n; j++)
            {
                mA[i][j] = 0.0f;
            }
            mA[i][i] = 4.0f;
            yw[i] = b[i][1] = f2(xw[i]);
        }
        b[1][1] += h/3.0f * alpha;
        b[n][1] -= h/3.0f * beta;

        for (size_t i = 1; i < n; i++)
        {
            mA[i+1][i] = 1.0f;
            mA[i][i+1] = 1.0f;
        }
        mA[1][2] = mA[n][n-1] = 2.0f;

        // Calculating c vector
        gaussj(mA, n, b, 1);
        for (size_t i = 1; i <= n; i++)
        {
            c[i] = b[i][1];
        }
        c[0] = c[2] - h/3.0f * alpha;
        c[n+1] = c[n-1] + h/3.0f * beta;
        

        for (double i = xmin; i <= xmax; i+=dx)
        {
            fprintf(fout2, "%g %g\n", i, s(h, n, xw, c, i));
        }
        fprintf(fout2, "\n\n");

        free_vector(xw, -2, n+3);
        free_vector(yw, 1, n);
        free_matrix(mA, 1, n, 1, n);
        free_matrix(b, 1, n, 1, 1);
        free_vector(c, 0, n+1);
    }

    fclose(fout1);
    fclose(fout2);

    return 0;
}

float f1(float x) 
{
    return 1.0f/(1.0f+x*x);
}

float f2(float x)
{
    return cos(2*x);
}

float s(float h, int n, float* xw, float* c, float x)
{
    float res = 0.0f;
    for (size_t i = 0; i <= n+1; i++)
    {
        res += c[i] * phi(xw, x, i, h);
    }
    return res;
}

float phi(float* xw, float x, int i, float h)
{
    if (x >= xw[i-2] && x < xw[i-1])
    {
        return pow(x-xw[i-2], 3) / pow(h, 3);
    }
    else if (x >= xw[i-1] && x < xw[i])
    {
        return (pow(h, 3) + 3*h*h*(x-xw[i-1]) + 3*h*(x-xw[i-1])*(x-xw[i-1]) - 3*pow(x-xw[i-1], 3)) / pow(h, 3);
    }
    else if (x >= xw[i] && x < xw[i+1])
    {
        return (pow(h, 3) + 3*h*h*(xw[i+1]-x) + 3*h*(xw[i+1]-x)*(xw[i+1]-x) - 3*pow(xw[i+1]-x, 3)) / pow(h, 3);
    }
    else if (x >= xw[i+1] && x <= xw[i+2])
    {
        return pow(xw[i+2]-x, 3) / pow(h, 3);
    }
    else return 0.0f;
}
