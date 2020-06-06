#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include "nrutil.c"
#include "gauleg.c"
#include "gaulag.c"
#include "gammln.c"
#include "gauher.c"

float f1(float x);
float c1a(float x1, float x2);

float f2(float x, int k);
float c2a(int k);

float f3(float x, float y);
#define C3_EXACT 0.1919832644

int main(int argc, char const *argv[])
{
    FILE* f_c1 = fopen("c1.dat", "w");
    FILE* f_c2 = fopen("c2.dat", "w");
    FILE* f_c3 = fopen("c3.dat", "w");

    for (int n = 2; n <= 20; n++)
    {
        const float x1 = 0.0f;
        const float x2 = 2.0f;
        float* x = vector(1, n);
        float* w = vector(1, n);

        gauleg(x1, x2, x, w, n);
        
        float sum = 0.0f;
        for (int i = 1; i <= n; i++)
        {
            sum += f1(x[i]) * w[i];
        }
        fprintf(f_c1, "%d %g\n", n, fabs(sum - c1a(0.0f, 2.0f)));
        
        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }


    for (int k = 5; k <= 10; k+=5)
    {
        for (int n = 2; n <= 20; n++)
        {
            float* x = vector(1, n);
            float* w = vector(1, n);
            gaulag(x, w, n, 0);
            
            float sum = 0.0f;
            for (int i = 1; i <= n; i++)
            {
                sum += f2(x[i], k) * w[i];
            }
            fprintf(f_c2, "%d %g\n", n, fabs(sum - c2a(k)));
        
            free_vector(x, 1, n);
            free_vector(w, 1, n);
        }
        fprintf(f_c2, "\n\n");
    }

    
    for (int n = 2; n <= 15; n++)
    {
        float* x = vector(1, n);
        float* w = vector(1, n);

        gauher(x, w, n);
        
        float sum = 0.0f;
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                sum += f3(x[i], x[j]) * w[i] * w[j];
            }
        }
        fprintf(f_c3, "%d %g\n", n, fabs(sum - C3_EXACT));
        
        free_vector(x, 1, n);
        free_vector(w, 1, n);
    }


    fclose(f_c1);
    fclose(f_c2);
    fclose(f_c3);
    return 0;
}

float f1(float x)
{
    return x/(4*x*x + 1);
}
float c1a(float x1, float x2)
{
    return 1.0f/8.0f * (log(fabs(4.0f * x2*x2 + 1.0f)) - log(fabs(4.0f * x1*x1 + 1.0f)));
}

float f2(float x, int k)
{
    return pow(x, k);
}
float c2a(int k)
{
    int res = 1;
    for (int i = 1; i < k; i++, res*=i);
    return res;
}

float f3(float x, float y)
{
    return pow(sin(x), 2) * pow(sin(y), 4);
}