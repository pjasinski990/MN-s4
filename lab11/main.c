#define _USE_MATH_DEFINES
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"
#include "four1.c"

float fun(float x, float period);
float g(float x, float sigma);
float randd();

int main()
{
    FILE* k8 = fopen("k8.dat", "w");
    FILE* k10 = fopen("k10.dat", "w");
    FILE* k12 = fopen("k12.dat", "w");
    FILE* filearray[] = {k8, k10, k12};
    int filecounter = 0;

    for (int k = 8; k <= 12; k += 2)
    {
        const int n = pow(2, k);
        float *f, *g1, *g2;

        f = vector(1, 2*n);
        g1 = vector(1, 2*n);
        g2 = vector(1, 2*n);

        const float T = 1.0f;
        const float tmax = 3.0f * T;
        float dt = tmax / n;
        float sigma = T / 20.0f;

        for (int i = 1; i <= n; i++)
        {
            float ti = dt * (i - 1.0f);
            f[2 * i - 1] = fun(ti, 1.0f) + randd();
            g1[2 * i - 1] = g2[2 * i - 1] = g(ti, sigma);
            f[2 * i] = g1[2 * i] = g2[2 * i] = 0.0f;
            fprintf(filearray[filecounter], "%g %g\n", ti, f[2*i - 1]);
        }
        fprintf(filearray[filecounter], "\n\n");

        four1(f, n, 1);
        four1(g1, n, 1);
        four1(g2, n, -1);

        for (int i = 1; i <= n; i++)
        {
            float a1 = f[2*i - 1];
            float b1 = f[2*i];
            float a2 = g1[2*i - 1] + g2[2*i - 1];
            float b2 = g1[2*i] + g2[2*i];
            f[2*i - 1] = a1*a2 - b1*b2;
            f[2*i] = a1*b2 + a2*b1;
        }

        four1(f, n, -1);

        float max_val = fabs(f[1]);
        for (int i = 2; i <= n; i++)
        {
            if (max_val < fabs(f[2*i - 1]))
            {
                max_val = fabs(f[2*i - 1]);
            }
        }

        for (int i = 1; i <= n; i++)
        {
            float ti = dt * (i - 1);
            fprintf(filearray[filecounter], "%g %g\n", ti, f[2*i - 1] * 2.5f / max_val);
        }

        free_vector(f, 1, 2 * n);
        free_vector(g1, 1, 2 * n);
        free_vector(g2, 1, 2 * n);
        ++filecounter;
    }

    fclose(k8);
    fclose(k10);
    fclose(k12);

    return 0;
}

float fun(float x, float period)
{
    return sin((2 * M_PI / period) * x) + sin(2 * (2 * M_PI / period) * x) + sin(3 * (2 * M_PI / period) * x);
}

float g(float x, float sigma)
{
    return (1.0f / (sigma * sqrt(2*M_PI))) * exp(-1.0f*(x*x) / (2.0f * sigma*sigma));
}

float randd()
{
    return (float)rand() / RAND_MAX - 0.5f;
}