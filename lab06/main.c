#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IT_MAX 30
double licz_r(double* a, double* b, int n, double xj);

int main(int argc, char const *argv[])
{
    const int N = 5;
    double *a, *b;
    double x0, x1, x2, R0, R1, R2;

    a = calloc(N+1, sizeof(double));
    b = calloc(N+1, sizeof(double));
    a[0] = 240.0;
    a[1] = -196.0;
    a[2] = -92.0;
    a[3] = 33.0;
    a[4] = 14.0;
    a[5] = 1.0;

    printf("L\tit\txit+1\tRit+1\n");
    for (int L = 1; L <= N; L++)
    {
        int n = N - L + 1;
        x0 = 0.0;
        x1 = 0.1;
        R0 = licz_r(a, b, n, x0);
        R1 = licz_r(a, b, n, x1);

        for (int i = 1; i < IT_MAX; i++)
        {
            x2 = x1 - R1*(x1-x0)/(R1-R0);
            R2 = licz_r(a, b, n, x2);
            R0 = R1;
            R1 = R2;

            x0 = x1;
            x1 = x2;
            printf("%d\t%d\t%.4g\t%g\n", L, i, x2, R2);
            if (fabs(x1-x0) < 10e-8) { break; }
        }

        printf("\n");
        for (int i = 0; i < n; i++)
        {
            a[i] = b[i];
        }
    }

    free(a);
    free(b);
    return 0;
}

double licz_r(double* a, double* b, int n, double xj)
{
    b[n] = 0;

    for (int i = n-1; i >= 0; i--)
    {
        b[i] = a[i+1] + xj*b[i+1];
    }
    return a[0] + xj*b[0];
}