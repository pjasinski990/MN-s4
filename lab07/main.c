#include <stdio.h>

double f(double x);
void diff_quotient(double* xm, double* ym, double* fm, int n);
double poly(double x, int n, double* xm, double* fm);

int main(int argc, char const *argv[])
{
    const int n = 6;
    const double step = 0.01;
    FILE* f_uneven = fopen("zad_1.dat", "w");
    FILE* f_even = fopen("zad_2.dat", "w");

    double xm[] = {-5.0, -2.0, -0.5, 0.0, 0.5, 2.0, 5.0};
    double ym[n+1];
    double fm[n+1];

    double xm_even[n+1];
    double ym_even[n+1];
    double fm_even[n+1];

    for (size_t i = 0; i < n+1; i++)
    {
        ym[i] = f(xm[i]);
        fm[i] = 0.0;

        xm_even[i] = -5.0 + i * (10.0/(n));
        ym_even[i] = f(xm_even[i]);
        fm_even[i] = 0.0;
    }
    
    diff_quotient(xm, ym, fm, n+1);
    diff_quotient(xm_even, ym_even, fm_even, n+1);
    for (double x = -5.0; x <= 5.0; x += step)
    {
        fprintf(f_uneven, "%.2f %f\n", x, poly(x, n, xm, fm));
        fprintf(f_even, "%.2f %f\n", x, poly(x, n, xm_even, fm_even));
    }
    
    printf("Ilorazy roznicowe dla nierowno rozmieszczonych wezlow:\n");
    for (int i = 0; i < n+1; i++) 
    {
    	printf("%g\n", fm[i]);
    }
    printf("\n");
    printf("Ilorazy roznicowe dla rowno rozmieszczonych wezlow:\n");
    for (int i = 0; i < n+1; i++) 
    {
    	printf("%g\n", fm_even[i]);
    }

    fclose(f_uneven);
    fclose(f_even);
    return 0;
}

double f(double x)
{
    return 1.0/(1.0+x*x);
}

void diff_quotient(double* xm, double* ym, double* fm, int n)
{
    for (size_t j = 0; j < n; j++)
    {
        for (size_t k = 0; k <= j; k++)
        {
            double factor = 1;
            for (size_t l = 0; l <= j; l++)
            {
                if (k == l) {continue;}
                factor *= 1.0/(xm[k]-xm[l]);
            }
            fm[j] += ym[k]*factor;
        }
    }
}

double poly(double x, int n, double* xm, double* fm)
{
    double res = 0;
    for (size_t j = 0; j <= n; j++)
    {
        double factor = 1;
        for (size_t i = 0; i < j; i++)
        {
            factor *= x-xm[i];
        }
        res += fm[j] * factor;
    }
    return res;
}
