#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>

double f1(double x);
double f1prim(double x);
double f2(double x);
double f3(double x);
double F(double x, double* a, double* b, int MS, int MC);

int main(int argc, char const *argv[])
{
    srand(time(0));

    FILE* f_points = fopen("pkt.dat", "w");
    FILE* f_data1 = fopen("approx1.dat", "w");
    FILE* f_data2 = fopen("approx2.dat", "w");
    FILE* f_data3 = fopen("approx3.dat", "w");
    FILE* f_data4 = fopen("approx4.dat", "w");
    FILE* f_ab = fopen("ab4.dat", "w");
    FILE* datafiles[] = {f_data1, f_data2, f_data3, f_data4};

    const int n = 100;
    const int MS[] = {5, 10, 30};
    const int MC[] = {5, 10, 30};
    const double xmin = 0.0;
    const double xmax = 2*M_PI;
    const double h = (xmax-xmin) / (n-1);
    const double dx = 0.01;
    
    double* a = calloc(sizeof(double), MS[0] + 2);
    double* b = calloc(sizeof(double), MC[0] + 2);
    double* aprim = calloc(sizeof(double), MS[1] + 2);
    double* bprim = calloc(sizeof(double), MC[1] + 2);
    double* abis = calloc(sizeof(double), MS[2] + 2);
    double* bbis = calloc(sizeof(double), MC[2] + 2);
    double* aparams[] = {a, aprim, abis};
    double* bparams[] = {b, bprim, bbis};

    double* x = calloc(sizeof(double), n);
    double* y1 = calloc(sizeof(double), n);
    double* y2 = calloc(sizeof(double), n);
    double* y3 = calloc(sizeof(double), n);
    double* y4 = calloc(sizeof(double), n);
    double (*functions[])(double) = {f1, f2, f3, f1prim};
    
    for (int i = 0; i < n; i++)
    {
        x[i] = i * h;
        y1[i] = f1(x[i]);
        y2[i] = f2(x[i]);
        y3[i] = f3(x[i]);
        y4[i] = f1prim(x[i]);
    }

    // Printing to points file
    for (int i = 0; i < sizeof(functions)/sizeof(*functions)-1; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf(f_points, "%f %f\n", x[j], functions[i](x[j]));
        }
        fprintf(f_points, "\n\n");
    }
    for (int j = 0; j < n; j++)
    {
        fprintf(f_points, "%f %f\n", x[j], y4[j]);
    }

    // Printing to approx1
    for (int i = 0; i <= MS[0] + 1; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += y1[j] * sin(i*x[j]);
        }
        a[i] = 2.0*sum / n;
    }
    for (int i = 0; i <= MC[0] + 1; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += y1[j] * cos(i*x[j]);
        }
        b[i] = 2.0*sum / n;
    }
    for (double d = xmin; d <= xmax; d+=dx)
    {
        fprintf(f_data1, "%f %f\n", d, F(d, a, b, MS[0], MC[0]));
    }

    // Printing to approx2
    for (int i = 0; i <= MS[0] + 1; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += y2[j] * sin(i*x[j]);
        }
        a[i] = 2.0*sum / n;
    }
    for (int i = 0; i <= MC[0] + 1; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += y2[j] * cos(i*x[j]);
        }
        b[i] = 2.0*sum / n;
    }
    for (double d = xmin; d <= xmax; d+=dx)
    {
        fprintf(f_data2, "%f %f\n", d, F(d, a, b, MS[0], MC[0]));
    }
    
    // Printing to approx3
    // MS = 5, MC = 0
    for (int i = 0; i <= MS[0] + 1; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += y3[j] * sin(i*x[j]);
        }
        a[i] = 2.0*sum / n;
    }
    double sum = 0.0;
    for (int j = 0; j < n; j++)
    {
        sum += y3[j]; // * cos(0);
    }
    b[0] = 2.0*sum / n;

    for (double d = xmin; d <= xmax; d+=dx)
    {
        fprintf(f_data3, "%f %f\n", d, F(d, a, b, MS[0], 0));
    }
    fprintf(f_data3, "\n\n");

    // Rest of MS, MC values   
    for (int m = 0; m < 2; m++)
    {
        for (int i = 0; i <= MS[m] + 1; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                sum += y3[j] * sin(i*x[j]);
            }
            aparams[m][i] = 2.0*sum / n;
        }
        for (int i = 0; i <= MC[m] + 1; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                sum += y3[j] * cos(i*x[j]);
            }
            bparams[m][i] = 2.0*sum / n;
        }
        for (double d = xmin; d <= xmax; d+=dx)
        {
            fprintf(f_data3, "%f %f\n", d, F(d, aparams[m], bparams[m], MS[m], MC[m]));
        }
        fprintf(f_data3, "\n\n");
    }

    // Printing to approx4
    for (int m = 0; m < 3; m+=2)
    {
        for (int i = 0; i <= MS[m] + 1; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                sum += y4[j] * sin(i*x[j]);
            }
            aparams[m][i] = 2.0*sum / n;
        }
        for (int i = 0; i <= MC[m] + 1; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                sum += y4[j] * cos(i*x[j]);
            }
            bparams[m][i] = 2.0*sum / n;
        }
        for (double d = xmin; d <= xmax; d+=dx)
        {
            fprintf(f_data4, "%f %f\n", d, F(d, aparams[m], bparams[m], MS[m], MC[m]));
        }
        fprintf(f_data4, "\n\n");
    }

    // Printing to ab file
    for (int m = 0; m < 3; m+=2)
    {
        for (int i = 0; i <= MS[m]; i++)
        {
            fprintf(f_ab, "%d, %f, %f\n", i, aparams[m][i], bparams[m][i]);
        }
        fprintf(f_ab, "\n\n");
    }
    
    free(x);
    free(y1);
    free(y2);
    free(y3);
    free(y4);
    free(a);
    free(aprim);
    free(abis);
    free(b);
    free(bprim);
    free(bbis);
    fclose(f_points);
    fclose(f_data1);
    fclose(f_data2);
    fclose(f_data3);
    fclose(f_data4);
    fclose(f_ab);
    return 0;
}


double f1(double x)
{
    return 2.0*sin(x) + sin(2.0*x) + 2.0*sin(3.0*x);
}

double f2(double x)
{
    return 2.0*sin(x) + sin(2.0*x) + 2.0*cos(x) + cos(2.0*x);
}

double f3(double x)
{
    return 2.0*sin(1.1*x) + sin(2.1*x) + 2.0*sin(3.1*x);
}

double f1prim(double x)
{
    double alpha = rand() / (RAND_MAX + 1.0) - 0.5;
    return 2.0*sin(x) + sin(2.0*x) + 2.0*sin(3.0*x) + alpha;
}

double F(double x, double* a, double* b, int MS, int MC)
{
    double sum = 0.0;
    for (int i = 1; i <= MS; i++)
    {
        sum += a[i]*sin(i*x);
    }
    for (int i = 1; i <= MC; i++)
    {
        sum += b[i]*cos(i*x);
    }
    return sum + b[0]/2.0;
}