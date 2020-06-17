#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gammln.c"
#include "gammp.c"
#include "gcf.c"
#include "gser.c"
#include "nrutil.h"
#include "nrutil.c"

double rand_linear(long a, long c, long m, int version);
double find_mu(double* numbers, int N);
double find_var(double* numbers, int N);
void make_histogram(double* numbers, int N, double min, double max, int* bins, int binN);

double normal_distribution(double x, double mu, double sigma);
double calc_pj(double xmin, double xmax, double mu, double sigma);

double normal_density(double x, double mu, double sigma);
double rand_normal(double mu, double sigma, double d);

int main(int argc, char const *argv[])
{
    const int N = 10000;

    long a1 = 123, c1 = 1, m1 = pow(2, 15);
    long a2 = 69069, c2 = 1, m2 = pow(2, 32);
    double* linear_arr1 = malloc(sizeof(double) * (N+1));
    double* linear_arr2 = malloc(sizeof(double) * (N+1));
    int linear_hist1[12];
    int linear_hist2[12];
    FILE* file_linear = fopen("U.dat", "w");
    FILE* file_linear_hist = fopen("U_hist.dat", "w");

    double last1 = 10.0 / (m1 + 1.0);
    double last2 = 10.0 / (m2 + 1.0);
    linear_arr1[0] = last1;
    linear_arr2[0] = last2;
    double curr1, curr2;

    // Uniform parameters
    for (int i = 1; i <= N; i++)
    {
        curr1 = rand_linear(a1, c1, m1, 1);
        curr2 = rand_linear(a2, c2, m2, 2);
        fprintf(file_linear, "%g %g\n", last1, curr1);
        last1 = curr1;
        last2 = curr2;
        
        linear_arr1[i] = last1;
        linear_arr2[i] = last2;
    }
    fprintf(file_linear, "\n\n");
    for (int i = 1; i <= N; i++)
    {
        fprintf(file_linear, "%g %g\n", linear_arr2[i-1], linear_arr2[i]);
    }

    double uniform_mu1 = find_mu(linear_arr1+1, N);
    double uniform_mu2 = find_mu(linear_arr2+1, N);
    double uniform_sigma1 = sqrt(find_var(linear_arr1+1, N));
    double uniform_sigma2 = sqrt(find_var(linear_arr2+1, N));

    printf("Uniform distribution\n");
    printf("mu values: %g, %g\n", uniform_mu1, uniform_mu2);
    printf("sigma values: %g, %g\n\n", uniform_sigma1, uniform_sigma2);

    // Uniform histogram 
    make_histogram(linear_arr1+1, N, 0.0, 1.0, linear_hist1, 12);
    make_histogram(linear_arr2+1, N, 0.0, 1.0, linear_hist2, 12);
    const double bin_width_uniform = 1.0 / 12.0;
    for (int i = 0; i < 12; i++)
    {
        fprintf(file_linear_hist, "%g %g\n", bin_width_uniform/2.0 + i*bin_width_uniform, (double)linear_hist1[i]/N);
    }
    fprintf(file_linear_hist, "\n\n");
    for (int i = 0; i < 12; i++)
    {
        fprintf(file_linear_hist, "%g %g\n", bin_width_uniform/2.0 + i*bin_width_uniform, (double)linear_hist2[i]/N);
    }


    // Normal parameters
    const double mu = 0.2;
    const double sigma = 0.5;
    const double d = 1.0;
    const double xmin = mu - 3*sigma;
    const double xmax = mu + 3*sigma;
    double* normal_arr = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++)
    {
        normal_arr[i] = rand_normal(mu, sigma, d);
    }
    double normal_mu = find_mu(normal_arr, N);
    double normal_var = find_var(normal_arr, N);
    double normal_sigma = sqrt(normal_var);

    printf("Normal distribution\n");
    printf("mu value: %g\n", normal_mu);
    printf("variation value: %g\n", normal_var);
    printf("sigma value: %g\n\n", normal_sigma);

    //Normal histogram
    FILE* file_normal_hist = fopen("N_hist.dat", "w");
    int normal_hist[12];
    make_histogram(normal_arr, N, xmin, xmax, normal_hist, 12);

    const double bin_width_normal = (xmax-xmin) / 12.0;
    for (int i = 0; i < 12; i++)
    {
        fprintf(file_normal_hist, "%g %g\n", xmin + bin_width_normal/2.0 + i*bin_width_normal, (double)normal_hist[i]/N);
    }

    // Testing
    double xtest = 0.0;
    for (int j = 0; j < 12; j++)
    {   
        double pj = calc_pj(xmin + j*bin_width_normal, xmin + (j+1)*bin_width_normal, mu, sigma);
        printf("j = %d, pj = %g\n", j, pj);
        xtest += pow((normal_hist[j] - N*pj), 2) / (N * pj);
    }
    double conf = gammp(9.0f/2.0f, xtest/2.0f);
    double sign = 1.0 - conf;

    printf("\n");
    printf("X^2 = %g\n", xtest); // = 14.982
                                 // with nu = 9, eps = 16.91 for a = 0.05 => not rejected
    printf("confidence: %g\n", conf);
    printf("significance: %g\n", sign);


    free(linear_arr1);
    free(linear_arr2);
    free(normal_arr);
    fclose(file_linear);
    fclose(file_linear_hist);
    fclose(file_normal_hist);
    return 0;
}

double rand_linear(long a, long c, long m, int version)
{
    static long x_last1 = 10;
    static long x_last2 = 10;
    static long x_last3 = 10;

    if (version == 1)
    {
        x_last1 = (a*x_last1 + c) % m;
        return x_last1 / (m + 1.0);
    }
    else if (version == 2)
    {
        x_last2 = (a*x_last2 + c) % m;
        return x_last2 / (m + 1.0);
    }
    else
    {
        x_last3 = (a*x_last3 + c) % m;
        return x_last3 / (m + 1.0);
    }
    
}

double find_mu(double* numbers, int N)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        sum += numbers[i];
    }
    return sum / N;
}

double find_var(double* numbers, int N)
{
    double mu = find_mu(numbers, N);
    double ssum = 0.0;
    for (int i = 0; i < N; i++)
    {
        ssum += (numbers[i] - mu) * (numbers[i] - mu);
    }
    return ssum / N;
}

void make_histogram(double* numbers, int N, double min, double max, int* bins, int binN)
{
    for (int i = 0; i < binN; i++)
    {
        bins[i] = 0;
    }
    double bin_width = (max-min) / binN;
    for (int i = 0; i < N; i++)
    {
        int index = (numbers[i] - min) / bin_width;
        ++bins[index];
    }
}

double normal_density(double x, double mu, double sigma)
{
    return 1 / (sigma * sqrt(2*M_PI)) * exp(-1 * (x-mu)*(x-mu) / (2*sigma*sigma));
}

double rand_normal(double mu, double sigma, double d)
{
    const long a = 69069;
    const long c = 1;
    const long m = pow(2, 32);
    
    double xmin = mu - 3*sigma;
    double xmax = mu + 3*sigma;
    double u1, u2;
    do
    {
        u1 = rand_linear(a, c, m, 3);
        u2 = rand_linear(a, c, m, 3);
        u1 = u1 * (xmax - xmin) + xmin;
        u2 *= d;
    } while (u2 > normal_density(u1, mu, sigma));
    return u1;
}

double normal_distribution(double x, double mu, double sigma)
{
    return (1.0 + erf((x - mu) / (sqrt(2) * sigma))) / 2.0;
}

double calc_pj(double xmin, double xmax, double mu, double sigma)
{
    return normal_distribution(xmax, mu, sigma) - normal_distribution(xmin, mu, sigma);
}