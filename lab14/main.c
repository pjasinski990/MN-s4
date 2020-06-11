#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double rand_linear(long a, long c, long m, int version);
double find_mu(double* numbers, int N);
double find_var(double* numbers, int N);
void make_histogram(double* numbers, int N, double min, double max, int* bins, int binN);

int main(int argc, char const *argv[])
{
    const int N = 10000;

    long a1 = 123, c1 = 1, m1 = pow(2, 15);
    long a2 = 69069, c2 = 1, m2 = pow(2, 32);
    double* linear_arr1 = malloc(sizeof(double) * (N+1));
    double* linear_arr2 = malloc(sizeof(double) * (N+1));
    int* linear_hist1 = calloc(12, sizeof(int));
    int* linear_hist2 = calloc(12, sizeof(int));
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

    double mu1 = find_mu(linear_arr1+1, N);
    double mu2 = find_mu(linear_arr2+1, N);
    double sigma1 = sqrt(find_var(linear_arr1+1, N));
    double sigma2 = sqrt(find_var(linear_arr2+1, N));

    printf("Uniform distribution\n");
    printf("mu values: %g, %g\n", mu1, mu2);
    printf("sigma values: %g, %g\n", sigma1, sigma2);

    // Uniform histogram 
    make_histogram(linear_arr1+1, N, 0.0, 1.0, linear_hist1, 12);
    make_histogram(linear_arr2+1, N, 0.0, 1.0, linear_hist2, 12);
    double bin_width = 1.0 / 12.0;
    for (int i = 0; i < 12; i++)
    {
        fprintf(file_linear_hist, "%g %g\n", bin_width/2.0 + i*bin_width, (double)linear_hist1[i]/N);
    }
    fprintf(file_linear_hist, "\n\n");
    for (int i = 0; i < 12; i++)
    {
        fprintf(file_linear_hist, "%g %g\n", bin_width/2.0 + i*bin_width, (double)linear_hist2[i]/N);
    }


    // Normal parameters


    //Normal histogram

    free(linear_arr1);
    free(linear_arr2);
    free(linear_hist1);
    free(linear_hist2);
    fclose(file_linear);
    return 0;
}

double rand_linear(long a, long c, long m, int version)
{
    static long x_last1 = 10;
    static long x_last2 = 10;

    if (version == 1)
    {
        x_last1 = (a*x_last1 + c) % m;
        return x_last1 / (m + 1.0);
    }
    else
    {
        x_last2 = (a*x_last2 + c) % m;
        return x_last2 / (m + 1.0);
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
    double bin_width = (max-min) / binN;
    for (int i = 0; i < N; i++)
    {
        int index = (numbers[i] - min) / bin_width;
        ++bins[index];
    }
}