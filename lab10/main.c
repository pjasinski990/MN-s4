#include <stdio.h>
#include <math.h>

double f(double x, double y);
double dfdx(double x, double y, double dx);
double dfdy(double x, double y, double dx);

int main(int argc, char const *argv[])
{
    const double eps1 = 1.0e-2;
    const double eps2 = 1.0e-3;
    const int it_max = 1000;
    const double h = 0.1;
    const double dx = 1.0e-4;

    double x0 = -0.75;
    double y0 = 1.75;
    double x1, y1;

    FILE* f_eps1 = fopen("eps1.dat", "w");
    FILE* f_eps2 = fopen("eps2.dat", "w");

    fprintf(f_eps1, "%g %g\n", x0, y0);
    for (int i = 0; i < it_max; i++)
    {
        x1 = x0 - h * dfdx(x0, y0, dx);
        y1 = y0 - h * dfdy(x0, y0, dx);
        fprintf(f_eps1, "%g %g\n", x1, y1);

        if (sqrt(pow((x1-x0), 2) + pow((y1-y0), 2)) < eps1)
        {
            break;
        }
        x0 = x1;
        y0 = y1;
    }

    x0 = -0.75;
    y0 = 1.75;
    fprintf(f_eps2, "%g %g\n", x0, y0);
    for (int i = 0; i < it_max; i++)
    {
        x1 = x0 - h * dfdx(x0, y0, dx);
        y1 = y0 - h * dfdy(x0, y0, dx);
        fprintf(f_eps2, "%g %g\n", x1, y1);

        if (sqrt(pow((x1-x0), 2) + pow((y1-y0), 2)) < eps2)
        {
            break;
        }
        x0 = x1;
        y0 = y1;
    }

    fclose(f_eps1);
    fclose(f_eps2);
    return 0;
}

double f(double x, double y)
{
    return 5.0/2.0 * (x*x - y)*(x*x - y) + (1 - x)*(1 - x);
}

double dfdx(double x, double y, double dx)
{
    return (f(x + dx, y) - f(x - dx, y)) / (2.0*dx);
}

double dfdy(double x, double y, double dy)
{
    return (f(x, y + dy) - f(x, y - dy)) / (2.0*dy);
}
