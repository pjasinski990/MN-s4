#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"

#define N 1000
#define argM 5

void conjgrad(double** A, double* x, double* b, int size);
double scalar_prod(double* v, double* w, int n);
double* vector_sum_scaled(double* v, double* w, double scale, int n);
double* matrix_vector_prod(double** m, double* v, int n);

int main(int argc, const char *argv[]) {
	double** m = dmatrix(0, N-1, 0, N-1);
	double* b = dvector(0, N-1);
	double* x = dvector(0, N-1);

	for (int i = 0; i < N; i++) {
		x[i] = 0.0;
		b[i] = i+1.0;
		for (int j = 0; j < N; j++) {
			if (abs(i-j) <= argM) {
				m[i][j] = 1.0 / (1+abs(i-j));
			} else {
				m[i][j] = 0;
			}
		}
	}

	conjgrad(m, x, b, N);

	free_dmatrix(m, 0, N-1, 0, N-1);
	free_dvector(x, 0, N-1);
	free_dvector(b, 0, N-1);
	return 0;
}

void conjgrad(double** A, double* x, double* b, int size) {
	double* v = matrix_vector_prod(A, x, size);
	double* r = dvector(0, size-1);
	for (int i = 0; i < size; i++) {
		r[i] = v[i] = b[i] - v[i];
	}

	int k = 0;
	double alpha, beta, old_scalar_rr;
	alpha = beta = old_scalar_rr = 0.0;

	printf("%d %f %f %f %f\n", k, sqrt(scalar_prod(r, r, size)), alpha, beta, sqrt(scalar_prod(x, x, size)));

	while (sqrt(scalar_prod(r, r, size)) > 1e-6) {
		old_scalar_rr = scalar_prod(r, r, size);

		double* avprod = matrix_vector_prod(A, v, size);
		alpha = old_scalar_rr / scalar_prod(v, avprod, size);
		
		double* temp = vector_sum_scaled(x, v, alpha, size);
		for (int i = 0; i < size; i++) {
			x[i] = temp[i];		
		}
		free_dvector(temp, 0, size-1);

		temp = vector_sum_scaled(r, avprod, -alpha, size);
		free_dvector(r, 0, size-1);
		free_dvector(avprod, 0, size-1);
		r = temp;

		beta = scalar_prod(r, r, size)/old_scalar_rr;

		temp = vector_sum_scaled(r, v, beta, size);
		free_dvector(v, 0, size-1);
		v = temp;

		++k;
		printf("%d %f %f %f %f\n", k, sqrt(scalar_prod(r, r, size)), alpha, beta, sqrt(scalar_prod(x, x, size)));
	}
	free_dvector(v, 0, size-1);
	free_dvector(r, 0, size-1);
}

double scalar_prod(double* v, double* w, int n) {
	double temp = 0;
	for (int i = 0; i < n; i++) {
		temp += v[i]*w[i];
	}
	return temp;
}

double* matrix_vector_prod(double** m, double* v, int n) {
	double* result = dvector(0, n-1);

	for (int i = 0; i < n; i++) {
		result[i] = 0;
		for (int j = 0; j < n; j++) {
			result[i] += m[i][j]*v[j];
		}
	}
	return result;
}

double* vector_sum_scaled(double* v, double* w, double scale, int n) {
	double* result = dvector(0, n-1);
	for (int i = 0; i < n; i++) {
		result[i] = v[i] + w[i]*scale;
	}
	return result;
}
