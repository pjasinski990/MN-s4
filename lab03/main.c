#include <stdio.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c"

#define N 1000
#define argM 5

void conjgrad(float** A, float* x, float* b, int size);
float scalar_prod(float* v, float* w, int n);
float* vector_sum_scaled(float* v, float* w, float scale, int n);
float* matrix_vector_prod(float** m, float* v, int n);

int main(int argc, const char *argv[]) {
	float** m = matrix(1, N, 1, N);
	float* b = vector(1, N);
	float* x = vector(1, N);

	for (int i = 1; i <= N; i++) {
		x[i] = 0.0;
		b[i] = i;
		for (int j = 1; j <= N; j++) {
			if (abs(i-j) <= argM) {
				m[i][j] = 1.0 / (1+abs(i-j));
			} else {
				m[i][j] = 0;
			}
		}
	}

	// Calculating times
	clock_t t1g, t2g, t1cg, t2cg;
	float tg, tcg;

	t1cg = clock();
	conjgrad(m, x, b, N);
	t2cg = clock();
	tcg = ((float)(t2cg-t1cg))/CLOCKS_PER_SEC;
	// Comment since there to have proper output for graphing
	printf("tcg is: %lf\n", tcg);

	float** tempM = matrix(1, N, 1, N);
	t1g = clock();
	gaussj(m, N, tempM, 1);
	t2cg = clock();
	tcg = ((float)(t2cg-t1cg))/CLOCKS_PER_SEC;
	printf("tcg is: %lf\n", tcg);
	free_matrix(tempM, 1, N, 1, N);
	// Comment till there to have proper output for graphing

	free_matrix(m, 1, N, 1, N);
	free_vector(x, 1, N);
	free_vector(b, 1, N);
	return 0;
}

void conjgrad(float** A, float* x, float* b, int size) {
	float* v = matrix_vector_prod(A, x, size);
	float* r = vector(1, size);
	for (int i = 1; i <= size; i++) {
		r[i] = v[i] = b[i] - v[i];
	}

	int k = 0;
	float alpha, beta, old_scalar_rr;
	alpha = beta = old_scalar_rr = 0.0;

	printf("%d %f %f %f %f\n", k, sqrt(scalar_prod(r, r, size)), alpha, beta, sqrt(scalar_prod(x, x, size)));

	while (sqrt(scalar_prod(r, r, size)) > 1e-6) {
		old_scalar_rr = scalar_prod(r, r, size);

		float* avprod = matrix_vector_prod(A, v, size);
		alpha = old_scalar_rr / scalar_prod(v, avprod, size);
		
		float* temp = vector_sum_scaled(x, v, alpha, size);
		for (int i = 1; i <= size; i++) {
			x[i] = temp[i];		
		}
		free_vector(temp, 1, size);

		temp = vector_sum_scaled(r, avprod, -alpha, size);
		free_vector(r, 1, size);
		free_vector(avprod, 1, size);
		r = temp;

		beta = scalar_prod(r, r, size)/old_scalar_rr;

		temp = vector_sum_scaled(r, v, beta, size);
		free_vector(v, 1, size);
		v = temp;

		++k;
		printf("%d %f %f %f %f\n", k, sqrt(scalar_prod(r, r, size)), alpha, beta, sqrt(scalar_prod(x, x, size)));
	}
	free_vector(v, 1, size);
	free_vector(r, 1, size);
}

float scalar_prod(float* v, float* w, int n) {
	float temp = 0;
	for (int i = 1; i <= n; i++) {
		temp += v[i]*w[i];
	}
	return temp;
}

float* matrix_vector_prod(float** m, float* v, int n) {
	float* result = vector(1, n);

	for (int i = 1; i <= n; i++) {
		result[i] = 0;
		for (int j = 1; j <= n; j++) {
			result[i] += m[i][j]*v[j];
		}
	}
	return result;
}

float* vector_sum_scaled(float* v, float* w, float scale, int n) {
	float* result = vector(1, n);
	for (int i = 1; i <= n; i++) {
		result[i] = v[i] + w[i]*scale;
	}
	return result;
}
