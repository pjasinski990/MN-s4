#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"

#define IT_MAX 12

void conjgrad(float** A, float* x, float* b, int size);
float scalar_prod(float* v, float* w, int n);
float* matrix_vector_prod(float** m, float* v, int n);
float** matrix_matrix_prod(float** m1, float** m2, int n);
float* vector_sum_scaled(float* v, float* w, float scale, int n);
void gs_ort(float* v, float* u, int n);
void fill_with_ones(float* vec, int n);
float* get_column(float** m, int index, int n);

int main(int argc, const char *argv[]) {
	const int n = 7;
	const int k_val = n;
	
	float** mA = matrix(0, n-1, 0, n-1);
	float** mX = matrix(0, n-1, 0, n-1);
	float** mXT = matrix(0, n-1, 0, n-1);
	float* x_k = vector(0, n-1);
	float* x_k1 = vector(0, n-1);

	FILE* lambda_f = fopen("lambda.dat", "w");

	// Filling matrix
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			mA[i][j] = 1/sqrt(2+abs(i-j));
		}
	}

	// Loop
	for (int k = 0; k < k_val; k++) {
		fill_with_ones(x_k, n);

		for (int i = 1; i <= IT_MAX; i++) {
			free_vector(x_k1, 0, n-1);
			x_k1 = matrix_vector_prod(mA, x_k, n);

			for (int j = 0; j < k; j++) {
				gs_ort(x_k1, mXT[j], n);
			}

			float lambda = scalar_prod(x_k1, x_k, n)/scalar_prod(x_k, x_k, n);
			fprintf(lambda_f, "%d %f\n", i, lambda);

			float norm = sqrt(scalar_prod(x_k1, x_k1, n));
			for (int i = 0; i < n; i++) {
				x_k[i] = x_k1[i]/norm;
			}
		}
		fprintf(lambda_f, "\n\n");

		for (int i = 0; i < n; i++) {
			mXT[k][i] = x_k[i];
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			mX[i][j] = mXT[j][i];
		}
	}


	float** temp = matrix_matrix_prod(mXT, mA, n);
	float** mD = matrix_matrix_prod(temp, mX, n);

	// Print mD matrix
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%f ", mD[i][j]);
		}
		printf("\n");
	}

	free_matrix(temp, 0, n-1, 0, n-1);

	fclose(lambda_f);
	free_matrix(mA, 0, n-1, 0, n-1);
	free_matrix(mX, 0, n-1, 0, n-1);
	free_matrix(mXT, 0, n-1, 0, n-1);
	free_matrix(mD, 0, n-1, 0, n-1);
	free_vector(x_k, 0, n-1);
	free_vector(x_k1, 0, n-1);
	return 0;
}

float scalar_prod(float* v, float* w, int n) {
	float temp = 0;
	for (int i = 0; i < n; i++) {
		temp += v[i]*w[i];
	}
	return temp;
}

float* matrix_vector_prod(float** m, float* v, int n) {
	float* result = vector(0, n-1);

	for (int i = 0; i < n; i++) {
		result[i] = 0;
		for (int j = 0; j < n; j++) {
			result[i] += m[i][j]*v[j];
		}
	}
	return result;
}

float* vector_sum_scaled(float* v, float* w, float scale, int n) {
	float* result = vector(0, n-1);
	for (int i = 0; i < n; i++) {
		result[i] = v[i] + w[i]*scale;
	}
	return result;
}

void gs_ort(float* v, float* u, int n) {
	float* temp = vector_sum_scaled(v, u, -1*scalar_prod(v, u, n), n);
	for (int i = 0; i < n; i++) {
		v[i] = temp[i];
	}
	free_vector(temp, 0, n-1);
}

void fill_with_ones(float* vec, int n) {
	for (int i = 0; i < n; i++) {
		vec[i] = 1;
	}
}

float* get_column(float** m, int index, int n) {
	float* res = vector(0, n-1);
	for (int i = 0; i < n; i++) {
		res[i] = m[i][index];
	}
	return res;
}

float** matrix_matrix_prod(float** m1, float** m2, int n) {
	float** res = matrix(0, n-1, 0, n-1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0;
			for (int k = 0; k < n; k++) {
				sum += m1[i][k]*m2[k][j];
			}
			res[i][j] = sum;
		}
	}
	return res;
}
