#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>

double rho(double alpha, int x);
int kdelta(int i, int j);

const int n = 200;
const double N = 1.0;
const double L = 10.0;

int main(int argc, const char *argv[]) {

	const double dx = L/(n+1);

	gsl_matrix* mA = gsl_matrix_calloc(n, n);
	gsl_matrix* mB = gsl_matrix_calloc(n, n);
	gsl_matrix* evec = gsl_matrix_calloc(n, n);
	gsl_vector* eval = gsl_vector_calloc(n);
	gsl_eigen_gensymmv_workspace* w = gsl_eigen_gensymmv_alloc(n);

	FILE* evec_f = fopen("evec.dat", "w");
	FILE* eval_f = fopen("eval.dat", "w");

	for (int alpha = 0; alpha <= 100; alpha+=2) {
		//Filling matrices
		for (int i = 0; i < n; i++) {
			double x = -L/2.0 + dx*(i+1);
			for (int j = 0; j < n; j++) {
				double vala = (-1*kdelta(i, j+1) + 2*kdelta(i, j) - kdelta(i, j-1)) / (dx*dx);
				double valb = rho(x, alpha)/N * kdelta(i, j);
				gsl_matrix_set(mA, i, j, vala);
				gsl_matrix_set(mB, i, j, valb);
			}
		}

		gsl_eigen_gensymmv(mA, mB, eval, evec, w);
		gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		//Printing eval values
		fprintf(eval_f, "%d ", alpha);
		for (int i = 0; i < 6; i++) {
			fprintf(eval_f, "%lf ", sqrt(gsl_vector_get(eval, i)));
		}
		fprintf(eval_f, "\n");

		//Printing evec values
		if (alpha == 0 || alpha == 100) {
			for (int i = 0; i < n; i++) {
				double x = -L/2.0 + dx*(i+1);
				fprintf(evec_f, "%lf ", x);

				for (int j = 0; j < 6; j++) {
					fprintf(evec_f, "%lf ", gsl_matrix_get(evec, i, j));
				}
				fprintf(evec_f, "\n");
			}
			if (alpha == 0) {
				fprintf(evec_f, "\n");
				fprintf(evec_f, "\n");
			}
		}
	}
	gsl_matrix_free(mA);
	gsl_matrix_free(mB);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);
	gsl_eigen_gensymmv_free(w);

	fclose(evec_f);
	fclose(eval_f);
	return 0;
}

double rho(double x, int alpha) {
	return 4*alpha*x*x + 1.0;
}

int kdelta(int i, int j) {
	return i==j? 1:0;
}
