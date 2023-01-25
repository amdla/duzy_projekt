#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "makespl.h"
#include "gaus/piv_ge_solver.h"

double laguerre(int n, double x) {
	if (n == 0) {
		return 1;
	}
	if (n == 1) {
		return -x+1;
	}
	return ((2*(n-1)+1-x) * laguerre(n-1, x) - (n-1) * laguerre(n-2, x)) / n; //przesuwamy indeksy h(x) o jeden w lewo, stąd n - 1 ; n - 2
}

double ld1(int n, double x) {
	if (n == 0) {
		return 0;
	}
	return ld1(n-1, x) - laguerre(n-1, x);
}

double ld2(int n, double x) {
	if (n == 0) {
		return 0;
	}
	return ld2(n-1, x) - ld1(n-1, x);
}

double ld3(int n, double x) {
	if (n == 0) {
		return 0;
	}     
	return ld3(n-1, x) - ld2(n-1, x);
}


void make_spl(points_t* pts, spline_t* spl) {
	matrix_t* eqs = NULL;
	double* x = pts->x;
	double* y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1]; 
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
	char* nbEnv = getenv("APPROX_BASE_SIZE");

	if (nbEnv != NULL && atoi(nbEnv) > 0)
		nb = atoi(nbEnv);

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, laguerre(i, x[k]) * laguerre(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * laguerre(j, x[k]));
	}

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;	
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i * (b - a) / (spl->n - 1);
			xx += 10.0 * DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i] += ck * laguerre(k, xx);
				spl->f1[i] += ck * ld1(k, xx);
				spl->f2[i] += ck * ld2(k, xx);
				spl->f3[i] += ck * ld3(k, xx);
			}
		}
	}
}