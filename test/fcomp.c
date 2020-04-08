#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../header/wavelib.h"

int main() {
	int i, N, J,subscale,a0,iter,nd,k, seq;
	double dt, dj,s0, param,mn;
	cwt_object wt;
	double voltage, fft_out, sum_all_alpha;
	int sample = 0;

	FILE *ifp;

	char *wave = "morlet";// Set Morlet wavelet. Other options "paul" and "dog"
	char *type = "pow";

	N = 256;
	param = 6.0;
	subscale = 4;
	dt = 0.0025;
	s0 = dt;
	dj = 0.125; //1.0 / (double)subscale;
	J = 80; //11 * subscale; // Total Number of scales
	a0 = 2;//power
	nd = N/2 - 1;
	double power_array[N];

	ifp = fopen("x.data", "r");
	i = 0;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}

	wt = cwt_init(wave, param, N,dt, J);

	setCWTScales(wt, s0, dj, type, a0);

	while (fscanf(ifp, "%d %lf %lf\n", &seq, &voltage, &fft_out)) {
		sample++;
		if (sample > N - 1) {
			memmove(&power_array[0], &power_array[1], sizeof(double) * (N - 1));
			sample = N - 1;
		}
		power_array[sample] = voltage;

		if (sample == N - 1) {
			cwt(wt, power_array);
		//	cwt_summary(wt);
			sum_all_alpha = 0;
			for (int k = 40; k <= 45; k++) {	// just the alpha scales
				iter = nd + k * N;
				sum_all_alpha += wt->output[iter].re * wt->output[iter].re + wt->output[iter].im * wt->output[iter].im;
			}
			printf("%d,%lf,%lf,%lf\n", seq, sum_all_alpha, fft_out, voltage);
		}
	}
	cwt_free(wt);
	fclose(ifp);
	return 0;
}
