#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "mhd.h"
#include "jsw_rand.h"


static void random_beltrami_field(double x[4],
				  double u[4], int model, int k2, double h, int rank);


void mhd_initial_data_kh
(struct mhd_sim *sim, double x[4], double u[4], double b[4])
{
  double p = sim->status.time_simulation > 0.0 ? 0.0 : sim->user.pert;
  double a = sqrt(sim->user.k2);
  u[1] = 0.5 * (tanh((x[2]-0.25) * a) - tanh((x[2]-0.75) * a) - 1);
  u[2] = p * sin(4 * M_PI * x[1]);
  u[3] = 0.0;
}


void mhd_initial_data_abc
(struct mhd_sim *sim, double x[4], double u[4], double b[4])
{
  double a1 = sim->user.abc[0];
  double a2 = sim->user.abc[1];
  double a3 = sim->user.abc[2];
  double alpha = sqrt(sim->user.k2) * 2 * M_PI;
  
  u[1] = 0.0;
  u[2] = 0.0;
  u[3] = 0.0;

  u[2] += a1 * cos(alpha * x[1]);
  u[3] -= a1 * sin(alpha * x[1]);
  u[1] += 0.0;

  u[3] += a2 * cos(alpha * x[2]);
  u[1] -= a2 * sin(alpha * x[2]);
  u[2] += 0.0;

  u[1] += a3 * cos(alpha * x[3]);
  u[2] -= a3 * sin(alpha * x[3]);
  u[3] += 0.0;
}


void mhd_initial_data_beltrami
(struct mhd_sim *sim, double x[4], double u[4], double b[4])
{
  u[1] = 0.0;
  u[2] = 0.0;
  u[3] = 0.0;

  int rank = 0;
  rank += sim->user.N[0] > 1;
  rank += sim->user.N[1] > 1;
  rank += sim->user.N[2] > 1;

  random_beltrami_field(x, u, 0, sim->user.k2, sim->user.helicity, rank);
}



typedef complex double Complex;
typedef struct fourier_mode
{
  double k[4];
  Complex A[4];
} fourier_mode;




void random_beltrami_field
(double x[4], double B[4], int model, int k2, double h, int rank)
{
#define RAND jsw_random_double(&rand, -1, 1)
  int m,d,i,j,k;
  jsw_rand_t rand;
  jsw_seed(&rand, model);


  int k2_sphere = k2;
  int k_cube = floor(sqrt(k2_sphere));
  fourier_mode *modes = NULL;
  int num_modes = 0;
  Complex A[4] = {0, 0, 0, 0};


  for (i=-k_cube*(rank >= 1); i<=k_cube*(rank >= 1); ++i) {
    for (j=-k_cube*(rank >= 2); j<=k_cube*(rank >= 2); ++j) {
      for (k=-k_cube*(rank >= 3); k<=k_cube*(rank >= 3); ++k) {

	//printf("checking [%d %d %d]\n", i, j, k);

  	fourier_mode M;
	double phase = RAND * M_PI;


  	if (i*i + j*j + k*k != k2_sphere) {
  	  continue;
  	}
  	else {
	  
  	  //printf("k[%d] = [%d %d %d] is on shell\n", num_modes, i, j, k);

  	  M.k[0] = 0.0;
  	  M.k[1] = i;
  	  M.k[2] = j;
  	  M.k[3] = k;


	  double e0[4] = {0, RAND, RAND, RAND}; /* any vector not parallel to k */
	  double e1[4] = CROSS(M.k, e0);
	  double A1 = sqrt(DOT(e1, e1));

	  for (d=1; d<=3; ++d) {
	    e1[d] /= A1;
	  }

	  M.A[0] = 0.0;
	  for (d=1; d<=3; ++d) {	  
	    M.A[d] = e1[d] * cexp(I * phase);
	  }

  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;


	  /* reality condition */
	  for (d=1; d<=3; ++d) {
	    M.k[d] = -M.k[d];
	    M.A[d] = conj(M.A[d]);
	  }
  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;
	}
      }
    }
  }

  for (m=0; m<num_modes; ++m) {
    double a = sqrt(k2);
    fourier_mode M = modes[m];
    Complex K[4] = {0, I*M.k[1], I*M.k[2], I*M.k[3]};
    Complex Ikx  = (K[1]*x[1] + K[2]*x[2] + K[3]*x[3]) * 2 * M_PI;
    Complex P[4] = {0, M.A[1], M.A[2], M.A[3]}; /* a times psi */
    Complex H[4] = CROSS(K, P);

    A[1] += (P[1] + h*H[1]/a) * cexp(Ikx);
    A[2] += (P[2] + h*H[2]/a) * cexp(Ikx);
    A[3] += (P[3] + h*H[3]/a) * cexp(Ikx);
  }


  free(modes);

  /* printf("re(A) = %+8.6e %+8.6e %+8.6e\n", */
  /* 	 creal(A[1]), creal(A[2]), creal(A[3])); */
  /* printf("im(A) = %+8.6e %+8.6e %+8.6e\n", */
  /* 	 cimag(A[1]), cimag(A[2]), cimag(A[3])); */

  B[1] = A[1];
  B[2] = A[2];
  B[3] = A[3];
  
#undef RAND
}
