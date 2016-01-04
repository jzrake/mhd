#ifndef MHD_HEADER
#define MHD_HEADER

#ifndef M_PI
#define M_PI 3.1415926535897926
#endif

#define MHD_NG 3		/* number of guard zones */
#define MHD_DIFFERENCE_ORDER 4
#define MHD_DISSIPATION_ORDER 6

#include <stdio.h> /* FILE */
#include "cow/cow.h"
#include "ser.h"


struct mhd_sim;
typedef void (*InitialDataFunc)
(struct mhd_sim *sim, double x[4], double u[4], double b[4]);


struct mhd_sim {
  struct mhd_user user;
  struct mhd_status status;
  InitialDataFunc initial_data;
  char problem_name[1024];
  cow_domain *domain;
  cow_dfield *velocity[6];
  cow_dfield *magnetic[6];
  double grid_spacing[4];
};


/* ini.c */
void mhd_initial_data_kh
(struct mhd_sim *sim, double x[4], double u[4], double b[4]);
void mhd_initial_data_abc
(struct mhd_sim *sim, double x[4], double u[4], double b[4]);
void mhd_initial_data_beltrami
(struct mhd_sim *sim, double x[4], double u[4], double b[4]);


/* ana.c */
int mhd_status_fscanf(struct mhd_status *stat, FILE *F);
int mhd_status_fprintf(struct mhd_status *stat, FILE *F);
void mhd_truncate_logfile(double t, const char *fname);
void mhd_sim_measure(struct mhd_sim *sim, struct mhd_status *stat);
void mhd_sim_analyze(struct mhd_sim *sim, struct mhd_status *stat, char *filename);


/*
 * Macro for a three-dimensional loop over all interior cells
 * =====================================================================
 */
#define FOR_ALL_INTERIOR(N1, N2, N3)				\
  for (int i=N1==1?0:MHD_NG; i<N1+(N1==1?0:MHD_NG); ++i)	\
    for (int j=N2==1?0:MHD_NG; j<N2+(N2==1?0:MHD_NG); ++j)	\
      for (int k=N3==1?0:MHD_NG; k<N3+(N3==1?0:MHD_NG); ++k)	\



/*
 * Macro to calculate linear index of (i,j,k,m) ... m goes from 1, not 0
 * =====================================================================
 */

#define INDV(i,j,k) ((i)*si + (j)*sj + (k)*sk - 1)
#define INDS(i,j,k) ((i)*ti + (j)*tj + (k)*tk)	/* for scalar field */



/*
 * Macros to calculate finite dimhdrences
 * =====================================================================
 */
#define DIFF1C2(F,s) ((-1*(F)[-1*s] +		\
		       +1*(F)[+1*s]) / 2.0)

#define DIFF1C4(F,s) ((+1*(F)[-2*s] +		\
		       -8*(F)[-1*s] +		\
		       +8*(F)[+1*s] +		\
		       -1*(F)[+2*s]) / 12.0)

#define DIFF2C2(F,s) ((+1*(F)[-1*s] +		\
		       -2*(F)[+0*s] +		\
		       +1*(F)[+1*s]) / 1.0)

#define DIFF2C4(F,s) ((-1 *(F)[-2*s] +		\
		       +16*(F)[-1*s] +		\
		       -30*(F)[+0*s] +		\
		       +16*(F)[+1*s] +		\
		       -1 *(F)[+2*s]) / 12.0)

#define DIFF4C2(F,s) ((+1*(F)[-2*s] +		\
		       -4*(F)[-1*s] +		\
		       +6*(F)[+0*s] +		\
		       -4*(F)[+1*s] +		\
		       +1*(F)[+2*s]) / 1.0)

#define DIFF4C4(F,s) ((-1 *(F)[-3*s] +		\
		       +12*(F)[-2*s] +		\
		       -39*(F)[-1*s] +		\
		       +56*(F)[+0*s] +		\
		       -39*(F)[+1*s] +		\
		       +12*(F)[+2*s] +		\
		       -1 *(F)[+3*s]) / 6.0)

#define DIFF6C2(F,s) (( 1 *(F)[-3*s] +		\
		       -6 *(F)[-2*s] +		\
		       +15*(F)[-1*s] +		\
		       -20*(F)[+0*s] +		\
		       +15*(F)[+1*s] +		\
		       -6 *(F)[+2*s] +		\
		        1 *(F)[+3*s]) / 1.0)


#define CROSS(E,B) {0.0,				\
      (E)[2]*(B)[3]-(E)[3]*(B)[2],			\
      (E)[3]*(B)[1]-(E)[1]*(B)[3],			\
      (E)[1]*(B)[2]-(E)[2]*(B)[1]}			\

#define DOT(E,B) ((E)[1]*(B)[1] + (E)[2]*(B)[2] + (E)[3]*(B)[3])

#define MAX3(a,b,c) (a>=b && b>=c ? a : (b >= c ? b : c))
#define MIN3(a,b,c) (a<=b && b<=c ? a : (b <= c ? b : c))


/* https://einsteintoolkit.org/documentation/ThornDoc/CactusNumerical/Dissipation */
#if (MHD_DISSIPATION_ORDER == 2)
#define KOV(F,c) ((Ni==1 ? 0.0 : DIFF2C2(F+m+c,si)) +	\
		  (Nj==1 ? 0.0 : DIFF2C2(F+m+c,sj)) +	\
		  (Nk==1 ? 0.0 : DIFF2C2(F+m+c,sk)))
#define KOS(F  ) ((Ni==1 ? 0.0 : DIFF2C2(F+n  ,ti)) +	\
		  (Nj==1 ? 0.0 : DIFF2C2(F+n  ,tj)) +	\
		  (Nk==1 ? 0.0 : DIFF2C2(F+n  ,tk)))
#elif (MHD_DISSIPATION_ORDER == 4)
#define KOV(F,c) ((Ni==1 ? 0.0 : DIFF4C2(F+m+c,si)) +	\
		  (Nj==1 ? 0.0 : DIFF4C2(F+m+c,sj)) +	\
		  (Nk==1 ? 0.0 : DIFF4C2(F+m+c,sk)))
#define KOS(F  ) ((Ni==1 ? 0.0 : DIFF4C2(F+n  ,ti)) +	\
		  (Nj==1 ? 0.0 : DIFF4C2(F+n  ,tj)) +	\
		  (Nk==1 ? 0.0 : DIFF4C2(F+n  ,tk)))
#elif (MHD_DISSIPATION_ORDER == 6)
#define KOV(F,c) ((Ni==1 ? 0.0 : DIFF6C2(F+m+c,si)) +	\
		  (Nj==1 ? 0.0 : DIFF6C2(F+m+c,sj)) +	\
		  (Nk==1 ? 0.0 : DIFF6C2(F+m+c,sk)))
#define KOS(F  ) ((Ni==1 ? 0.0 : DIFF6C2(F+n  ,ti)) +	\
		  (Nj==1 ? 0.0 : DIFF6C2(F+n  ,tj)) +	\
		  (Nk==1 ? 0.0 : DIFF6C2(F+n  ,tk)))
#endif


#if (MHD_DIFFERENCE_ORDER == 2)

#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C2(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C2(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C2(F+m+c,sk)/dz)

#define S1(F  )  (Ni==1 ? 0.0 : DIFF1C2(F+n+0,ti)/dx)
#define S2(F  )  (Nj==1 ? 0.0 : DIFF1C2(F+n+0,tj)/dy)
#define S3(F  )  (Nk==1 ? 0.0 : DIFF1C2(F+n+0,tk)/dz)

#elif (MHD_DIFFERENCE_ORDER == 4)

#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C4(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C4(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C4(F+m+c,sk)/dz)

#define S1(F  )  (Ni==1 ? 0.0 : DIFF1C4(F+n+0,ti)/dx)
#define S2(F  )  (Nj==1 ? 0.0 : DIFF1C4(F+n+0,tj)/dy)
#define S3(F  )  (Nk==1 ? 0.0 : DIFF1C4(F+n+0,tk)/dz)

#endif



#endif				/* MHD_HEADER */
