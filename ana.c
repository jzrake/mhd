#include <stdio.h>
#include <math.h>
#include "mhd.h"



int mhd_status_fscanf(struct mhd_status *stat, FILE *F)
{
  return fscanf(F, "%lf "
		"%lf %lf %lf %lf "
		"%lf %lf %lf %lf "
		"\n",
		&stat->time_simulation,
		&stat->velocity_energy,
		&stat->velocity_helicity,
		&stat->velocity_monopole,
		&stat->velocity_L2,
		&stat->magnetic_energy,
		&stat->magnetic_helicity,
		&stat->magnetic_monopole,
		&stat->magnetic_L2);
}



int mhd_status_fprintf(struct mhd_status *stat, FILE *F)
{
  return fprintf(F, "%+12.10e "
		 "%+12.10e %+12.10e %+12.10e %+12.10e "
		 "%+12.10e %+12.10e %+12.10e %+12.10e "
		 "\n",
		 stat->time_simulation,
		 stat->velocity_energy,
		 stat->velocity_helicity,
		 stat->velocity_monopole,
		 stat->velocity_L2,
		 stat->magnetic_energy,
		 stat->magnetic_helicity,
		 stat->magnetic_monopole,
		 stat->magnetic_L2);
}



/*
 * Truncate log files on restart
 * =====================================================================
 */
void mhd_truncate_logfile(double t, const char *fname)
{
  FILE *logf = fopen(fname, "r");

  if (logf == NULL) {
    printf("[mhd] warning: no existing mhd.dat file\n");
    return;
  }

  size_t S = sizeof(struct mhd_status);
  int N = 1;
  struct mhd_status *measure = (struct mhd_status *) malloc(S);

  while (mhd_status_fscanf(&measure[N-1], logf) != EOF) {
    measure = (struct mhd_status *) realloc(measure, ++N * S);
  }

  fclose(logf);
  logf = fopen(fname, "w");

  for (int n=0; n<N-1; ++n) {
    if (measure[n].time_simulation <= t) {
      mhd_status_fprintf(&measure[n], logf);
    }
  }

  fclose(logf);
  free(measure);
}



/*
 * Carry out measurement diagnostic
 * =====================================================================
 */
void mhd_sim_measure
(struct mhd_sim *sim, struct mhd_status *stat)
{
#define GLB_AVG(x) x = cow_domain_dblsum(sim->domain, x) / Nt

  long long Nt = cow_domain_getnumglobalzones(sim->domain, COW_ALL_DIMS);
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];

  double *u = cow_dfield_getdatabuffer(sim->velocity[0]);
  double *b = cow_dfield_getdatabuffer(sim->magnetic[0]);

  stat->time_simulation = sim->status.time_simulation;
  stat->velocity_energy = 0.0;
  stat->magnetic_energy = 0.0;
  stat->velocity_monopole = 0.0;
  stat->magnetic_monopole = 0.0;
  stat->velocity_L2 = 0.0;
  stat->magnetic_L2 = 0.0;

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    double divu = D1(u,1) + D2(u,2) + D3(u,3);
    double divb = D1(b,1) + D2(b,2) + D3(b,3);

    double uu = DOT(&u[m], &u[m]);
    double bb = DOT(&b[m], &b[m]);

    stat->velocity_energy += 0.5 * uu;
    stat->magnetic_energy += 0.5 * bb;
    stat->velocity_monopole += fabs(divu);
    stat->magnetic_monopole += fabs(divb);
    
    double x[4] = {0,
    		   cow_domain_positionatindex(sim->domain, 0, i),
    		   cow_domain_positionatindex(sim->domain, 1, j),
    		   cow_domain_positionatindex(sim->domain, 2, k)};

    double u0[4];
    double b0[4];

    sim->initial_data(sim, x, u0, b0);

    double du[4] = {0, u[m+1] - u0[1], u[m+2] - u0[2], u[m+3] - u0[3]};
    double db[4] = {0, b[m+1] - b0[1], b[m+2] - b0[2], b[m+3] - b0[3]};

    stat->velocity_L2 += DOT(du, du);
    stat->magnetic_L2 += DOT(db, db); 
  }
  
  GLB_AVG(stat->velocity_energy);
  GLB_AVG(stat->magnetic_energy);
  GLB_AVG(stat->velocity_monopole);
  GLB_AVG(stat->magnetic_monopole);
  GLB_AVG(stat->velocity_L2);
  GLB_AVG(stat->magnetic_L2);
#undef GLB_AVG
}


void mhd_sim_analyze
(struct mhd_sim *sim, struct mhd_status *stat, char *filename)
{
  cow_domain *domain = sim->domain;
  cow_dfield *velocity = sim->velocity[0];
  cow_dfield *magnetic = sim->magnetic[0];
  cow_dfield *vecpoten = cow_dfield_new();
  cow_dfield *vorticit = cow_dfield_new();

  char gname[1024];
  char nname[1024];


  snprintf(gname, 1024, "spectra-%06d", sim->status.iteration);
  snprintf(nname, 1024, "%12.10e", sim->status.time_simulation);


  /* Data fields setup */
  /* ---------------------------------------------------- */
  cow_dfield_setdomain(vecpoten, domain);
  cow_dfield_setname(vecpoten, "vector_potential");
  cow_dfield_addmember(vecpoten, "a1");
  cow_dfield_addmember(vecpoten, "a2");
  cow_dfield_addmember(vecpoten, "a3");
  cow_dfield_commit(vecpoten);

  cow_dfield_setdomain(vorticit, domain);
  cow_dfield_setname(vorticit, "vorticity");
  cow_dfield_addmember(vorticit, "w1");
  cow_dfield_addmember(vorticit, "w2");
  cow_dfield_addmember(vorticit, "w3");
  cow_dfield_commit(vorticit);



  /* Histograms setup */
  /* ---------------------------------------------------- */
  cow_histogram *Pu = cow_histogram_new();
  cow_histogram_setlower(Pu, 0, 1);
  cow_histogram_setupper(Pu, 0, sim->user.max_pspec_bin);
  cow_histogram_setnbins(Pu, 0, sim->user.num_pspec_bin);
  cow_histogram_setspacing(Pu, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Pu, "velocity");
  cow_histogram_setfullname(Pu, nname);


  cow_histogram *Pb = cow_histogram_new();
  cow_histogram_setlower(Pb, 0, 1);
  cow_histogram_setupper(Pb, 0, sim->user.max_pspec_bin);
  cow_histogram_setnbins(Pb, 0, sim->user.num_pspec_bin);
  cow_histogram_setspacing(Pb, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Pb, "magnetic");
  cow_histogram_setfullname(Pb, nname);


  cow_histogram *Hu = cow_histogram_new();
  cow_histogram_setlower(Hu, 0, 1);
  cow_histogram_setupper(Hu, 0, sim->user.max_pspec_bin);
  cow_histogram_setnbins(Hu, 0, sim->user.num_pspec_bin);
  cow_histogram_setspacing(Hu, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Hu, "velocity_helicity");
  cow_histogram_setfullname(Hu, nname);

  
  cow_histogram *Hb = cow_histogram_new();
  cow_histogram_setlower(Hb, 0, 1);
  cow_histogram_setupper(Hb, 0, sim->user.max_pspec_bin);
  cow_histogram_setnbins(Hb, 0, sim->user.num_pspec_bin);
  cow_histogram_setspacing(Hb, COW_HIST_SPACING_LINEAR);
  cow_histogram_setnickname(Hb, "magnetic_helicity");
  cow_histogram_setfullname(Hb, nname);

  cow_fft_curl       (velocity, vorticit);
  cow_fft_inversecurl(magnetic, vecpoten);

  cow_fft_pspecvecfield(velocity, Pu);
  cow_fft_pspecvecfield(magnetic, Pb);
  cow_fft_helicityspec(magnetic, Hu, 'u');
  cow_fft_helicityspec(magnetic, Hb, 'b');

  long long Nt = cow_domain_getnumglobalzones(sim->domain, COW_ALL_DIMS);
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);
  double *u = (double*) cow_dfield_getdatabuffer(velocity);
  double *w = (double*) cow_dfield_getdatabuffer(vorticit);
  double *a = (double*) cow_dfield_getdatabuffer(vecpoten);
  double *b = (double*) cow_dfield_getdatabuffer(magnetic);

  double uutot = 0.0;
  double bbtot = 0.0;
  double uwtot = 0.0;
  double abtot = 0.0;
  
  FOR_ALL_INTERIOR_NO_THREAD(Ni, Nj, Nk) {
    int m = INDV(i,j,k);

    double uu = DOT(&u[m], &u[m]);
    double bb = DOT(&b[m], &b[m]);
    double uw = DOT(&u[m], &w[m]);
    double ab = DOT(&a[m], &b[m]);

    uutot += uu;
    bbtot += bb;
    uwtot += uw;
    abtot += ab;
  }

  uutot = cow_domain_dblsum(domain, uutot) / Nt;
  bbtot = cow_domain_dblsum(domain, bbtot) / Nt;
  uwtot = cow_domain_dblsum(domain, uwtot) / Nt;
  abtot = cow_domain_dblsum(domain, abtot) / Nt;


  if (stat) {
    stat->velocity_helicity = uwtot;
    stat->magnetic_helicity = abtot;
  }


  if (filename) {

    cow_histogram_dumphdf5(Pu, filename, gname);
    cow_histogram_dumphdf5(Pb, filename, gname);
    cow_histogram_dumphdf5(Hu, filename, gname);
    cow_histogram_dumphdf5(Hb, filename, gname);

    if (0) { /* write derived fields */
      cow_dfield_write(vecpoten, filename);
      cow_dfield_write(vorticit, filename);
    }

  }

  cow_histogram_del(Pu);
  cow_histogram_del(Pb);
  cow_histogram_del(Hu);
  cow_histogram_del(Hb);

  cow_dfield_del(vecpoten);
  cow_dfield_del(vorticit);
}
