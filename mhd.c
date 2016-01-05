#include <stdio.h>
#include <string.h>
#include <sys/stat.h> /* mkdir */
#include "jsw_rand.h"
#include "mhd.h"



#define DEL2(F,c) ((Ni==1 ? 0.0 : DIFF2C2(F+m+c,si))/(dx*dx) +	\
		   (Nj==1 ? 0.0 : DIFF2C2(F+m+c,sj))/(dy*dy) +	\
		   (Nk==1 ? 0.0 : DIFF2C2(F+m+c,sk))/(dz*dz))


/*
 * Sim initializer
 * =====================================================================
 */
void mhd_sim_init(struct mhd_sim *sim)
{
  int *N = sim->user.N;
  sim->domain = cow_domain_new();

  sim->grid_spacing[0] = 0.0;
  sim->grid_spacing[1] = 1.0 / N[0];
  sim->grid_spacing[2] = 1.0 / N[1];
  sim->grid_spacing[3] = 1.0 / N[2];

  cow_domain_setndim(sim->domain, (N[0]>1) + (N[1]>1) + (N[2]>1));
  cow_domain_setsize(sim->domain, 0, N[0]);
  cow_domain_setsize(sim->domain, 1, N[1]);
  cow_domain_setsize(sim->domain, 2, N[2]);
  cow_domain_setguard(sim->domain, MHD_NG);
  cow_domain_commit(sim->domain);
 


  for (int n=0; n<6; ++n) {
    sim->velocity[n] = cow_dfield_new();
    sim->magnetic[n] = cow_dfield_new();

    cow_dfield_setname(sim->velocity[n], "velocity");
    cow_dfield_setname(sim->magnetic[n], "magnetic");
    cow_dfield_setdomain(sim->velocity[n], sim->domain);
    cow_dfield_setdomain(sim->magnetic[n], sim->domain);

    cow_dfield_addmember(sim->velocity[n], "u1");
    cow_dfield_addmember(sim->velocity[n], "u2");
    cow_dfield_addmember(sim->velocity[n], "u3");
    cow_dfield_addmember(sim->magnetic[n], "b1");
    cow_dfield_addmember(sim->magnetic[n], "b2");
    cow_dfield_addmember(sim->magnetic[n], "b3");

    cow_dfield_commit(sim->velocity[n]);
    cow_dfield_commit(sim->magnetic[n]);
  }
}



/*
 * Sim destructor
 * =====================================================================
 */
void mhd_sim_free(struct mhd_sim *sim)
{
  for (int n=0; n<6; ++n) {
    cow_dfield_del(sim->velocity[n]);
    cow_dfield_del(sim->magnetic[n]);
  }
  cow_domain_del(sim->domain);
}



/*
 * Set up the problem
 * =====================================================================
 */
int mhd_problem_setup(struct mhd_sim *sim, const char *problem_name)
{
  if (problem_name == NULL) {
    printf("\nproblems are:\n");
    printf("1. kh\n");
    printf("2. abc\n");
    if (sim) {
      printf("\nuser options are:\n");
      mhd_user_report(&sim->user);
    }
    return 0;
  }
  else if (!strcmp(problem_name, "kh")) {
    sim->initial_data = mhd_initial_data_kh;
    return 0;
  }
  else if (!strcmp(problem_name, "abc")) {
    sim->initial_data = mhd_initial_data_abc;
    return 0;
  }
  else {
    return 1;
  }
}



/*
 * Evaluate initial data
 * =====================================================================
 */
void mhd_sim_initial_data(struct mhd_sim *sim)
{
  jsw_rand_t R;
  jsw_seed(&R, cow_domain_getcartrank(sim->domain));


  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);
  double *u = cow_dfield_getdatabuffer(sim->velocity[0]);
  double *b = cow_dfield_getdatabuffer(sim->magnetic[0]);

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    double x[4] = {0,
    		   cow_domain_positionatindex(sim->domain, 0, i),
    		   cow_domain_positionatindex(sim->domain, 1, j),
    		   cow_domain_positionatindex(sim->domain, 2, k)};

    sim->initial_data(sim, x, &u[m], &b[m]);
    
    /* We add a white-noise perturbation to the velocity if pert < 0.0 */
    if (sim->user.pert < 0.0) {
      double du1 = jsw_random_double(&R, -1.0, 1.0);
      double du2 = jsw_random_double(&R, -1.0, 1.0);
      double du3 = jsw_random_double(&R, -1.0, 1.0);
      u[m+1] += du1 * sim->user.pert;
      u[m+2] += du2 * sim->user.pert;
      u[m+3] += du3 * sim->user.pert;
    }
    else {

    }
  }

  cow_dfield_syncguard(sim->velocity[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
  cow_fft_helmholtzdecomp(sim->velocity[0], COW_PROJECT_OUT_DIV);
  cow_fft_helmholtzdecomp(sim->magnetic[0], COW_PROJECT_OUT_DIV);
}



/*
 * Apply Kreiss-Oliger operator to subtract high frequencies
 * =====================================================================
 */
void mhd_sim_kreiss_oliger(struct mhd_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);
  double *u = cow_dfield_getdatabuffer(sim->velocity[0]);
  double *b = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *du = cow_dfield_getdatabuffer(sim->velocity[1]);
  double *db = cow_dfield_getdatabuffer(sim->magnetic[1]);

  double KO_const = 0.0;

  switch (MHD_DISSIPATION_ORDER) {
  case 2: KO_const = -1./4 ; break;
  case 4: KO_const = -1./16; break;
  case 6: KO_const = -1./64; break;
  }

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    for (int d=1; d<=3; ++d) {
      du[m+d] = KOV(u,d);
      db[m+d] = KOV(b,d);
    }
  }

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    double eps = sim->user.eps;

    for (int d=1; d<=3; ++d) {
      u[m+d] -= eps * KO_const * du[m+d];
      b[m+d] -= eps * KO_const * db[m+d];
    }
  }

  cow_dfield_syncguard(sim->velocity[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
}



/*
 * Advance the simulation by one Runge-Kutta substep
 * =====================================================================
 */
void mhd_sim_advance_rk(struct mhd_sim *sim, int RKstep)
{
  double RKparam_array[5] = {0.0, 0.5, 0.5, 1.0};
  double RKparam = RKparam_array[RKstep];

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);

  double *u0 = cow_dfield_getdatabuffer(sim->velocity[0]);
  double *b0 = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *u  = cow_dfield_getdatabuffer(sim->velocity[1]);
  double *b  = cow_dfield_getdatabuffer(sim->magnetic[1]);
  double *dtu = cow_dfield_getdatabuffer(sim->velocity[RKstep+2]);
  double *dtb = cow_dfield_getdatabuffer(sim->magnetic[RKstep+2]);

  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt = sim->status.time_step;

  double KO_const = 0.0;

  switch (MHD_DISSIPATION_ORDER) {
  case 2: KO_const = -1./4 ; break;
  case 4: KO_const = -1./16; break;
  case 6: KO_const = -1./64; break;
  }


  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {
    int m = INDV(i,j,k);
    u[m+1] = u0[m+1] + dt * RKparam * dtu[m+1];
    u[m+2] = u0[m+2] + dt * RKparam * dtu[m+2];
    u[m+3] = u0[m+3] + dt * RKparam * dtu[m+3];

    b[m+1] = b0[m+1] + dt * RKparam * dtb[m+1];
    b[m+2] = b0[m+2] + dt * RKparam * dtb[m+2];
    b[m+3] = b0[m+3] + dt * RKparam * dtb[m+3];
  }
  cow_dfield_syncguard(sim->velocity[1]);
  cow_dfield_syncguard(sim->magnetic[1]);
  
  
  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) time-derivative register, reading from the
   * [1] field register
   * ===========================================================================
   */

  double nu = sim->user.nu;

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {
    int m = INDV(i,j,k);

    double Du1[4] = { 0, D1(u,1), D2(u,1), D3(u,1) };
    double Du2[4] = { 0, D1(u,2), D2(u,2), D3(u,2) };
    double Du3[4] = { 0, D1(u,3), D2(u,3), D3(u,3) };
    double d2u[4] = { 0, DEL2(u,1), DEL2(u,2), DEL2(u,3) };

    dtu[m+1] = -(u[m+1]*Du1[1] + u[m+2]*Du1[2] + u[m+3]*Du1[3]) + nu*d2u[1];
    dtu[m+2] = -(u[m+1]*Du2[1] + u[m+2]*Du2[2] + u[m+3]*Du2[3]) + nu*d2u[2];
    dtu[m+3] = -(u[m+1]*Du3[1] + u[m+2]*Du3[2] + u[m+3]*Du3[3]) + nu*d2u[3];

    dtb[m+1] = 0;
    dtb[m+2] = 0;
    dtb[m+3] = 0;
  }
  cow_fft_helmholtzdecomp(sim->velocity[RKstep+2], COW_PROJECT_OUT_DIV);
}



/*
 * Average Runge-Kutta substeps to complete a full time step
 * =====================================================================
 */
void mhd_sim_average_rk(struct mhd_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);

  int si = cow_dfield_getstride(sim->velocity[0], 0);
  int sj = cow_dfield_getstride(sim->velocity[0], 1);
  int sk = cow_dfield_getstride(sim->velocity[0], 2);

  double *u = cow_dfield_getdatabuffer(sim->velocity[0]);
  double *b = cow_dfield_getdatabuffer(sim->magnetic[0]);

  double *dtu0 = cow_dfield_getdatabuffer(sim->velocity[2]);
  double *dtu1 = cow_dfield_getdatabuffer(sim->velocity[3]);
  double *dtu2 = cow_dfield_getdatabuffer(sim->velocity[4]);
  double *dtu3 = cow_dfield_getdatabuffer(sim->velocity[5]);

  double *dtb0 = cow_dfield_getdatabuffer(sim->magnetic[2]);
  double *dtb1 = cow_dfield_getdatabuffer(sim->magnetic[3]);
  double *dtb2 = cow_dfield_getdatabuffer(sim->magnetic[4]);
  double *dtb3 = cow_dfield_getdatabuffer(sim->magnetic[5]);

  double dt = sim->status.time_step;


  /* ===========================================================================
   * Average the RK substeps, write result into register 0
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = INDV(i,j,k);

    u[m+1] += dt/6 * (dtu0[m+1] + 2*dtu1[m+1] + 2*dtu2[m+1] + dtu3[m+1]);
    u[m+2] += dt/6 * (dtu0[m+2] + 2*dtu1[m+2] + 2*dtu2[m+2] + dtu3[m+2]);
    u[m+3] += dt/6 * (dtu0[m+3] + 2*dtu1[m+3] + 2*dtu2[m+3] + dtu3[m+3]);

    b[m+1] += dt/6 * (dtb0[m+1] + 2*dtb1[m+1] + 2*dtb2[m+1] + dtb3[m+1]);
    b[m+2] += dt/6 * (dtb0[m+2] + 2*dtb1[m+2] + 2*dtb2[m+2] + dtb3[m+2]);
    b[m+3] += dt/6 * (dtb0[m+3] + 2*dtb1[m+3] + 2*dtb2[m+3] + dtb3[m+3]);
  }


  cow_dfield_syncguard(sim->velocity[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
}



/*
 * Advance the simulation by one full iteration
 * =====================================================================
 */
void mhd_sim_advance(struct mhd_sim *sim)
{
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt_light = MIN3(dx, dy, dz);

  /* !!! will need to revisit this to put in sound speed !!! */
  sim->status.time_step = sim->user.cfl * dt_light;

  mhd_sim_advance_rk(sim, 0);
  mhd_sim_advance_rk(sim, 1);
  mhd_sim_advance_rk(sim, 2);
  mhd_sim_advance_rk(sim, 3);
  mhd_sim_average_rk(sim);

  if (sim->user.eps > 0.0) {
    mhd_sim_kreiss_oliger(sim);
  }

  sim->status.iteration += 1;
  sim->status.time_simulation += sim->status.time_step;
}



void mhd_sim_write_checkpoint(struct mhd_sim *sim, const char *base_name)
{
  char chkpt_name[1024];
  if (base_name == NULL) {
    snprintf(chkpt_name, 1024, "%s/chkpt.%04d.h5",
	     sim->user.outdir,
	     sim->status.checkpoint_number);
  }
  else {
    snprintf(chkpt_name, 1024, "%s/chkpt.%s.h5",
	     sim->user.outdir,
	     base_name);
  }

  cow_dfield_write(sim->velocity[0], chkpt_name);
  cow_dfield_write(sim->magnetic[0], chkpt_name);

  if (cow_domain_getcartrank(sim->domain) == 0) {
    mhd_read_write_status(&sim->status, chkpt_name, 'w');
    mhd_read_write_user(&sim->user, chkpt_name, 'w');
  }
}



/*
 * Main function
 * =====================================================================
 */
int main(int argc, char **argv)
{
  cow_init(0, NULL, 0);


  
  int norun_main = 0;
  int restarted_run = 0;
  char logfile_name[1024];
  char anlfile_name[1024];
  char slcfile_name[1024];
  struct mhd_sim sim;

  mhd_user_set_defaults(&sim.user);

  
  /*
   * Print a help message
   * ===================================================================
   */
  printf("\nIncompressible MHD solver\n");
  printf("Jonathan Zrake, Stanford University (2016)\n");

  if (argc == 1) {
    printf("usage: mhd <problem-name> [opt1=val1] [opt2=val2]\n");
    mhd_problem_setup(&sim, NULL);
    cow_finalize();
    return 0;
  }
  else {
    strncpy(sim.problem_name, argv[1], 1024);
  }



  /*
   * Set up the problem defaults
   * ===================================================================
   */

  if (strstr(argv[1], ".h5") != 0) {
   
    norun_main += mhd_read_write_status(&sim.status, argv[1], 'r');
    norun_main += mhd_read_write_user(&sim.user, argv[1], 'r');

    if (norun_main == 0) {
      restarted_run = 1;
    }
  }
  else {
    mhd_status_set_defaults(&sim.status);
  }

  
  /*
   * Scan command line arguments (override restart user restart)
   * ===================================================================
   */
  for (int n=2; n<argc; ++n) {
    norun_main += mhd_user_set_from_arg(&sim.user, argv[n]);
  }

  
  mhd_user_report(&sim.user);


  if (norun_main) {
    cow_finalize();
    return 0;
  }
  else if (mhd_problem_setup(&sim, sim.problem_name)) {
    printf("[ffe] error: unkown problem name: '%s', choose one of\n",
	   sim.problem_name);
    mhd_problem_setup(NULL, NULL);
    cow_finalize();
    return 0;
  }


  mhd_sim_init(&sim);

  /* cow_domain_setcollective(sim.domain, sim.io_use_collective); */
  /* cow_domain_setchunk(sim.domain, sim.io_use_chunked); */
  /* cow_domain_setalign(sim.domain, */
  /* 		      sim.io_align_threshold * 1024, */
  /* 		      sim.io_disk_block_size * 1024); */


  if (restarted_run) {

    cow_dfield_read(sim.magnetic[0], argv[1]);
    cow_dfield_read(sim.velocity[0], argv[1]);

  }

  else {

    sim.status.time_last_checkpoint = -sim.user.cpi;
    sim.status.checkpoint_number = -1;
    mhd_sim_initial_data(&sim);

  }


  int local_grid_size = cow_domain_getnumlocalzonesinterior(sim.domain,
							    COW_ALL_DIMS);


  snprintf(logfile_name, 1024, "%s/mhd.dat"    , sim.user.outdir);
  snprintf(anlfile_name, 1024, "%s/analysis.h5", sim.user.outdir);
  snprintf(slcfile_name, 1024, "%s/slices.h5",   sim.user.outdir);


  if (cow_domain_getcartrank(sim.domain) == 0) {

    /*
     * Set up the problem directory and output log
     */

    FILE *logf = NULL;

    mkdir(sim.user.outdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (restarted_run) {
      mhd_truncate_logfile(sim.status.time_simulation, logfile_name);
      logf = fopen(logfile_name, "a");
    }
    else {
      logf = fopen(logfile_name, "w");
    }

    
    if (logf == NULL) {
      printf("[mhd] error: could not open log file '%s'\n", logfile_name);
      norun_main += 1;
    }
    else {
      fclose(logf);
    }
  }



  /* Propagate any errors to all procs */
  norun_main += cow_domain_intsum(sim.domain, norun_main);  

  if (norun_main) {
    sim.user.tmax = 0.0;
  }


  
  while (sim.status.time_simulation < sim.user.tmax) {


    /*
     * Write a checkpoint if it's time
     * =================================================================
     */
    if (sim.status.time_simulation - sim.status.time_last_checkpoint >=
	sim.user.cpi && sim.user.cpi > 0.0) {

      sim.status.time_last_checkpoint += sim.user.cpi;
      sim.status.checkpoint_number += 1;

      mhd_sim_write_checkpoint(&sim, NULL);
    }


    /*
     * Handle post-processing and reductions
     * =================================================================
     */
    int iter = sim.status.iteration;
    if (iter % sim.user.measure_cadence == 0) mhd_sim_measure(&sim, &sim.status);
    if (iter % sim.user.analyze_cadence == 0) mhd_sim_analyze(&sim, &sim.status,
							      anlfile_name);


    if (sim.user.slice_cadence != 0 &&
	iter % sim.user.slice_cadence == 0) {
      char group[256];
      snprintf(group, 256, "%06d", sim.status.iteration);
      cow_dfield_write_slice(sim.velocity[0], slcfile_name, group, 1);
      cow_dfield_write_slice(sim.velocity[0], slcfile_name, group, 2);
    }

    if (cow_domain_getcartrank(sim.domain) == 0) {
      FILE *logf = fopen(logfile_name, "a");
      mhd_status_fprintf(&sim.status, logf);
      fclose(logf);
    }



    /*
     * Evolve the system
     * =================================================================
     */
    void *start_cycle = cow_start_clock();

    mhd_sim_advance(&sim);

    double seconds = cow_stop_clock(start_cycle);

    sim.status.kzps = 1e-3 * local_grid_size / seconds;

    if (sim.status.iteration % 1 == 0) {
      printf("[ffe] n=%06d t=%6.4e dt=%6.4e %3.2f kzps\n",
	     sim.status.iteration,
	     sim.status.time_simulation,
	     sim.status.time_step,
	     sim.status.kzps);
      fflush(stdout);
    }
  }


  if (sim.user.cpi > 0.0) {
    if (norun_main == 0) {
      if (sim.user.tmax == 0.0) {
	sim.status.checkpoint_number = 0;
	mhd_sim_write_checkpoint(&sim, NULL);
      }
      else {
	mhd_sim_write_checkpoint(&sim, "final");
      }
    }
  }


  mhd_sim_free(&sim);
  cow_finalize();
  return 0;
}
