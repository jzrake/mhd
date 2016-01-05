

struct mhd_user {
    char problem_name[1024];
    char outdir[256];
    double tmax;
    double cpi;
    int N[3];
    int k2;
    double abc[3];
    double helicity;
    double pert;
    double cfl;
    double eps;
    double nu;
    int measure_cadence;
    int analyze_cadence;
    int slice_cadence;
    int num_pspec_bin;
    int max_pspec_bin;
};


struct mhd_status {
    int iteration;
    int checkpoint_number;
    double time_simulation;
    double time_step;
    double time_last_checkpoint;
    double kzps;
    double velocity_energy;
    double magnetic_energy;
    double velocity_monopole;
    double magnetic_monopole;
    double velocity_helicity;
    double magnetic_helicity;
    double velocity_L2;
    double magnetic_L2;
};


int mhd_user_set_from_arg(struct mhd_user *user, char *arg);
void mhd_user_report(struct mhd_user *user);
void mhd_user_set_defaults(struct mhd_user *user);
void mhd_status_set_defaults(struct mhd_status *status);


int mhd_read_write_status(struct mhd_status *status,
			  const char *chkpt_name, char mode);
int mhd_read_write_user(struct mhd_user *user,
			const char *chkpt_name, char mode);
