#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mhd.h"


int mhd_user_set_from_arg(struct mhd_user *user, char *arg)
{
    int error = 0;
    if (0) {
    } else if (!strncmp(arg, "outdir=", 7)) {
	sscanf(arg, "outdir=%256s", user->outdir);
    } else if (!strncmp(arg, "tmax=", 5)) {
	sscanf(arg, "tmax=%lf", &user->tmax);
    } else if (!strncmp(arg, "cpi=", 4)) {
	sscanf(arg, "cpi=%lf", &user->cpi);
    } else if (!strncmp(arg, "N=", 2)) {
	int num = sscanf(arg, "N=%d,%d,%d",
			 &user->N[0], &user->N[1], &user->N[2]);
	if (num != 3) {
	    printf("[mhd] error: badly formatted option 'N'\n");
	    error += 1;
	}
    } else if (!strncmp(arg, "k2=", 3)) {
	sscanf(arg, "k2=%d", &user->k2);
    } else if (!strncmp(arg, "abc=", 4)) {
	int num = sscanf(arg, "abc=%lf,%lf,%lf",
			 &user->abc[0], &user->abc[1], &user->abc[2]);
	if (num != 3) {
	    printf("[mhd] error: badly formatted option 'abc'\n");
	    error += 1;
	}
    } else if (!strncmp(arg, "helicity=", 9)) {
	sscanf(arg, "helicity=%lf", &user->helicity);
    } else if (!strncmp(arg, "pert=", 5)) {
	sscanf(arg, "pert=%lf", &user->pert);
    } else if (!strncmp(arg, "cfl=", 4)) {
	sscanf(arg, "cfl=%lf", &user->cfl);
    } else if (!strncmp(arg, "eps=", 4)) {
	sscanf(arg, "eps=%lf", &user->eps);
    } else if (!strncmp(arg, "nu=", 3)) {
	sscanf(arg, "nu=%lf", &user->nu);
    } else if (!strncmp(arg, "measure_cadence=", 16)) {
	sscanf(arg, "measure_cadence=%d", &user->measure_cadence);
    } else if (!strncmp(arg, "analyze_cadence=", 16)) {
	sscanf(arg, "analyze_cadence=%d", &user->analyze_cadence);
    } else if (!strncmp(arg, "num_pspec_bin=", 14)) {
	sscanf(arg, "num_pspec_bin=%d", &user->num_pspec_bin);
    } else if (!strncmp(arg, "max_pspec_bin=", 14)) {
	sscanf(arg, "max_pspec_bin=%d", &user->max_pspec_bin);
    } else {
	printf("[mhd] error: no such option '%s'\n", arg);
	error = 1;
    }
    return error;
}




void mhd_user_report(struct mhd_user *user)
{
    printf("------------------------------------------------\n");
    printf("outdir .................. %s\n", user->outdir);
    printf("tmax .................... %4.3lf\n", user->tmax);
    printf("cpi ..................... %4.3lf\n", user->cpi);
    printf("N ....................... %d %d %d \n",
	   user->N[0], user->N[1], user->N[2]);
    printf("k2 ...................... %d\n", user->k2);
    printf("abc ..................... %4.3lf %4.3lf %4.3lf \n",
	   user->abc[0], user->abc[1], user->abc[2]);
    printf("helicity ................ %4.3lf\n", user->helicity);
    printf("pert .................... %4.3lf\n", user->pert);
    printf("cfl ..................... %4.3lf\n", user->cfl);
    printf("eps ..................... %4.3lf\n", user->eps);
    printf("nu ...................... %4.3lf\n", user->nu);
    printf("measure_cadence ......... %d\n", user->measure_cadence);
    printf("analyze_cadence ......... %d\n", user->analyze_cadence);
    printf("num_pspec_bin ........... %d\n", user->num_pspec_bin);
    printf("max_pspec_bin ........... %d\n", user->max_pspec_bin);
    printf("------------------------------------------------\n");
}



void mhd_user_set_defaults(struct mhd_user *user)
{
    strcpy(user->outdir, ".");
    user->tmax = 1.0;
    user->cpi = 1.0;
    user->N[0] = 128;
    user->N[1] = 128;
    user->N[2] = 1;
    user->k2 = 1;
    user->abc[0] = 1;
    user->abc[1] = 1;
    user->abc[2] = 0;
    user->helicity = 1;
    user->pert = 0.0;
    user->cfl = 0.15;
    user->eps = 0.0;
    user->nu = 0.001;
    user->measure_cadence = 1;
    user->analyze_cadence = 128;
    user->num_pspec_bin = 256;
    user->max_pspec_bin = 1024;
}



void mhd_status_set_defaults(struct mhd_status *status)
{
    status->iteration = 0;
    status->checkpoint_number = 0;
    status->time_simulation = 0.0;
    status->time_step = 0.0;
    status->time_last_checkpoint = 0.0;
    status->kzps = 0.0;
    status->velocity_energy = 0.0;
    status->magnetic_energy = 0.0;
    status->velocity_monopole = 0.0;
    status->magnetic_monopole = 0.0;
    status->velocity_helicity = 0.0;
    status->magnetic_helicity = 0.0;
    status->velocity_L2 = 0.0;
    status->magnetic_L2 = 0.0;
}








#include <hdf5.h>
static int read_old_h5type(hid_t src_dset, hid_t dst_type, void *dst_data,
			   const char *type_name);
static int read_write_kernel(hid_t h5f, hid_t h5s, hid_t h5t,
			     char mode, char *name, void *data);



int mhd_read_write_status(struct mhd_status *status,
			  const char *chkpt_name, char mode)
{
    int error = 1;

#define ADD_MEM(nm, tp) H5Tinsert(h5t, #nm, HOFFSET(struct mhd_status, nm), tp)

    if (!H5Fis_hdf5(chkpt_name)) {
	printf("[mhd] error: '%s' does not exist\n", chkpt_name);
	return 1;
    }

    hid_t h5f = H5Fopen(chkpt_name, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t h5s = H5Screate(H5S_SCALAR);
    hid_t h5t = H5Tcreate(H5T_COMPOUND, sizeof(struct mhd_status));

    ADD_MEM(iteration, H5T_NATIVE_INT);
    ADD_MEM(checkpoint_number, H5T_NATIVE_INT);
    ADD_MEM(time_simulation, H5T_NATIVE_DOUBLE);
    ADD_MEM(time_step, H5T_NATIVE_DOUBLE);
    ADD_MEM(time_last_checkpoint, H5T_NATIVE_DOUBLE);
    ADD_MEM(kzps, H5T_NATIVE_DOUBLE);
    ADD_MEM(velocity_energy, H5T_NATIVE_DOUBLE);
    ADD_MEM(magnetic_energy, H5T_NATIVE_DOUBLE);
    ADD_MEM(velocity_monopole, H5T_NATIVE_DOUBLE);
    ADD_MEM(magnetic_monopole, H5T_NATIVE_DOUBLE);
    ADD_MEM(velocity_helicity, H5T_NATIVE_DOUBLE);
    ADD_MEM(magnetic_helicity, H5T_NATIVE_DOUBLE);
    ADD_MEM(velocity_L2, H5T_NATIVE_DOUBLE);
    ADD_MEM(magnetic_L2, H5T_NATIVE_DOUBLE);

    error = read_write_kernel(h5f, h5s, h5t, mode, "status", status);

    H5Tclose(h5t);
    H5Sclose(h5s);
    H5Fclose(h5f);

#undef ADD_MEM

    return error;
}


int mhd_read_write_user(struct mhd_user *user,
			const char *chkpt_name, char mode)
{
    int error = 1;

#define ADD_MEM(nm, tp) H5Tinsert(h5t, #nm, HOFFSET(struct mhd_user, nm), tp)

    if (!H5Fis_hdf5(chkpt_name)) {
	printf("[mhd] error: '%s' does not exist\n", chkpt_name);
	return 1;
    }

    hid_t h5f = H5Fopen(chkpt_name, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t h5s = H5Screate(H5S_SCALAR);
    hid_t h5t = H5Tcreate(H5T_COMPOUND, sizeof(struct mhd_user));

    if (1) {
	hid_t h5t_string = H5Tcopy(H5T_C_S1);
	H5Tset_size(h5t_string, 256);
	ADD_MEM(outdir, h5t_string);
	H5Tclose(h5t_string);
    }
    ADD_MEM(tmax, H5T_NATIVE_DOUBLE);
    ADD_MEM(cpi, H5T_NATIVE_DOUBLE);
    ADD_MEM(N[0], H5T_NATIVE_INT);
    ADD_MEM(N[1], H5T_NATIVE_INT);
    ADD_MEM(N[2], H5T_NATIVE_INT);
    ADD_MEM(k2, H5T_NATIVE_INT);
    ADD_MEM(abc[0], H5T_NATIVE_DOUBLE);
    ADD_MEM(abc[1], H5T_NATIVE_DOUBLE);
    ADD_MEM(abc[2], H5T_NATIVE_DOUBLE);
    ADD_MEM(helicity, H5T_NATIVE_DOUBLE);
    ADD_MEM(pert, H5T_NATIVE_DOUBLE);
    ADD_MEM(cfl, H5T_NATIVE_DOUBLE);
    ADD_MEM(eps, H5T_NATIVE_DOUBLE);
    ADD_MEM(nu, H5T_NATIVE_DOUBLE);
    ADD_MEM(measure_cadence, H5T_NATIVE_INT);
    ADD_MEM(analyze_cadence, H5T_NATIVE_INT);
    ADD_MEM(num_pspec_bin, H5T_NATIVE_INT);
    ADD_MEM(max_pspec_bin, H5T_NATIVE_INT);

    error = read_write_kernel(h5f, h5s, h5t, mode, "user", user);
    H5Tclose(h5t);
    H5Sclose(h5s);
    H5Fclose(h5f);

#undef ADD_MEM

    return error;
}



int read_old_h5type(hid_t src_dset, hid_t dst_type, void *dst_data,
		    const char *type_name)
{
    hid_t src_type = H5Dget_type(src_dset);
    int dst_nmembers = H5Tget_nmembers(dst_type);
    int src_nmembers = H5Tget_nmembers(src_type);
    size_t src_size = H5Tget_size(src_type);
    void *src_data = malloc(src_size);

    int error = 0;

    H5Dread(src_dset, src_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_data);


    for (int n = 0; n < dst_nmembers; ++n) {

	char *member_name = H5Tget_member_name(dst_type, n);
	int found = 0;

	for (int m = 0; m < src_nmembers; ++m) {

	    if (!strcmp(member_name, H5Tget_member_name(src_type, m))) {

		size_t src_member_offset =
		    H5Tget_member_offset(src_type, m);
		size_t dst_member_offset =
		    H5Tget_member_offset(dst_type, n);
		hid_t src_member_type = H5Tget_member_type(src_type, m);
		hid_t dst_member_type = H5Tget_member_type(dst_type, n);

		if (H5Tequal(src_member_type, dst_member_type) <= 0) {

		    printf
			("[mhd] error: incompatible type for member %s:%s\n",
			 type_name, member_name);
		    error += 1;

		} else {

		    size_t member_size = H5Tget_size(dst_member_type);
		    memcpy(dst_data + dst_member_offset,
			   src_data + src_member_offset, member_size);

		}

		H5Tclose(src_member_type);
		H5Tclose(dst_member_type);

		found = 1;
		break;
	    }
	}

	if (found == 0) {
	    printf("[mhd] warning: %s:%s not found in source data\n",
		   type_name, member_name);
	}

    }

    free(src_data);
    H5Tclose(src_type);

    return error;
}



int read_write_kernel(hid_t h5f, hid_t h5s, hid_t h5t,
		      char mode, char *name, void *data)
{
    hid_t h5d = -1;
    int error = 1;

    if (mode == 'w') {

	if (H5Lexists(h5f, name, H5P_DEFAULT)) {
	    H5Ldelete(h5f, name, H5P_DEFAULT);
	}

	h5d =
	    H5Dcreate(h5f, name, h5t, h5s, H5P_DEFAULT, H5P_DEFAULT,
		      H5P_DEFAULT);
	H5Dwrite(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	error = 0;
    }

    else if (mode == 'r') {

	h5d = H5Dopen(h5f, name, H5P_DEFAULT);

	hid_t src_type = H5Dget_type(h5d);

	if (H5Tequal(src_type, h5t) <= 0) {

	    printf("[mhd] warning: checkpoint '%s' out of date\n", name);

	    if (read_old_h5type(h5d, h5t, data, "status")) {
		printf("[mhd] error: checkpoint '%s' could not be read\n",
		       name);
	    } else {
		error = 0;	/* success */
	    }

	} else {
	    H5Dread(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	    error = 0;		/* success */
	}

	H5Tclose(src_type);

    }

    H5Dclose(h5d);
    return error;
}
