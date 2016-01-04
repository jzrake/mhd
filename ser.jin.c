#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mhd.h"


int {{app_name}}_user_set_from_arg(struct {{app_name}}_user *user, char *arg)
{
  int error = 0;
  if (0) {
  }
  {% for m in user_struct %}
  {% if m.dtype == 'double' and not m.array %}
  else if (!strncmp(arg, "{{m.name}}=", {{m.name|length+1}})) {
    sscanf(arg, "{{m.name}}=%lf", &user->{{m.name}});
  }
  {% endif %}
  {% if m.dtype == 'int' and not m.array %}
  else if (!strncmp(arg, "{{m.name}}=", {{m.name|length+1}})) {
    sscanf(arg, "{{m.name}}=%d", &user->{{m.name}});
  }
  {% endif %}
  {% if m.dtype == 'char' and not m.array %}
  else if (!strncmp(arg, "{{m.name}}=", {{m.name|length+1}})) {
    sscanf(arg, "{{m.name}}=%c", &user->{{m.name}});
  }
  {% endif %}
  {% if m.dtype == 'double' and m.array %}
  else if (!strncmp(arg, "{{m.name}}=", {{m.name|length+1}})) {
    int num = sscanf(arg, "{{m.name}}={{(['%lf']*m.array)|join(',')}}",
                     &user->{{m.printf_arg|join(', &user->')}});
    if (num != {{m.array}}) {
      printf("[{{app_name}}] error: badly formatted option '{{m.name}}'\n");
      error += 1;
    }
  }
  {% endif %}
  {% if m.dtype == 'int' and m.array %}
  else if (!strncmp(arg, "{{m.name}}=", {{m.name|length+1}})) {
    int num = sscanf(arg, "{{m.name}}={{(['%d']*m.array)|join(',')}}",
                     &user->{{m.printf_arg|join(', &user->')}});
    if (num != {{m.array}}) {
      printf("[{{app_name}}] error: badly formatted option '{{m.name}}'\n");
      error += 1;
    }
  }
  {% endif %}
  {% if m.dtype == 'char' and m.array %}
  else if (!strncmp(arg, "{{m.name}}=", {{m.name|length+1}})) {
    sscanf(arg, "{{m.name}}=%{{m.array}}s", user->{{m.name}});
  }
  {% endif %}
  {% endfor %}
  else {
    printf("[{{app_name}}] error: no such option '%s'\n", arg);
    error = 1;
  }
  return error;
}




void {{app_name}}_user_report(struct {{app_name}}_user *user)
{
  printf("{{48 * '-'}}\n");
  {% for m in user_struct %}
  printf("{{m.name}} {{'.' * (24 - m.name|length)}} {{m.printf_fmt}}\n",
         user->{{m.printf_arg|join(', user->')}});
  {% endfor %}
  printf("{{48 * '-'}}\n");
}



void {{app_name}}_user_set_defaults(struct {{app_name}}_user *user)
{
  {% for m in user_struct %}
  {% if m.dtype == 'char' and not m.array %}
  user->{{m.name}} = '{{m.default_value}}';
  {% endif %}
  {% if m.dtype == 'char' and m.array %}
  strcpy(user->{{m.name}}, "{{m.default_value}}");
  {% endif %}
  {% if m.dtype in ['int', 'double'] and not m.array %}
  user->{{m.name}} = {{m.default_value}};
  {% endif %}
  {% if m.dtype in ['int', 'double'] and m.array %}
  {% for i in range(m.array) %}
  user->{{m.name}}[{{i}}] = {{m.default_value[i]}};
  {% endfor %}
  {% endif %}
  {% endfor %}
}



void {{app_name}}_status_set_defaults(struct {{app_name}}_status *status)
{
  {% for m in status_struct %}
  {% if m.dtype == 'char' and not m.array %}
  status->{{m.name}} = '{{m.default_value}}';
  {% endif %}
  {% if m.dtype == 'char' and m.array %}
  strcpy(status->{{m.name}}, "{{m.default_value}}");
  {% endif %}
  {% if m.dtype in ['int', 'double'] and not m.array %}
  status->{{m.name}} = {{m.default_value}};
  {% endif %}
  {% if m.dtype in ['int', 'double'] and m.array %}
  {% for i in range(m.array) %}
  status->{{m.name}}[{{i}}] = {{m.default_value[i]}};
  {% endfor %}
  {% endif %}
  {% endfor %}
}







{% if app_use_hdf5 %}

#include <hdf5.h>
static int read_old_h5type(hid_t src_dset, hid_t dst_type, void *dst_data,
			   const char *type_name);
static int read_write_kernel(hid_t h5f, hid_t h5s, hid_t h5t,
			     char mode, char *name, void *data);



int {{app_name}}_read_write_status(struct {{app_name}}_status *status,
				   const char *chkpt_name, char mode)
{
  int error = 1;

#define ADD_MEM(nm, tp) H5Tinsert(h5t, #nm, HOFFSET(struct {{app_name}}_status, nm), tp)

  if (!H5Fis_hdf5(chkpt_name)) {
    printf("[{{app_name}}] error: '%s' does not exist\n", chkpt_name);
    return 1;
  }

  hid_t h5f = H5Fopen(chkpt_name, H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t h5s = H5Screate(H5S_SCALAR);
  hid_t h5t = H5Tcreate(H5T_COMPOUND, sizeof(struct {{app_name}}_status));

  {% for m in status_struct %}
  {% if m.dtype == 'char' and not m.array %}
  ADD_MEM({{m.name}}, H5T_C_S1);
  {% endif %}
  {% if m.dtype == 'char' and m.array %}
  if (1) {
    hid_t h5t_string = H5Tcopy(H5T_C_S1); H5Tset_size(h5t_string, {{m.array}});
    ADD_MEM({{m.name}}, h5t_string);
    H5Tclose(h5t_string);
  }
  {% endif %}
  {% if m.dtype in ['int', 'double'] and not m.array %}
  ADD_MEM({{m.name}}, H5T_NATIVE_{{m.dtype|upper}});
  {% endif %}
  {% if m.dtype in ['int', 'double'] and m.array %}
  {% for i in range(m.array) %}
  ADD_MEM({{m.name}}[{{i}}], H5T_NATIVE_{{m.dtype|upper}}); 
  {% endfor %}
  {% endif %}
  {% endfor %}
  
  error = read_write_kernel(h5f, h5s, h5t, mode, "status", status);

  H5Tclose(h5t);
  H5Sclose(h5s);
  H5Fclose(h5f);

#undef ADD_MEM

  return error;
}


int {{app_name}}_read_write_user(struct {{app_name}}_user *user,
				 const char *chkpt_name, char mode)
{
  int error = 1;

#define ADD_MEM(nm, tp) H5Tinsert(h5t, #nm, HOFFSET(struct {{app_name}}_user, nm), tp)

  if (!H5Fis_hdf5(chkpt_name)) {
    printf("[{{app_name}}] error: '%s' does not exist\n", chkpt_name);
    return 1;
  }

  hid_t h5f = H5Fopen(chkpt_name, H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t h5s = H5Screate(H5S_SCALAR);
  hid_t h5t = H5Tcreate(H5T_COMPOUND, sizeof(struct {{app_name}}_user));

  {% for m in user_struct %}
  {% if m.dtype == 'char' and not m.array %}
  ADD_MEM({{m.name}}, H5T_C_S1);
  {% endif %}
  {% if m.dtype == 'char' and m.array %}
  if (1) {
    hid_t h5t_string = H5Tcopy(H5T_C_S1); H5Tset_size(h5t_string, {{m.array}});
    ADD_MEM({{m.name}}, h5t_string);
    H5Tclose(h5t_string);
  }
  {% endif %}
  {% if m.dtype in ['int', 'double'] and not m.array %}
  ADD_MEM({{m.name}}, H5T_NATIVE_{{m.dtype|upper}});
  {% endif %}
  {% if m.dtype in ['int', 'double'] and m.array %}
  {% for i in range(m.array) %}
  ADD_MEM({{m.name}}[{{i}}], H5T_NATIVE_{{m.dtype|upper}}); 
  {% endfor %}
  {% endif %}
  {% endfor %}
  
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
      

  for (int n=0; n<dst_nmembers; ++n) {

    char *member_name = H5Tget_member_name(dst_type, n);
    int found = 0;

    for (int m=0; m<src_nmembers; ++m) {

      if (!strcmp(member_name, H5Tget_member_name(src_type, m))) {

	size_t src_member_offset = H5Tget_member_offset(src_type, m);
	size_t dst_member_offset = H5Tget_member_offset(dst_type, n);
	hid_t src_member_type = H5Tget_member_type(src_type, m);
	hid_t dst_member_type = H5Tget_member_type(dst_type, n);

	if (H5Tequal(src_member_type, dst_member_type) <= 0) {

	  printf("[{{app_name}}] error: incompatible type for member %s:%s\n",
		 type_name, member_name);
	  error += 1;

	}
	else {

	  size_t member_size = H5Tget_size(dst_member_type);
	  memcpy(dst_data + dst_member_offset, src_data + src_member_offset,
		 member_size);

	}

	H5Tclose(src_member_type);
	H5Tclose(dst_member_type);

	found = 1;
	break;
      }
    }

    if (found == 0) {
      printf("[{{app_name}}] warning: %s:%s not found in source data\n",
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

    h5d = H5Dcreate(h5f, name, h5t, h5s, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    error = 0;
  }

  else if (mode == 'r') {

    h5d = H5Dopen(h5f, name, H5P_DEFAULT);

    hid_t src_type = H5Dget_type(h5d);

    if (H5Tequal(src_type, h5t) <= 0) {

      printf("[{{app_name}}] warning: checkpoint '%s' out of date\n", name);

      if (read_old_h5type(h5d, h5t, data, "status")) {
	printf("[{{app_name}}] error: checkpoint '%s' could not be read\n", name);
      }
      else {
	error = 0; /* success */
      }

    }
    else {
      H5Dread(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      error = 0; /* success */
    }

    H5Tclose(src_type);

  }

  H5Dclose(h5d);
  return error;
}

{% endif %}
