

struct {{app_name}}_user
{
  {% for n in user_struct %}
  {{n.dtype}} {{n.name}}{{n.array and [n.array] or ''}};
  {% endfor %}
} ;


struct {{app_name}}_status
{
  {% for n in status_struct %}
  {{n.dtype}} {{n.name}}{{n.array and [n.array] or ''}};
  {% endfor %}
} ;


int {{app_name}}_user_set_from_arg(struct {{app_name}}_user *user, char *arg);
void {{app_name}}_user_report(struct {{app_name}}_user *user);
void {{app_name}}_user_set_defaults(struct {{app_name}}_user *user);
void {{app_name}}_status_set_defaults(struct {{app_name}}_status *status);


int {{app_name}}_read_write_status(struct {{app_name}}_status *status,
				   const char *chkpt_name, char mode);
int {{app_name}}_read_write_user(struct {{app_name}}_user *user,
				 const char *chkpt_name, char mode);


