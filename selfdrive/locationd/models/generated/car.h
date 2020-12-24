/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4545313164873459153);
void inv_err_fun(double *nom_x, double *true_x, double *out_5710345084490788900);
void H_mod_fun(double *state, double *out_1761468809653102083);
void f_fun(double *state, double dt, double *out_7925444319858141107);
void F_fun(double *state, double dt, double *out_3327887193246759214);
void h_25(double *state, double *unused, double *out_8192360643023472400);
void H_25(double *state, double *unused, double *out_2394081245720893807);
void h_24(double *state, double *unused, double *out_716181223335845285);
void H_24(double *state, double *unused, double *out_2074104887891188763);
void h_30(double *state, double *unused, double *out_7590942485410108846);
void H_30(double *state, double *unused, double *out_7219817665228289473);
void h_26(double *state, double *unused, double *out_2723045651308121618);
void H_26(double *state, double *unused, double *out_8098674400660103080);
void h_27(double *state, double *unused, double *out_5252710163025233571);
void H_27(double *state, double *unused, double *out_8507399653064914785);
void h_29(double *state, double *unused, double *out_1335117763451061592);
void H_29(double *state, double *unused, double *out_7403127639352430190);
void h_28(double *state, double *unused, double *out_3742549332889463446);
void H_28(double *state, double *unused, double *out_9104659100379877385);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
