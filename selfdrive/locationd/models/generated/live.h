/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8591145062780468058);
void inv_err_fun(double *nom_x, double *true_x, double *out_5366121289682465607);
void H_mod_fun(double *state, double *out_3761673658166548496);
void f_fun(double *state, double dt, double *out_1635295267032129457);
void F_fun(double *state, double dt, double *out_592885697240666499);
void h_3(double *state, double *unused, double *out_8281074584219591627);
void H_3(double *state, double *unused, double *out_6651748594458866258);
void h_4(double *state, double *unused, double *out_7182696806422862661);
void H_4(double *state, double *unused, double *out_3714349636798514341);
void h_9(double *state, double *unused, double *out_3505899231045441646);
void H_9(double *state, double *unused, double *out_5654185195707198454);
void h_10(double *state, double *unused, double *out_5841938846673429364);
void H_10(double *state, double *unused, double *out_3098799257860960459);
void h_12(double *state, double *unused, double *out_3528178390144380618);
void H_12(double *state, double *unused, double *out_1607794125059304772);
void h_31(double *state, double *unused, double *out_3723675284790689541);
void H_31(double *state, double *unused, double *out_4282883900202724593);
void h_32(double *state, double *unused, double *out_2796902278867199432);
void H_32(double *state, double *unused, double *out_6233385621850476240);
void h_13(double *state, double *unused, double *out_466711105087916565);
void H_13(double *state, double *unused, double *out_3331843704135987412);
void h_14(double *state, double *unused, double *out_3505899231045441646);
void H_14(double *state, double *unused, double *out_5654185195707198454);
void h_19(double *state, double *unused, double *out_7228132108948637853);
void H_19(double *state, double *unused, double *out_8359879026899921981);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);