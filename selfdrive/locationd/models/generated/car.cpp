
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4545313164873459153) {
   out_4545313164873459153[0] = delta_x[0] + nom_x[0];
   out_4545313164873459153[1] = delta_x[1] + nom_x[1];
   out_4545313164873459153[2] = delta_x[2] + nom_x[2];
   out_4545313164873459153[3] = delta_x[3] + nom_x[3];
   out_4545313164873459153[4] = delta_x[4] + nom_x[4];
   out_4545313164873459153[5] = delta_x[5] + nom_x[5];
   out_4545313164873459153[6] = delta_x[6] + nom_x[6];
   out_4545313164873459153[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5710345084490788900) {
   out_5710345084490788900[0] = -nom_x[0] + true_x[0];
   out_5710345084490788900[1] = -nom_x[1] + true_x[1];
   out_5710345084490788900[2] = -nom_x[2] + true_x[2];
   out_5710345084490788900[3] = -nom_x[3] + true_x[3];
   out_5710345084490788900[4] = -nom_x[4] + true_x[4];
   out_5710345084490788900[5] = -nom_x[5] + true_x[5];
   out_5710345084490788900[6] = -nom_x[6] + true_x[6];
   out_5710345084490788900[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_1761468809653102083) {
   out_1761468809653102083[0] = 1.0;
   out_1761468809653102083[1] = 0.0;
   out_1761468809653102083[2] = 0.0;
   out_1761468809653102083[3] = 0.0;
   out_1761468809653102083[4] = 0.0;
   out_1761468809653102083[5] = 0.0;
   out_1761468809653102083[6] = 0.0;
   out_1761468809653102083[7] = 0.0;
   out_1761468809653102083[8] = 0.0;
   out_1761468809653102083[9] = 1.0;
   out_1761468809653102083[10] = 0.0;
   out_1761468809653102083[11] = 0.0;
   out_1761468809653102083[12] = 0.0;
   out_1761468809653102083[13] = 0.0;
   out_1761468809653102083[14] = 0.0;
   out_1761468809653102083[15] = 0.0;
   out_1761468809653102083[16] = 0.0;
   out_1761468809653102083[17] = 0.0;
   out_1761468809653102083[18] = 1.0;
   out_1761468809653102083[19] = 0.0;
   out_1761468809653102083[20] = 0.0;
   out_1761468809653102083[21] = 0.0;
   out_1761468809653102083[22] = 0.0;
   out_1761468809653102083[23] = 0.0;
   out_1761468809653102083[24] = 0.0;
   out_1761468809653102083[25] = 0.0;
   out_1761468809653102083[26] = 0.0;
   out_1761468809653102083[27] = 1.0;
   out_1761468809653102083[28] = 0.0;
   out_1761468809653102083[29] = 0.0;
   out_1761468809653102083[30] = 0.0;
   out_1761468809653102083[31] = 0.0;
   out_1761468809653102083[32] = 0.0;
   out_1761468809653102083[33] = 0.0;
   out_1761468809653102083[34] = 0.0;
   out_1761468809653102083[35] = 0.0;
   out_1761468809653102083[36] = 1.0;
   out_1761468809653102083[37] = 0.0;
   out_1761468809653102083[38] = 0.0;
   out_1761468809653102083[39] = 0.0;
   out_1761468809653102083[40] = 0.0;
   out_1761468809653102083[41] = 0.0;
   out_1761468809653102083[42] = 0.0;
   out_1761468809653102083[43] = 0.0;
   out_1761468809653102083[44] = 0.0;
   out_1761468809653102083[45] = 1.0;
   out_1761468809653102083[46] = 0.0;
   out_1761468809653102083[47] = 0.0;
   out_1761468809653102083[48] = 0.0;
   out_1761468809653102083[49] = 0.0;
   out_1761468809653102083[50] = 0.0;
   out_1761468809653102083[51] = 0.0;
   out_1761468809653102083[52] = 0.0;
   out_1761468809653102083[53] = 0.0;
   out_1761468809653102083[54] = 1.0;
   out_1761468809653102083[55] = 0.0;
   out_1761468809653102083[56] = 0.0;
   out_1761468809653102083[57] = 0.0;
   out_1761468809653102083[58] = 0.0;
   out_1761468809653102083[59] = 0.0;
   out_1761468809653102083[60] = 0.0;
   out_1761468809653102083[61] = 0.0;
   out_1761468809653102083[62] = 0.0;
   out_1761468809653102083[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_7925444319858141107) {
   out_7925444319858141107[0] = state[0];
   out_7925444319858141107[1] = state[1];
   out_7925444319858141107[2] = state[2];
   out_7925444319858141107[3] = state[3];
   out_7925444319858141107[4] = state[4];
   out_7925444319858141107[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7925444319858141107[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7925444319858141107[7] = state[7];
}
void F_fun(double *state, double dt, double *out_3327887193246759214) {
   out_3327887193246759214[0] = 1;
   out_3327887193246759214[1] = 0;
   out_3327887193246759214[2] = 0;
   out_3327887193246759214[3] = 0;
   out_3327887193246759214[4] = 0;
   out_3327887193246759214[5] = 0;
   out_3327887193246759214[6] = 0;
   out_3327887193246759214[7] = 0;
   out_3327887193246759214[8] = 0;
   out_3327887193246759214[9] = 1;
   out_3327887193246759214[10] = 0;
   out_3327887193246759214[11] = 0;
   out_3327887193246759214[12] = 0;
   out_3327887193246759214[13] = 0;
   out_3327887193246759214[14] = 0;
   out_3327887193246759214[15] = 0;
   out_3327887193246759214[16] = 0;
   out_3327887193246759214[17] = 0;
   out_3327887193246759214[18] = 1;
   out_3327887193246759214[19] = 0;
   out_3327887193246759214[20] = 0;
   out_3327887193246759214[21] = 0;
   out_3327887193246759214[22] = 0;
   out_3327887193246759214[23] = 0;
   out_3327887193246759214[24] = 0;
   out_3327887193246759214[25] = 0;
   out_3327887193246759214[26] = 0;
   out_3327887193246759214[27] = 1;
   out_3327887193246759214[28] = 0;
   out_3327887193246759214[29] = 0;
   out_3327887193246759214[30] = 0;
   out_3327887193246759214[31] = 0;
   out_3327887193246759214[32] = 0;
   out_3327887193246759214[33] = 0;
   out_3327887193246759214[34] = 0;
   out_3327887193246759214[35] = 0;
   out_3327887193246759214[36] = 1;
   out_3327887193246759214[37] = 0;
   out_3327887193246759214[38] = 0;
   out_3327887193246759214[39] = 0;
   out_3327887193246759214[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3327887193246759214[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3327887193246759214[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3327887193246759214[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3327887193246759214[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3327887193246759214[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3327887193246759214[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3327887193246759214[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3327887193246759214[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3327887193246759214[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3327887193246759214[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3327887193246759214[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3327887193246759214[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3327887193246759214[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3327887193246759214[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3327887193246759214[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3327887193246759214[56] = 0;
   out_3327887193246759214[57] = 0;
   out_3327887193246759214[58] = 0;
   out_3327887193246759214[59] = 0;
   out_3327887193246759214[60] = 0;
   out_3327887193246759214[61] = 0;
   out_3327887193246759214[62] = 0;
   out_3327887193246759214[63] = 1;
}
void h_25(double *state, double *unused, double *out_8192360643023472400) {
   out_8192360643023472400[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2394081245720893807) {
   out_2394081245720893807[0] = 0;
   out_2394081245720893807[1] = 0;
   out_2394081245720893807[2] = 0;
   out_2394081245720893807[3] = 0;
   out_2394081245720893807[4] = 0;
   out_2394081245720893807[5] = 0;
   out_2394081245720893807[6] = 1;
   out_2394081245720893807[7] = 0;
}
void h_24(double *state, double *unused, double *out_716181223335845285) {
   out_716181223335845285[0] = state[4];
   out_716181223335845285[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2074104887891188763) {
   out_2074104887891188763[0] = 0;
   out_2074104887891188763[1] = 0;
   out_2074104887891188763[2] = 0;
   out_2074104887891188763[3] = 0;
   out_2074104887891188763[4] = 1;
   out_2074104887891188763[5] = 0;
   out_2074104887891188763[6] = 0;
   out_2074104887891188763[7] = 0;
   out_2074104887891188763[8] = 0;
   out_2074104887891188763[9] = 0;
   out_2074104887891188763[10] = 0;
   out_2074104887891188763[11] = 0;
   out_2074104887891188763[12] = 0;
   out_2074104887891188763[13] = 1;
   out_2074104887891188763[14] = 0;
   out_2074104887891188763[15] = 0;
}
void h_30(double *state, double *unused, double *out_7590942485410108846) {
   out_7590942485410108846[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7219817665228289473) {
   out_7219817665228289473[0] = 0;
   out_7219817665228289473[1] = 0;
   out_7219817665228289473[2] = 0;
   out_7219817665228289473[3] = 0;
   out_7219817665228289473[4] = 1;
   out_7219817665228289473[5] = 0;
   out_7219817665228289473[6] = 0;
   out_7219817665228289473[7] = 0;
}
void h_26(double *state, double *unused, double *out_2723045651308121618) {
   out_2723045651308121618[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8098674400660103080) {
   out_8098674400660103080[0] = 0;
   out_8098674400660103080[1] = 0;
   out_8098674400660103080[2] = 0;
   out_8098674400660103080[3] = 0;
   out_8098674400660103080[4] = 0;
   out_8098674400660103080[5] = 0;
   out_8098674400660103080[6] = 0;
   out_8098674400660103080[7] = 1;
}
void h_27(double *state, double *unused, double *out_5252710163025233571) {
   out_5252710163025233571[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8507399653064914785) {
   out_8507399653064914785[0] = 0;
   out_8507399653064914785[1] = 0;
   out_8507399653064914785[2] = 0;
   out_8507399653064914785[3] = 1;
   out_8507399653064914785[4] = 0;
   out_8507399653064914785[5] = 0;
   out_8507399653064914785[6] = 0;
   out_8507399653064914785[7] = 0;
}
void h_29(double *state, double *unused, double *out_1335117763451061592) {
   out_1335117763451061592[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7403127639352430190) {
   out_7403127639352430190[0] = 0;
   out_7403127639352430190[1] = 1;
   out_7403127639352430190[2] = 0;
   out_7403127639352430190[3] = 0;
   out_7403127639352430190[4] = 0;
   out_7403127639352430190[5] = 0;
   out_7403127639352430190[6] = 0;
   out_7403127639352430190[7] = 0;
}
void h_28(double *state, double *unused, double *out_3742549332889463446) {
   out_3742549332889463446[0] = state[5];
   out_3742549332889463446[1] = state[6];
}
void H_28(double *state, double *unused, double *out_9104659100379877385) {
   out_9104659100379877385[0] = 0;
   out_9104659100379877385[1] = 0;
   out_9104659100379877385[2] = 0;
   out_9104659100379877385[3] = 0;
   out_9104659100379877385[4] = 0;
   out_9104659100379877385[5] = 1;
   out_9104659100379877385[6] = 0;
   out_9104659100379877385[7] = 0;
   out_9104659100379877385[8] = 0;
   out_9104659100379877385[9] = 0;
   out_9104659100379877385[10] = 0;
   out_9104659100379877385[11] = 0;
   out_9104659100379877385[12] = 0;
   out_9104659100379877385[13] = 0;
   out_9104659100379877385[14] = 1;
   out_9104659100379877385[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
