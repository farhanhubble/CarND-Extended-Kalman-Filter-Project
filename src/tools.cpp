#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  int len = estimations.size();
  VectorXd rmse  = VectorXd::Zero(estimations[0].size());

  for(int i = 0; i< len; i++){
    VectorXd error  = estimations[i] - ground_truth[i];
    error = error.array() * error.array();
    rmse  += error;
  }

  rmse = rmse / len;
  rmse = rmse.array().sqrt();

  return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd J = Eigen::MatrixXd(3,4);

  auto px = x_state[0];
  auto py = x_state[1];
  auto vx = x_state[2];
  auto vy = x_state[3];

  auto p_sqrd = (px*px + py*py);
  auto p = sqrt(p_sqrd);
  auto J11 = px / p;
  auto J12 = py / p;
  auto J13 = 0.0;
  auto J14 = 0.0;

  auto J21 = -py / p_sqrd;
  auto J22 =  px / p_sqrd;
  auto J23 = 0.0;
  auto J24 = 0.0;

  auto J31 = py*(vx*py-vy*px)/pow(p,3);
  auto J32 = px*(vy*px-vx*py)/pow(p,3);
  auto J33 = J11;
  auto J34 = J12;

  J << J11,J12,J13,J14,
      J21,J22,J23,J24,
      J31,J32,J33,J34;

  return J;
}



VectorXd Tools::cart2polar(const VectorXd& x) {
auto px = x[0];
auto py = x[1];
auto vx = x[2];
auto vy = x[3];

auto polar = Eigen::VectorXd(3);

auto rho =  sqrt(px*px+py*py);
if(rho == 0){
  rho = 0.000001;
}
auto phi =  (px == 0) ? atan(INFINITY) : atan2(py,px);
auto rho_dot = (rho == 0) ? INFINITY :  ((px * vx + py * vy) / rho);

polar << rho, phi, rho_dot;

return polar;
}


VectorXd Tools::polar2cart(const VectorXd& x) {
auto rho = x[0];
auto phi = x[1];
auto rho_dot = x[2];

auto cartesian = Eigen::VectorXd(4);

auto px =  rho * cos(phi);
auto py =  rho * sin(phi);
auto vx =  rho * cos(phi);
auto vy =  rho * sin(phi);

cartesian << px , py, vx, vy;

return cartesian;
}


void Tools::printMat(MatrixXd &m){
  cout << "\n" << m << "\n";
}

void Tools::printVec(VectorXd &v){
  cout << "\n" << v << "\n";
}
