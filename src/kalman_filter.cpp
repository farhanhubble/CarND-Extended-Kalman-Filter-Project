#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /* Update state mean vector x_ and state covariance matrix P_,
     assuming process noise 'nu' to be zero mean, Gaussian. */
  x_ = F_ * x_ ;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /* Compute intermediates. */
  auto y = z - H_ * x_;
  auto S = H_ * P_* H_.transpose() + R_;
  auto K = P_ * H_.transpose() * S.inverse();
  
  auto I = MatrixXd::Identity(P_.rows(),P_.cols());

  /* Update (correct) state means and variances. */
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  std::function cart2polar = [const VectorXd x_](void) {
    px = x_[0];
    py = x_[1];
    vx = x_[2];
    vy = x_[3];

    polar = Eigen::VectorXd(3);

    auto rho =  cmath::sqrt(px*px+py*py);
    auto phi =  (px == 0) ? cmath::atan(cmath::INFINITY) : cmath::atan(py/px);
    auto rho_dot = (rho == 0) ? cmath::INFINITY :  ((px * vx + py * vy) / rho);

    polar << rho << phi << rho_dot;

    return polar;
  };


  std::function jacobian = [const VectorXd x_[](void){
    J = Eigen::MatrixXd(3,4);

    px = x_[0];
    py = x_[1];
    vx = x_[2];
    vy = x_[3];

    p_sqrd = (px*px + py*py);
    p = cmath::sqrt(p_sqrd);
    auto J11 = px / p;
    auto J12 = py / p;
    auto J13 = 0.0;
    auto J14 = 0.0;

    auto J21 = -py / p_sqrd;
    auto J22 =  px / p_sqrd;
    auto J13 = 0.0;
    auto J14 = 0.0;

    auto J31 = py*(vx*py-vy*px)/cmath::cube(p);
    auto J32 = px*(vy*px-vx*py)/cmath::cube(p);
    auto J33 = J11;
    auto J34 = J12;

    J << J11,J12,J13,J14,
        J21,J22,J23,J24,
        J31,J32,J33,J34;

    return J;
  };

  /* Compute intermediates. */
  // Replace H_ * x_ with cart2polar(x_), which is the non-linear measuremnet function.
  auto y = z - cart2polar(); 
  // Replace H_ with Jacobian (Hj) of the measurement function cart2polar() evaluated at previous.
  // values of x_;
  auto Hj = jacobian();
  auto S = H_ * P_* Hj.transpose() + R_;
  auto K = P_ * Hj.transpose() * S.inverse();

  auto I = MatrixXd::Identity(P_.rows(),P_.cols());

  /* Update (correct) state means and variances. */
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
  
}
