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
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_* H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  
  MatrixXd I = MatrixXd::Identity(P_.rows(),P_.cols());

  /* Update (correct) state means and variances. */
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  /* Compute intermediates. */
  VectorXd x_polar = Tools::cart2polar(x_);
  VectorXd y = z - x_polar;
  MatrixXd S = H_ * P_* H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  MatrixXd I = MatrixXd::Identity(P_.rows(),P_.cols());

  /* Update (correct) state means and variances. */
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
  
}
