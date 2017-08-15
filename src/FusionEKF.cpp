#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  /* Initialize matrices. */

  // State vector with horizontal and vertical positions and velocities.
  x = VectorXd(4);
  
  // State transition matrix, maps old states to new.
  F = MatrixXd(4,4); 

  // Measurement matrices, map states to measurements.
  // The laser measurement matrix is constant but the radar
  // matrix is actually a jacobian that is re-evaluated after
  // each prediction step.
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  Hj_ = MatrixXd(3, 4);

  // State covariance matrix, models Gaussian uncertainity present in
  // state estimates.
  P = MatrixXd::Identity(4,4);

  // Process noise covariance matrix. Covariance of inherent noise pesent
  // in the process itself. For example small accelerations in constant
  // speed motion process.
  Q = MatrixXd(4,4);

  // Measurement noise covariance matrices.
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    cout << "Initializing the filter..." << endl;
    ekf_.x_ = x;
    ekf_.x_ << 1.0, 1.0 ,1.0, 1.0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      auto x_state = Tools::polar2cart(measurement_pack.raw_measurements_);
      ekf_.x_ << x_state;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      auto x_pos = measurement_pack.raw_measurements_;
      ekf_.x_[0] = x_pos[0];
      ekf_.x_[1] = x_pos[1];
    }

    ekf_.F_ = F;
    ekf_.P_ = P;
    ekf_.Q_ = Q;

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "Initialization done!" << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  auto delta_t  = measurement_pack.timestamp_ - previous_timestamp_;

  if( delta_t > 0.0){
    // Update state transition.
    F << 1, 0, delta_t, 0,
         0, 1, 0, delta_t,
         0, 0, 1, 0,
         0, 0, 0, 1;
    
    
    // Update process noise
    double noise_ax = 9.0;
    double noise_ay = 9.0;
    auto delta_t_2 = delta_t * delta_t;
    auto delta_t_3 = delta_t * delta_t_2;
    auto delta_t_4 = delta_t * delta_t_3;
    
    Q << (0.25 * delta_t_4 * noise_ax), 0, (0.50 * delta_t_3 * noise_ax), 0,
         0, (0.25 * delta_t_4, noise_ay), 0, (0.50* delta_t_3* noise_ay),
         (0.50 * delta_t_3 * noise_ax), 0, (delta_t_2 * noise_ax), 0,
         0, (0.50* delta_t_3* noise_ay), 0, (delta_t_2 * noise_ay);

    
    ekf_.Predict();
  }


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

  previous_timestamp_ = measurement_pack.timestamp_;
}
