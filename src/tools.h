#ifndef TOOLS_H_
#define TOOLS_H_
#include <cmath>
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * A helper method to transform measurements from cartesian to polar system.
  */
  static VectorXd cart2polar(const VectorXd& x);

  /**
  * A helper method to transform measurements from polar to cartesian system.
  */

  static VectorXd  polar2cart(const VectorXd& x);

};

#endif /* TOOLS_H_ */
