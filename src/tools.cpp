#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() > 0 && ground_truth.size()==estimations.size() ) {
    VectorXd sum(4);
    for (int i=0; i < estimations.size(); ++i) {
      VectorXd diff = estimations[i]-ground_truth[i];
      VectorXd sq = diff.array()*diff.array();
      sum+=sq;
    }

    sum = sum / estimations.size();
    rmse = sum.array().sqrt();
  } else {
    cout << "Invalid estimation or ground_truth data" << endl;
  }
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float pSum= px*px+py*py;
  if (pSum>0) {
  // compute the Jacobian matrix
    float pRoot = sqrt(pSum);
    float pv=py*vx-px*vy;
    Hj << px/pRoot, py/pRoot,0,0,
          -py/pSum,px/pSum,0,0,
          py*pv/pSum/pRoot,-px*pv/pSum/pRoot, px/pRoot, py/pRoot;

  } else {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
  }
  return Hj;
}
