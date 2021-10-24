#include "tools.h"
#include <iostream>
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian.
   */

    MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  
  // pre-compute a set of terms to avoid repeated calculation and prevent division by 0
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = c1*c2;

    // check division by zero
  if (fabs(c1) < 0.00001) {
    c1 += 0.0001;
//     cout << "CalculateJacobian () - Error - Division by Zero" << endl;
//     return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;   
}

VectorXd Tools::Cartesian2Polar(const VectorXd& cartesian) {
  
  VectorXd polar(3);

  float px = cartesian[0];
  float py = cartesian[1];
  float vx = cartesian[2];
  float vy = cartesian[3];
  
  double displacement = sqrt(px*px + py*py);
  double heading_in_rad = atan2(py , px);
  // prevent division by 0
  if (displacement < 0.00001) {
  	displacement += 0.0001;
  }
  double rhodot = (px*vx + py*vy) / displacement;

   // (rads must be between +/- pi)
  if (heading_in_rad > M_PI) {
    heading_in_rad = heading_in_rad - 2 * M_PI;
  }
  if (heading_in_rad < -M_PI) {
    heading_in_rad = heading_in_rad + 2 * M_PI;
  }
	
  float heading_in_degrees = heading_in_rad; //* 180 / M_PI;
    
  polar << displacement, heading_in_degrees, rhodot;

  return polar;
  
}

VectorXd Tools::Polar2Cartesian(const VectorXd& polar) {
  
  VectorXd cartesian(4);
  
    // recover polar difference parameters
  float displacement = polar(0); // rho (distance)
  float heading_in_degrees = polar(1); // phi 
  float rhodot = polar(2); // rho dot (velocity)
  
  // convert degrees to radians 
  float rads = (heading_in_degrees * M_PI / 180); 

  // (rads must be between +/- pi)
  if (rads > M_PI) {
    rads = rads - 2 * M_PI;
  }
  if (rads < -M_PI) {
    rads = rads + 2 * M_PI;
  }

  float px = displacement * cos(rads);
  float py = displacement * sin(rads);
  // not enough info in rhodot to calculate velocity??
  float vx = rhodot * cos(rads);
  float vy = rhodot * cos(rads);

  //  assign calculated position state parameters and return
  cartesian << px, py, vx, vy;
  return cartesian;
  
}
