#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse<<0,0,0,0;

  /** the estimation vector size should not be zero
    * the estimation vector size should equal ground truth vector size
    */
  if(estimations.size()!=ground_truth.size()||estimations.size()==0)
  {
    std::cout<<"Error: Tools.cpp, CalculateRMSE, Invalid vector size.\n";
    return rmse;
  }
  
  // accumulate squared residuals
  for(unsigned int i=0;i<estimations.size();++i)
  {
    VectorXd residual = estimations[i]-ground_truth[i];
    
    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  
  //calculate the mean
  rmse = rmse / estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  return rmse;  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);

  if (x_state.size() != 4) {
    std::cout <<"Error: Tools.cpp, CalculateJacobian(), The state vector must have size 4\n";
    return Hj;
  }
  
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  double c1 = px * px + py * py;
  double c2 = sqrt(c1);
  double c3 = (c1*c2);
  
  //division by zero error
  if (fabs(c1) < 0.0001) {
    std::cout <<"Error: Tools.cpp, CalculateJacobian, divided by 0\n";
    return Hj;
  }

  /** compute the Jacobian matrix
  *update Radar measurement, To Calculate the P and X*/
  
  Hj<<px/c2,py/c2,0,0,
      -(py/c1),px/c1,0,0,
      py*(vx*py-vy*px)/c3,px*(vy*px-vx*py)/c3,px/c2,py/c2;
  
  return Hj;
  
}
