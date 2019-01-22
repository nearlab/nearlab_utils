#include <Eigen/Dense>
#include "params.h"
//#include <ros/ros.h>

#ifndef ATTITUDE_PROPAGATOR_H
#define ATTITUDE_PROPAGATOR_H

void attProp(Eigen::MatrixXd& stateHist, const Eigen::Vector4d& q0, const Eigen::Vector3d& a0, const Eigen::MatrixXd& control, const double& tf, const int& intervals, const AttitudeParams& p);

Eigen::VectorXd attDeriv(const double& t, const Eigen::VectorXd& state, const Eigen::VectorXd& u, const Params& p);

#endif
