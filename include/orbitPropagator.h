#include <Eigen/Dense>
//#include <ros/ros.h>

#ifndef ORBIT_PROPAGATOR_H
#define ORBIT_PROPAGATOR_H

void cwProp(Eigen::MatrixXd& stateHist, const Eigen::Vector3d& r0, const Eigen::Vector3d& v0, const Eigen::MatrixXd& control, const double& tf, const int& intervals, const OrbitalParams& p);

Eigen::VectorXd cwDeriv(const double& t, const Eigen::VectorXd& state, const Eigen::VectorXd& u, const OrbitalParams& p);

#endif
