#include <Eigen/Dense>
#include "orbitalParams.h"

#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

// Modified from http://mathfaculty.fullerton.edu/mathews/n2003/rungekuttafehlbergmod.html
// Note that y is modified to provide the answer rather than returning anything
void rungeKutta(Eigen::VectorXd& y, const double& t0, const double& tf, const double& dt, const Eigen::VectorXd& u, const OrbitalParams& p,
                Eigen::VectorXd (*dydt)(const double&, const Eigen::VectorXd&, const Eigen::VectorXd&, const OrbitalParams&), int order);

#endif