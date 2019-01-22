#ifndef QUATMATH_H
#define QUATMATH_H

#include <math.h>
#include <Eigen/Dense>

Eigen::Matrix3d quat2rot(const Eigen::Vector4d& q);
Eigen::Vector4d quatRot(const Eigen::Vector4d& q, const Eigen::Vector4d& dq);
Eigen::Matrix4d skew(const Eigen::Vector4d& q);
Eigen::Matrix3d crossProductEquivalent(const Eigen::Vector3d& a);
Eigen::Vector4d inverse(const Eigen::Vector4d& q);
Eigen::Vector4d quatDot(const Eigen::Vector4d& q, const Eigen::Vector3d& w);
#endif