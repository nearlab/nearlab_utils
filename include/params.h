#include <Eigen/Dense>

#ifndef PARAMS_H
#define PARAMS_H

class Params{
  Params();
};

class OrbitalParams : public Params{
public:
  double m;//mass, kg
  double F;//Thrust, kg*m/s^2 (N)
  double tau;//time constant, s
  double nu;//distance constant, m
  double mu;//Gravitational Parameter for Earth kg*m^2/s^3
  Eigen::Vector3d r;//Orbital radius
  double n;//Mean rate
  OrbitalParams():m(100),F(.1),tau(100),nu(300),mu(3.986e14){
    this->r = Eigen::VectorXd::Zero(3);
    this->r << 42164000,0,0;
    this->n = sqrt(3.986e14/pow(42164000,3));
  }
  OrbitalParams(const double& m, const double& F, const double& tau, const double& nu, const Eigen::Vector3d& r, const double& mu=3.986e14):m(m),F(F),tau(tau),nu(nu),mu(mu),r(r){
    this->n = sqrt(mu/pow(r.norm(),3));
  }
};

class AttitudeParams : public Params{
public:
  // This will eventually need to include other things, presumably, for various attitude effects in space
  Eigen::Matrix3d J;//Orbital radius
  AttitudeParams(){
    this->J = Eigen::MatrixXd::Identity(3,3);
  }
  AttitudeParams(const Eigen::Matrix3d& J):J(J){}
};


#endif