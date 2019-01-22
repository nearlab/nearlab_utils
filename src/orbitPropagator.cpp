#include "orbitPropagator.h"
#include "rungeKutta.h"
#include <math.h>



void cwProp(Eigen::MatrixXd& stateHist, const Eigen::Vector3d& r0, const Eigen::Vector3d& v0, const Eigen::MatrixXd& control, const double& tf, const int& intervals, const OrbitalParams& p){
  Eigen::VectorXd state = Eigen::VectorXd::Zero(6);
  state.head(3) = r0.head(3)/p.nu;
  state.tail(3) = v0.head(3)*p.tau/p.nu;
  double dt = tf/(intervals)/100;
  double dtConst = tf/intervals;
  double t = 0;
  double err = 1e-7;
  int iter = 0;
  double s = 1;
  while(iter<intervals){
      // double tIterLow = std::floor(t/dtConst);
      // double tIterHigh = std::ceil(t/dtConst);
      // Eigen::VectorXd u;
      // if(tIterHigh<intervals){
      //   u = (t/dtConst-tIterLow)*control.col((int)tIterLow) + (tIterHigh-t/dtConst)*control.col((int)tIterHigh);
      // }else{
      //   u = control.col((int)tIterLow);
      // }
      Eigen::VectorXd u = control.col(iter);
      Eigen::VectorXd state4 = state;
      Eigen::VectorXd state5 = state;
      rungeKutta(state4,t,t+dt,dt,u,p,cwDeriv,4);
      rungeKutta(state5,t,t+dt,dt,u,p,cwDeriv,5);
      
      //ROS_INFO_STREAM("Time: " << t <<"\titer:"<<iter<< "\tdelta-t: "<< dt <<"\tDifference: " << (state5-state4).norm() <<"\nState4\n"<<state4<<"\nState5\n"<<state5);
      
      double sLast = s;  
      s = pow(err*dt/2/(state5-state4).norm(),.25); 
      if(std::isinf(s)){
        s = 1.2;//Heuristic
      } 
      
      if((state5-state4).cwiseAbs().minCoeff()>err){
        dt = s*dt;
        continue;
      }
      
      double tStar = (iter)*dtConst;
      if(t>tStar){
        Eigen::VectorXd stateFixed = ((t-tStar)*state5 + (tStar-(t-dt/sLast))*state)/(dt/sLast);
        stateHist.col(iter++) << stateFixed;
      }
      state = state5;

      t += dt;
      
      dt = std::min(dtConst,s*dt);
  }
  stateHist.col(intervals-1) << state;
  stateHist *= p.nu;
  stateHist.block(3,0,3,intervals) /= p.tau;
}

Eigen::VectorXd cwDeriv(const double& t, const Eigen::VectorXd& y, const Eigen::VectorXd& u, const Params& params){

    const OrbitalParams& p = dynamic_cast<const OrbitalParams&>(params);

    Eigen::VectorXd yDot = Eigen::VectorXd::Zero(y.size());
    
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6,6);
    A << 0, 0, 0, 1, 0, 0,
       0, 0, 0, 0, 1, 0,
       0, 0, 0, 0, 0, 1,
       3*p.n*p.n, 0, 0, 0, 2*p.n, 0,
       0, 0, 0, -2*p.n, 0, 0,
       0, 0, -p.n*p.n, 0, 0, 0;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6,6);
    B.block(3,3,3,3) = Eigen::MatrixXd::Identity(3,3)*( p.F*p.tau*p.tau/p.m/p.nu);

    Eigen::VectorXd v = u;
    v(3) = -v(3);
    yDot = A*y + B*v;
    return yDot;
}
