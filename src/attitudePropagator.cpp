#include "attitudePropagator.h"
#include "rungeKutta.h"
#include "quatMath.h"
#include <math.h>



void attProp(Eigen::MatrixXd& stateHist, const Eigen::Vector4d& q0, const Eigen::Vector3d& w0, const Eigen::MatrixXd& control, const double& tf, const int& intervals, const AttitudeParams& p){
  Eigen::VectorXd state = Eigen::VectorXd::Zero(6);
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
      rungeKutta(state4,t,t+dt,dt,u,p,attDeriv,4);
      rungeKutta(state5,t,t+dt,dt,u,p,attDeriv,5);
      
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
        // TODO: THIS IS NOT OPTIMAL FOR QUATERNIONS
        Eigen::VectorXd stateFixed = ((t-tStar)*state5 + (tStar-(t-dt/sLast))*state)/(dt/sLast);
        stateFixed.head(4).normalize();
        stateHist.col(iter++) << stateFixed;
      }
      state = state5;

      t += dt;
      
      dt = std::min(dtConst,s*dt);
  }
  stateHist.col(intervals-1) << state;
}

Eigen::VectorXd attDeriv(const double& t, const Eigen::VectorXd& y, const Eigen::VectorXd& u, const Params& params){

    const AttitudeParams& p = dynamic_cast<const AttitudeParams&>(params);
    Eigen::Vector4d q = y.head(4);
    Eigen::Vector3d w = y.tail(3);

    Eigen::VectorXd yDot = Eigen::VectorXd::Zero(y.size());

    yDot.tail(3) = p.J.inverse()*u;
    yDot.head(4) = quatDot(q,w);

    return yDot;
}
