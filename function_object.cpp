
#include "function_object.h"

#include <Eigen/Dense>

#define DEG_TO_RAD (M_PI / 180.0)

double FunctionObject::xp()
{
  return _xp;
}

double FunctionObject::yp()
{
  return _yp;
}

double FunctionObject::x3d()
{
  return _x3d;
}

double FunctionObject::y3d()
{
  return _y3d;
}

double FunctionObject::z3d()
{
  return _z3d;
}

double FunctionObject::operator()(Eigen::VectorXd& beta)
{
  double ret;

  if(_indx % 2 == 0)
  {
    ret = _xp - EvalBackProject(beta);
  }
  else
  {
    ret = _yp - EvalBackProject(beta);
  }

  return ret;
}

// If _indx is even, calculate the x screen coordinate of point number
// _indx/2, given the 3d coordinates (_x3d, _y3d, _z3d) of the point
// and the camera parameters (beta).
// If _indx is odd, calculate the y screen coordinate of point number
// (_indx - 1)/2, given the 3d coordinates (_x3d, _y3d, _z3d) of the
// point and the camera parameters (beta).
// This is called back projection.

double FunctionObject::EvalBackProject(Eigen::VectorXd& beta)
{
  double k = beta[0];
  double tx = beta[1];
  double ty = beta[2];
  double tz = beta[3];
  double rx = beta[4] * DEG_TO_RAD; 
  double ry = beta[5] * DEG_TO_RAD;
  double rz = beta[6] * DEG_TO_RAD;

  Eigen::Affine3d rya = Eigen::Affine3d(Eigen::AngleAxisd(-ry,
      Eigen::Vector3d::UnitY()));
  Eigen::Affine3d rxa = Eigen::Affine3d(Eigen::AngleAxisd(-rx,
  Eigen::Vector3d::UnitX()));
  Eigen::Affine3d rza = Eigen::Affine3d(Eigen::AngleAxisd(-rz,
      Eigen::Vector3d::UnitZ()));
  Eigen::Affine3d t1(Eigen::Translation3d(-tx, -ty, -tz));

  Eigen::Matrix4d m4 = (rza * rxa * rya * t1).matrix();
  Eigen::Vector3d pn(_x3d, _y3d, _z3d);
  Eigen::Vector4d pnh = pn.homogeneous();
  Eigen::Vector4d pnh2 = m4 * pnh;
  Eigen::Vector3d pn2 = pnh2.hnormalized();

  // handle case when point is behind camera
  if(pn2[2] > 0.0)
  {
    return 1.e10;
  }

  double ret;
  if(_indx % 2 == 0)
  {
    ret = k * pn2[0] / -pn2[2];
  }
  else
  {
    ret = k * pn2[1] / -pn2[2];
  }

  return ret;
}

