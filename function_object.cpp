
/***************************************************************************

The MIT License (MIT)

Copyright (c) 2013-2017 Joel E. Merritt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

***************************************************************************/

#include "function_object.h"

#include <vector>
#include <Eigen/Dense>

#define DEG_TO_RAD (M_PI / 180.0)

double FunctionObject::xp() const
{
  return _xp;
}

double FunctionObject::yp() const
{
  return _yp;
}

double FunctionObject::x3d() const
{
  return _x3d;
}

double FunctionObject::y3d() const
{
  return _y3d;
}

double FunctionObject::z3d() const
{
  return _z3d;
}

double FunctionObject::operator()(Eigen::VectorXd const &beta) const
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

double FunctionObject::EvalBackProject(Eigen::VectorXd const &beta) const
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

