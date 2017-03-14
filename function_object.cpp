
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

//  This module defines the FunctionObject class whose objects represent
//  each function r in the Gauss-Newton non-linear least squares solver
//  implementation described in main.cpp.
//
//  Each function object is constructed to return either the screen x or
//  screen y error term for a given point. This error term is the
//  distance in pixels between the observed screen position of the point
//  and the computed screen position of the point, given the camera
//  described by the beta vector. These error terms are used as feedback
//  for the Gauss-Newton solver.

#include "function_object.h"

#include <vector>
#include <Eigen/Dense>

#define DEG_TO_RAD (M_PI / 180.0)


//
//  Member Function: FunctionObject::xp
//
//  Return the screen coordinate x (_xp) stored in the object by the
//  constructor.
//
//  No input parameters
//

double FunctionObject::xp() const
{
  return _xp;
}


//
//  Member Function: FunctionObject::yp
//
//  Return the screen coordinate y (_yp) stored in the object by the
//  constructor.
//
//  No input parameters
//

double FunctionObject::yp() const
{
  return _yp;
}


//
//  Member Function: FunctionObject::x3d
//
//  Return the 3d coordinate x (_x3d) stored in the object by the
//  constructor.
//
//  No input parameters
//

double FunctionObject::x3d() const
{
  return _x3d;
}


//
//  Member Function: FunctionObject::y3d
//
//  Return the 3d coordinate y (_y3d) stored in the object by the
//  constructor.
//
//  No input parameters
//

double FunctionObject::y3d() const
{
  return _y3d;
}


//
//  Member Function: FunctionObject::z3d
//
//  Return the 3d coordinate z (_z3d) stored in the object by the
//  constructor.
//
//  No input parameters
//

double FunctionObject::z3d() const
{
  return _z3d;
}


//
//  Member Function: FunctionObject::operator()
//
//  If an even valued index (_indx) has been stored in the object by the
//  constructor, return the difference between the point's x screen
//  coordinate and the x screen backprojection coordinate of the point's
//  3d coordinates backprojected by the camera specified by the input
//  beta parameter vector. Otherwise, return the difference between the
//  point's y screen coordinate and the y screen backprojection coordinate.
//
//  Input parameters:
//  beta -- vector containing camera parameters to optimize
//

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


//
// Member Function: FunctionObject::EvalBackProject
//

//  If an even valued index (_indx) has been stored in the object by the
//  constructor, return the x screen backprojection coordinate of the
//  point's 3d coordinates backprojected by the camera specified by the
//  input beta parameter vector. Otherwise, return the y screen
//  backprojection coordinate.
//
//  Input parameters:
//  beta -- vector containing camera parameters to optimize
//

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

  //  Handle case when point is behind camera.
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

