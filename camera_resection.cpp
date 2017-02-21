
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

/*  This is a Work in Progress. It is a single file which sets up and runs  */
/*  its own single test case.                                               */

/*  To build:  g++ -Wall -O2 camera_resection.cpp -o camera_resection       */

/*  This repo contains a function to perform camera resection. That is,     */
/*  determining the camera's position, orientation and pixel multiplier     */
/*  given the 2d screen coordinates and 3d cartesian coordinates of a       */
/*  number of points.                                                       */

/*  The algorithm uses a Gauss-Newton non-linear least squares solver.      */

/*  https://en.wikipedia.org/wiki/Gauss-Newton_algorithm                    */

/*  This implementation uses zero-based arrays unlike the Wikipedia form.   */

/*  Given there are m functions r[0],r[1],r[2],...r[m - 1], each a          */
/*  function of n variables beta[0],beta[1],beta[2],...beta[n - 1],         */

/*             m-1                                                          */
/*   minimize  Sum  ( r[i](beta[0],beta[1]...) )^2                          */
/*             i=0                                                          */

/*  In this case the r vector of functions is defined as follows:           */

/*  for r[i] where i is even, r[i] is the distance between the measured x   */
/*  screen coordinate of point number i/2  and the x coordinate of its      */
/*  back-projection                                                         */
/*  for r[i] where i is odd, r[i] is the distance between the measured y    */
/*  screen coordinate of point (i - 1)/2 and the y coordinate of its        */
/*  back-projection                                                         */

/*  The beta variables correspond to the camera parameters to be            */
/*  determined.                                                             */

/*  In this case, we adjust the camera parameters to minimize the sum of    */
/*  the square of the distances between each 2d point in the image and      */
/*  its calculated 2d position based on its 3d coordinates and the camera   */
/*  parameters.                                                             */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

#define DEG_TO_RAD (M_PI / 180.0)

#define NUM_POINTS 6
#define NUM_VARS 7


typedef struct point_data
{
  double xp;
  double yp;
  double x;
  double y;
  double z;
} point_data;

typedef std::vector<point_data> point_data_list;

// x_screen y_screen x_3d y_3d z_3d for each test point

point_data points_raw[] = {
  { 1464.0, 924.0, 0.017, 0.0, 0.017},
  { 1431.0, 2346.0, 0.0, 2.682, 0.0},
  { 431.0, 651.0, 0.247, 0.0, 2.089 },
  { 2850.0, 613.0, 2.704, 0.0, 0.017},
  { 2839.0, 1448.0, 2.800, 1.065, 0.310},
  { 2265.0, 174.0, 2.822, 0.0, 1.486 }
};

// Print x_screen y_screen x_3d y_3d z_3d for each point

void PointDataListPrint3(point_data_list& a)
{
  int n;

  std::cout << std::endl << "type: Point_Data_List2" << std::endl;
  n = a.size();
  std::cout << "num elements: " << n << std::endl;

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(6);
  for(int i = 0;i < n;i++)
  {
    std::cout << std::setw(12) << a[i].xp << " " << \
      std::setw(12) << a[i].yp << "   " << \
      std::setw(12) << a[i].x << " " << \
      std::setw(12) << a[i].y << " " << \
      std::setw(12) << a[i].z << std::endl;
  }

  return;
}


// Return the length squared of vector a

double GetErr2(Eigen::VectorXd& a)
{
  double ret = 0.0;
  for(int i = 0;i < a.size();i++)
  {
    ret += a[i] * a[i];
  }

  return ret;
}


// If i is even, calculate the x screen coordinate of point number i/2,
// given the 3d coordinates of the point and the camera parameters.
// If i is odd, calculate the y screen coordinate of point number (i-1)/2,
// given the 3d coordinates of the point and the camera parameters.
// This is called back projection.

double EvalBackProject(point_data_list& points, Eigen::VectorXd& beta, int i)
{
  double ret;

  int i2 = i / 2;

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
  Eigen::Affine3d t1(Eigen::Translation3d(-tx,-ty,-tz));

  Eigen::Matrix4d m4 = (rza * rxa * rya * t1).matrix();
  Eigen::Vector3d pn(points[i2].x, points[i2].y, points[i2].z);
  Eigen::Vector4d pnh = pn.homogeneous();
  Eigen::Vector4d pnh2 = m4 * pnh;
  Eigen::Vector3d pn2 = pnh2.hnormalized();

  /* handle case when point is behind camera */
  if(pn2[2] > 0.0)
  {
    return 1.e10;
  }

  if(i % 2 == 0)
  {
    ret = k * pn2[0] / -pn2[2];
  }
  else
  {
    ret = k * pn2[1] / -pn2[2];
  }

  return ret;
}

// Prints data

void PointDataListPrint4(point_data_list& a, Eigen::VectorXd& vbeta)
{
  std::cout << std::endl << "type: Point_Data_List2" << std::endl;
  int size = a.size();
  std::cout << "num elements: " << size << std::endl;

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(6);
  for(int i = 0;i < size;i++)
  {
    std::cout << std::setw(12) << a[i].xp << " " << \
      std::setw(12) << a[i].yp << "   (" << \
      std::setw(12) << EvalBackProject(a, vbeta, 2 * i) << \
      " " << std::setw(12) << \
      EvalBackProject(a, vbeta, 2 * i + 1) << ")    " << \
      std::setw(12) << a[i].x << " " << \
      std::setw(12) << a[i].y << " " << \
      std::setw(12) << a[i].z << std::endl;
  }

  return;
}

// Evaluate function r[i]

double EvalRFunction(int i, point_data_list& points, Eigen::VectorXd& vbeta)
{
  double ret;

  int i2 = i / 2;
  if(i % 2 == 0)
  {
    ret = points[i2].xp - EvalBackProject(points, vbeta, i);
  }
  else
  {
    ret = points[i2].yp - EvalBackProject(points, vbeta, i);
  }

  return ret;
}


// Evaluate r function vector

Eigen::VectorXd EvalRFunctionVector(point_data_list& points, Eigen::VectorXd& vbeta)
{
  int num_conds = 2 * points.size();
  Eigen::VectorXd ret(num_conds);

  for(int i = 0;i < num_conds;i++)
  {
    ret[i] = EvalRFunction(i, points, vbeta);
  }

  return ret;
}


// Get partial derivative

double GetPartialD(int i, int indx, point_data_list& points,
    Eigen::VectorXd& vbeta)
{
  double epsilon = 0.001;

  Eigen::VectorXd vbeta_epsilon = vbeta;
  vbeta_epsilon[indx] += epsilon;

  double ret = EvalBackProject(points, vbeta, i);
  double ret_epsilon = EvalBackProject(points, vbeta_epsilon, i);

  double partd = (ret_epsilon - ret) / epsilon;

  return partd;
}


// Calculate Jacobian

Eigen::MatrixXd Jacobian(point_data_list& points, Eigen::VectorXd& vbeta)
{
  int num_conds = 2 * points.size();
  Eigen::MatrixXd ret(num_conds, vbeta.size());

  for(int i = 0;i < num_conds;i++)
  {
    for(int j = 0;j < vbeta.size();j++)
    {
      ret(i, j) = GetPartialD(i, j, points, vbeta);
    }
  }

  return ret;
}


// Calculates the error of the function vector.

double CalcErr2(point_data_list& points, Eigen::VectorXd vbeta)
{
  double err2;
  Eigen::VectorXd r_vec;

  r_vec = EvalRFunctionVector(points, vbeta);
  err2 = GetErr2(r_vec);

  return err2;
}


//  The Gauss-Newton algorithm calculates successive approximations of the
//  beta vector by the following equation:
//
//  beta_next = beta_current - (J^T*J)^-1*J^T * r(beta_current)
//
//  where:
//  beta_next is the next iteration of the beta vector
//  beta_current is the current iteration of the beta vector
//  J is the Jacobian matrix
//  J^T is the transpose of the Jacobian matrix
//  r is the r vector of functions to be minimized
//
//  CalcDiff calculates and returns this part:
//
//  (J^T*J)^-1*J^T * r(beta_current)
//

Eigen::VectorXd CalcDiff(point_data_list& points2,
  Eigen::VectorXd& vbeta, double *err2)
{
  Eigen::VectorXd vcol;
  Eigen::MatrixXd ejf;
  Eigen::MatrixXd ejft;
  Eigen::MatrixXd ejf2;
  Eigen::MatrixXd ejf_inv;
  Eigen::MatrixXd ejf3;
  Eigen::VectorXd ediff;

  ejf = Jacobian(points2, vbeta);
  ejft = ejf.transpose();
  ejf2 = ejft * ejf;

  ejf_inv = ejf2.inverse();

  ejf3 = ejf_inv * ejft;

  vcol = EvalRFunctionVector(points2, vbeta);
  *err2 = GetErr2(vcol);

  ediff = ejf3 * vcol;

  return ediff;
}


//  Solve, the main solver, iterates to improve the beta vector, including
//  "dialling back" the next iteration of the beta vector until it's less
//  than the current beta vector.

Eigen::VectorXd Solve(point_data_list& points2, Eigen::VectorXd& vbeta)
{
  double atten;
  double err2;
  double new_err2;
  Eigen::VectorXd ediff;
  Eigen::VectorXd ediff1;
  Eigen::VectorXd vbeta_next;

  std::cout << "Here in Solve" << std::endl;

  /* iterate loop */
  for(int j = 0;j < 100;j++)
  {
    std::cout << "iteration " << j << std::endl;

    std::cout << "VBETA" << std::endl;
    std::cout << vbeta << std::endl;
    std::cout << std::endl;

    ediff = CalcDiff(points2, vbeta, &err2);

    std::cout << "err2 " << err2 << std::endl;

                ediff1 = ediff;

    atten = 1.0;
    for(int j1 = 0;j1 < 10;j1++)
    {
      for(int j2 = 0;j2 < ediff.rows();j2++)
      {
        ediff1(j2, 0) = atten * ediff(j2, 0);
      }

      vbeta_next = ediff1 + vbeta;

      new_err2 = CalcErr2(points2, vbeta_next);
      if(new_err2 < err2)
      {
        break;
      }
      else
      {
        atten *= 0.5;
      }
    }

    vbeta = vbeta_next;

    // stop iterating if the error is almost zero
    if(new_err2 < 0.0000001)
    {
      std::cout << "solved" << std::endl;

      break;
    }

    // stop iterating if the error has decreased minimally
    if(new_err2 / err2 > 0.99999)
    {
      std::cout << "converged" << std::endl;

      break;
    }
  }

  std::cout << "**********OUTPUT************" << std::endl;
  std::cout << "VBETA" << std::endl;
  std::cout << vbeta << std::endl;
  std::cout << std::endl;
  std::cout << "new_err2 " << new_err2 << std::endl;

        return vbeta;
}


// Main sets up the current hard-wired test

int main(int argc, char **argv)
{
  int i;
  point_data cur_point;
  point_data_list points2;

  double k = 1024.0;
  double tx = 10.0;
  double ty = 0.0;
  double tz = 10.0;
  double rx = 0.0;
  double ry = 45.0;
  double rz = 0.0;
  Eigen::VectorXd vbeta(7);
  Eigen::VectorXd vbeta_solved;

  if(argc < 1)
  {
    std::cerr << "usage: " << argv[0] << std::endl;

    return -1;
  }

  // Set up beta vector with initial estimate of camera parameters.
        vbeta << k, tx, ty, tz, rx, ry, rz;

  for(i = 0;i < NUM_POINTS;i++)
  {
    cur_point.x = points_raw[i].x;
    cur_point.y = points_raw[i].y;
    cur_point.z = points_raw[i].z;
    // Ad-hoc subtractions to make the optical center of the image (0,0)
    cur_point.xp = points_raw[i].xp - 1632.0;
    cur_point.yp = points_raw[i].yp - 1224.0;

    points2.push_back(cur_point);
  }

  PointDataListPrint3(points2);

  std::cout << std::endl;

  std::cout << "**********INPUT************" << std::endl;
  std::cout << "VBETA" << std::endl;
  std::cout << vbeta << std::endl;
  std::cout << std::endl;

  vbeta_solved = Solve(points2, vbeta);

  std::cout << "VBETA SOLVED" << std::endl;
  std::cout << vbeta_solved << std::endl;
  std::cout << std::endl;

  PointDataListPrint4(points2, vbeta_solved);

  return 0;
}

