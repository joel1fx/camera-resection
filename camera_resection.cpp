
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

#include "function_object.h"

#define DEG_TO_RAD (M_PI / 180.0)

#define NUM_POINTS 6


typedef struct PointData
{
  double xp;
  double yp;
  double x;
  double y;
  double z;
} PointData;

typedef std::vector<PointData> PointDataList;

// x_screen y_screen x_3d y_3d z_3d for each test point

PointData points_raw[] = {
  { 1464.0, 924.0, 0.017, 0.0, 0.017},
  { 1431.0, 2346.0, 0.0, 2.682, 0.0},
  { 431.0, 651.0, 0.247, 0.0, 2.089 },
  { 2850.0, 613.0, 2.704, 0.0, 0.017},
  { 2839.0, 1448.0, 2.800, 1.065, 0.310},
  { 2265.0, 174.0, 2.822, 0.0, 1.486 }
};

// Print x_screen y_screen x_3d y_3d z_3d for each point

void PointDataListPrint3(PointDataList& a)
{
  std::cout << std::endl << "type: Point_Data_List2" << std::endl;
  int size = a.size();
  std::cout << "num elements: " << size << std::endl;

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(6);
  for(int i = 0;i < size;i++)
  {
    std::cout << std::setw(12) << a[i].xp << " " << \
      std::setw(12) << a[i].yp << "   " << \
      std::setw(12) << a[i].x << " " << \
      std::setw(12) << a[i].y << " " << \
      std::setw(12) << a[i].z << std::endl;
  }

  return;
}

typedef std::vector<FunctionObject> FunctionObjectList;


// Prints data

void PointDataListPrint4(FunctionObjectList& r, Eigen::VectorXd& beta)
{
  std::cout << std::endl << "type: Point_Data_List2" << std::endl;
  int size = r.size() >> 1;
  std::cout << "num elements: " << size << std::endl;

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(6);
  for(int i = 0;i < size;i++)
  {
    int i2 = 2 * i;
    std::cout << std::setw(12) << r[i2].xp() << " " << \
      std::setw(12) << r[i2].yp() << "   (" << \
      std::setw(12) << r[i2].EvalBackProject(beta) << \
      " " << std::setw(12) << \
      r[i2 + 1].EvalBackProject(beta) << ")    " << \
      std::setw(12) << r[i2].x3d() << " " << \
      std::setw(12) << r[i2].y3d() << " " << \
      std::setw(12) << r[i2].z3d() << std::endl;
  }

  return;
}

class GNSolver
{
  public:
    GNSolver() {}
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

    // Evaluate r function vector
    Eigen::VectorXd EvalRFunctionVector(FunctionObjectList& r,
        Eigen::VectorXd& beta)
    {
      int num_conds = r.size();
      Eigen::VectorXd ret(num_conds);

      for(int i = 0;i < num_conds;i++)
      {
        ret[i] = r[i](beta);
      }

      return ret;
    }


    // Get partial derivative
    double GetPartialD(int i, int indx, FunctionObjectList& r,
        Eigen::VectorXd& beta)
    {
      double epsilon = 0.001;

      Eigen::VectorXd beta_epsilon = beta;
      beta_epsilon[indx] += epsilon;

      double ret = r[i].EvalBackProject(beta);
      double ret_epsilon = r[i].EvalBackProject(beta_epsilon);

      double partd = (ret_epsilon - ret) / epsilon;

      return partd;
    }

    // Calculate Jacobian
    Eigen::MatrixXd Jacobian(FunctionObjectList& r, Eigen::VectorXd& beta)
    {
      int num_conds = r.size();
      Eigen::MatrixXd ret(num_conds, beta.size());

      for(int i = 0;i < num_conds;i++)
      {
        for(int j = 0;j < beta.size();j++)
        {
          ret(i, j) = GetPartialD(i, j, r, beta);
        }
      }

      return ret;
    }

    // Calculates the error of the function vector.
    double CalcErr2(FunctionObjectList& r, Eigen::VectorXd beta)
    {
      Eigen::VectorXd vcol = EvalRFunctionVector(r, beta);

      return GetErr2(vcol);
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

    Eigen::VectorXd CalcDiff(FunctionObjectList& r, Eigen::VectorXd& beta)
    {
      Eigen::MatrixXd ejf = Jacobian(r, beta);

      Eigen::MatrixXd ejft = ejf.transpose();
      Eigen::MatrixXd ejf2 = ejft * ejf;
      Eigen::MatrixXd ejf3 = ejf2.inverse() * ejft;

      Eigen::MatrixXd vcol = EvalRFunctionVector(r, beta);

      return ejf3 * vcol;
    }

    //  The main solver, iterates to improve the beta vector, including
    //  "dialling back" the next iteration of the beta vector until it's less
    //  than the current beta vector.
    Eigen::VectorXd operator()(FunctionObjectList& r, Eigen::VectorXd& beta)
    {
      std::cout << "Entering Solve" << std::endl;

      double err2_min;
      // iterate loop
      for(int j = 0;j < 100;j++)
      {
        std::cout << "iteration " << j << std::endl;

        std::cout << "BETA" << std::endl;
        std::cout << beta << std::endl;
        std::cout << std::endl;

        double err2 = CalcErr2(r, beta);
        Eigen::VectorXd ediff = CalcDiff(r, beta);

        std::cout << "err2 (orig) " << err2 << std::endl;

        double atten = 1.0;
        err2_min = err2;
        Eigen::VectorXd beta_min = beta;
        for(int j1 = 0;j1 < 20;j1++)
        {
          Eigen::VectorXd beta_next = atten * ediff + beta;

          double err2_next = CalcErr2(r, beta_next);
          std::cout << "j1 " << j1 << " atten " << atten << " err2_next " <<
              err2_next << std::endl;
          if(err2_next < err2_min)
          {
            err2_min = err2_next;
            beta_min = beta_next;
            atten *= 0.7;
          }
          else
          {
            atten *= 0.7;
          }
        }

        beta = beta_min;

        // stop iterating if the error is almost zero
        if(err2_min < 0.0000001)
        {
          std::cout << "solved" << std::endl;

          break;
        }

        // stop iterating if the error has decreased minimally
        if(err2_min / err2 > 0.99999)
        {
          std::cout << "converged" << std::endl;

          break;
        }
      } // iterate

      std::cout << "**********OUTPUT************" << std::endl;
      std::cout << "BETA" << std::endl;
      std::cout << beta << std::endl;
      std::cout << std::endl;
      std::cout << "err2_min " << err2_min << std::endl;

      return beta;
    }
};


// Main sets up the current hard-wired test

int main(int argc, char **argv)
{
  PointDataList points;

  double k = 1024.0;
  double tx = 10.0;
  double ty = 0.0;
  double tz = 10.0;
  double rx = 0.0;
  double ry = 45.0;
  double rz = 0.0;
  Eigen::VectorXd beta(7);
  Eigen::VectorXd beta_solved;

  if(argc < 1)
  {
    std::cerr << "usage: " << argv[0] << std::endl;

    return -1;
  }

  // Set up beta vector with initial estimate of camera parameters.
        beta << k, tx, ty, tz, rx, ry, rz;

  FunctionObjectList r;
  for(int i = 0;i < NUM_POINTS;i++)
  {
    PointData cur_point;
    cur_point.x = points_raw[i].x;
    cur_point.y = points_raw[i].y;
    cur_point.z = points_raw[i].z;
    // Ad-hoc subtractions to make the optical center of the image (0,0)
    cur_point.xp = points_raw[i].xp - 1632.0;
    cur_point.yp = points_raw[i].yp - 1224.0;
    FunctionObject r0(cur_point.xp, cur_point.yp, cur_point.x,
        cur_point.y, cur_point.z, 2 * i);
    FunctionObject r1(cur_point.xp, cur_point.yp, cur_point.x,
        cur_point.y, cur_point.z, 2 * i + 1);

    points.push_back(cur_point);
    r.push_back(r0);
    r.push_back(r1);
  }

  PointDataListPrint3(points);

  std::cout << std::endl;

  std::cout << "**********INPUT************" << std::endl;
  std::cout << "BETA" << std::endl;
  std::cout << beta << std::endl;
  std::cout << std::endl;

  beta_solved = GNSolver()(r, beta);

  std::cout << "BETA SOLVED" << std::endl;
  std::cout << beta_solved << std::endl;
  std::cout << std::endl;

  printf("r %f\n", r[11](beta_solved));
  // printf("eval %f\n", EvalRFunction(11, points, beta_solved));

  PointDataListPrint4(r, beta_solved);

  return 0;
}

