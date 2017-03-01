
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

#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

#include "function_object.h"
#include "g_n_solver.h"

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

void PointDataListPrint3(PointDataList const &a)
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


// Prints data

void PointDataListPrint4(FunctionObjectList const &r,
                         Eigen::VectorXd const &beta)
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

