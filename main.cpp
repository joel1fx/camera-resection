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

//  TODO: Move testing to test.cpp.

//  This repo implements a Gauss-Newton non-linear least squares solver
//  to do camera resection. That is, the solver determines the camera's
//  position, orientation and pixel scale given the 2d screen
//  coordinates and 3d cartesian coordinates of a number of points.
//
//  The algorithm uses a Gauss-Newton non-linear least squares solver.
//
//  https://en.wikipedia.org/wiki/Gauss-Newton_algorithm
//
//  This implementation uses zero-based arrays unlike the Wikipedia form.
//
//  Given there are m functions r[0],r[1],r[2],...r[m - 1], each a
//  function of n variables beta[0],beta[1],beta[2],...beta[n - 1],
//
//             m-1
//   minimize  Sum  ( r[i](beta[0],beta[1]...) )^2
//             i=0
//
//  To implement this two classes are defined. The first,
//  FunctionObject, represents each element of the r function vector.
//  The second, GNSolver (Gauss-Newton Solver) implements the actual
//  solver.
//
//  In this implementation the r vector of functions is defined as
//  follows:
//
//  for r[i] where i is even, r[i] is the distance between the measured x
//  screen coordinate of point number i/2  and the x coordinate of its
//  back-projection
//  for r[i] where i is odd, r[i] is the distance between the measured y
//  screen coordinate of point (i - 1)/2 and the y coordinate of its
//  back-projection
//
//  The beta variables correspond to the camera parameters to be
//  determined.
//
//  In this case, we adjust the camera parameters to minimize the sum of
//  the square of the distances between each 2d point in the image and
//  its calculated 2d position based on its 3d coordinates and the
//  camera parameters.
//
//  This main function tests the FunctionObject and GNSolver classes by
//  constructing objects to test a set of points with known screen
//  coordinates and 3d coordinates to see if the solver produces a
//  reasonable solution.
//

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


//
//  Function: PointDataListPrint3
//
//  Print point data.
//  screen x, y; 3d coords x, y, z.
//
//  Input parameters:
//
//  a -- list of points
//
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


//
//  Function: PointDataListPrint4
//
//  Print point data.
//  screen x, y; backprojected x, y; 3d coords x, y, z.
//
//  Input parameters:
//
//  r    -- vector containing the functions to evaluate
//  beta -- vector containing camera parameters to optimize, fed into
//          each r function
//

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


//
//  Function: main
//
//  Execute main program function.
//
//  main sets up the current hard-wired test
//
//  Input parameters:
//
//  argc -- argument count
//  argv -- argument vector
//

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

  if(argc < 2)
  {
    std::cerr << "usage: " << argv[0] << " pointDataFile" << std::endl;

    return -1;
  }

  FILE *file = fopen(argv[1], "r");
  if (file == NULL)
  {
    std::cerr << "can't open " << argv[1] << std::endl;

    return -1;
  }

  //  Initialize beta vector with initial estimate of camera parameters.
  beta << k, tx, ty, tz, rx, ry, rz;

  //  Initialize r function vector, constructing the function of each
  //  element.
  FunctionObjectList r;
  int i = 0;
  while(1)
  {
    PointData cur_point;
    int ret = fscanf(file, "%lf %lf %lf %lf %lf", &cur_point.xp,
        &cur_point.yp, &cur_point.x, &cur_point.y, &cur_point.z);
    if(ret != 5)
    {
      break;
    }
    cur_point.xp -= 1632.0;
    cur_point.yp -= 1224.0;
    FunctionObject r0(cur_point.xp, cur_point.yp, cur_point.x,
        cur_point.y, cur_point.z, 2 * i);
    FunctionObject r1(cur_point.xp, cur_point.yp, cur_point.x,
        cur_point.y, cur_point.z, 2 * i + 1);

    points.push_back(cur_point);
    r.push_back(r0);
    r.push_back(r1);

    i++;
  }
  fclose(file);
 
  PointDataListPrint3(points);

  std::cout << std::endl;

  // Print initial beta vector.
  std::cout << "**********INPUT************" << std::endl;
  std::cout << "BETA" << std::endl;
  std::cout << beta << std::endl;
  std::cout << std::endl;

  // Solve the problem.
  beta_solved = GNSolver()(r, beta);

  // Print solved beta vector.
  std::cout << "BETA SOLVED" << std::endl;
  std::cout << beta_solved << std::endl;
  std::cout << std::endl;

  printf("r %f\n", r[11](beta_solved));

  PointDataListPrint4(r, beta_solved);

  return 0;
}

