#ifndef CAMERA_RESECTION_G_N_SOLVER_H_
#define CAMERA_RESECTION_G_N_SOLVER_H_

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

#include <Eigen/Dense>

#include "function_object.h"


// Gauss-Newton algorithm solver.
//
// https://en.wikipedia.org/wiki/Gauss-Newton_algorithm
//
// Example:
//    Eigen::VectorXd beta(7);
//    beta << k, tx, ty, tz, rx, ry, rz;
//    FunctionObjectList r;
//    intialize r...
//    Eigen::VectorXd beta_solved = GNSolver()(r, beta);

class GNSolver
{
  public:
    GNSolver() {}

    Eigen::VectorXd operator()(FunctionObjectList const &r,
                               Eigen::VectorXd &beta) const;
  private:
    double GetErr2(Eigen::VectorXd const &a) const;

    Eigen::VectorXd EvalRFunctionVector(FunctionObjectList const &r,
                                        Eigen::VectorXd const &beta) const;

    double GetPartialD(int i, int indx, FunctionObjectList const &r,
                       Eigen::VectorXd const &beta) const;

    Eigen::MatrixXd Jacobian(FunctionObjectList const &r,
                             Eigen::VectorXd const &beta) const;

    double CalcErr2(FunctionObjectList const &r,
                    Eigen::VectorXd const &beta) const;

    Eigen::VectorXd CalcDiff(FunctionObjectList const &r,
                             Eigen::VectorXd const &beta) const;
};

#endif // CAMERA_RESECTION_G_N_SOLVER_H_
