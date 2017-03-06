
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

#include "g_n_solver.h"

#include <iostream>
#include <Eigen/Dense>

// Return the length squared of vector a
double GNSolver::GetErr2(Eigen::VectorXd const &a) const
{
  double ret = 0.0;
  for(int i = 0;i < a.size();i++)
  {
    ret += a[i] * a[i];
  }

  return ret;
}

// Evaluate r function vector
Eigen::VectorXd GNSolver::EvalRFunctionVector(FunctionObjectList const &r,
                                              Eigen::VectorXd const &beta)
                                              const
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
double GNSolver::GetPartialD(int i, int indx, FunctionObjectList const &r,
                             Eigen::VectorXd const &beta) const
{
  double epsilon = 0.001;

  Eigen::VectorXd beta_epsilon = beta;
  beta_epsilon[indx] += epsilon;

  double ret = r[i](beta);
  double ret_epsilon = r[i](beta_epsilon);

  double partd = (ret - ret_epsilon) / epsilon;

  return partd;
}

// Calculate Jacobian
Eigen::MatrixXd GNSolver::Jacobian(FunctionObjectList const &r,
                                   Eigen::VectorXd const &beta) const
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
double GNSolver::CalcErr2(FunctionObjectList const &r,
                          Eigen::VectorXd const &beta) const
{
  Eigen::VectorXd vcol = EvalRFunctionVector(r, beta);

  return GetErr2(vcol);
}

//  The Gauss-Newton algorithm calculates successive approximations of
//  the beta vector by the following equation:
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

Eigen::VectorXd GNSolver::CalcDiff(FunctionObjectList const &r,
                                   Eigen::VectorXd const &beta) const
{
  Eigen::MatrixXd ejf = Jacobian(r, beta);

  Eigen::MatrixXd ejft = ejf.transpose();
  Eigen::MatrixXd ejf2 = ejft * ejf;
  Eigen::MatrixXd ejf3 = ejf2.inverse() * ejft;

  Eigen::MatrixXd vcol = EvalRFunctionVector(r, beta);

  return ejf3 * vcol;
}

//  The main solver, iterates to improve the beta vector, including
//  "dialing back" the next iteration of the beta vector until it's
//  less than the current beta vector.
Eigen::VectorXd GNSolver::operator()(FunctionObjectList const &r,
                                     Eigen::VectorXd &beta) const
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

