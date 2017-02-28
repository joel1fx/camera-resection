#ifndef CAMERA_RESECTION_FUNCTION_OBJECT_H_
#define CAMERA_RESECTION_FUNCTION_OBJECT_H_

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

class FunctionObject
{
  public:
    // Constructor
    FunctionObject(double xp, double yp, double x3d, double y3d, double z3d,
                   int indx) : _xp(xp), _yp(yp), _x3d(x3d), _y3d(y3d),
                   _z3d(z3d), _indx(indx) {}
    // Getter functions
    double xp();
    double yp();
    double x3d();
    double y3d();
    double z3d();
    // Evaluate function
    double operator()(Eigen::VectorXd& beta);
    // Evaluate back projection
    double EvalBackProject(Eigen::VectorXd& beta);
  private:
    double _xp, _yp;
    double _x3d, _y3d, _z3d;
    int _indx;
};

#endif // CAMERA_RESECTION_FUNCTION_OBJECT_H_
