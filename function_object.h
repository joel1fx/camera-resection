#ifndef CAMERA_RESECTION_FUNCTION_OBJECT_H_
#define CAMERA_RESECTION_FUNCTION_OBJECT_H_

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
