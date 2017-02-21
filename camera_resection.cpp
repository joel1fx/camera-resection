
/* To build:  g++ -Wall -O2 camera_resection.cpp -o camera_resection */

/* based on http://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm */

/* This implementation uses zero-based arrays unlike the Wikipedia form. */

/* Given there are m functions r[0],r[1],r[2],...r[m - 1], each a function */
/* of n variables beta[0],beta[1],beta[2],...beta[n - 1], */
/* minimize Sum(i = 0 to m - 1) r[i](beta[0],beta[1]...)^2 */

/* In this case the r vector is defined as follows: */
/* for r[i] where i is even, r[i] is the distance between the measured x */
/* screen coordinate and the x coordinate of its back-projection */
/* for r[i] where i is odd, r[i] is the distance between the measured y */
/* screen coordinate and the y coordinate of its back-projection */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

#define DEG_TO_RAD (M_PI / 180.0)

#define NUM_POINTS 6
#define NUM_VARS 7

#define NUM_CONDS (2 * NUM_POINTS)

typedef struct point_data
{
	double xp;
	double yp;
	double x;
	double y;
	double z;
} point_data;

typedef std::vector<point_data> point_data_list2;

point_data points_raw[] = {
	{ 1464.0, 924.0, 0.017, 0.0, 0.017},
	{ 1431.0, 2346.0, 0.0, 2.682, 0.0},
	{ 431.0, 651.0, 0.247, 0.0, 2.089 },
	{ 2850.0, 613.0, 2.704, 0.0, 0.017},
	{ 2839.0, 1448.0, 2.800, 1.065, 0.310},
	{ 2265.0, 174.0, 2.822, 0.0, 1.486 }
};


void PointDataListPrint3(point_data_list2& a)
{
	int i;
	int n;

	std::cout << std::endl << "type: Point_Data_List2" << std::endl;
	n = a.size();
	std::cout << "num elements: " << n << std::endl;

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(6);
	for(i = 0;i < n;i++)
	{
		std::cout << std::setw(12) << a[i].xp << " " << \
			std::setw(12) << a[i].yp << "   " << \
			std::setw(12) << a[i].x << " " << \
			std::setw(12) << a[i].y << " " << \
			std::setw(12) << a[i].z << std::endl;
	}

	return;
}


double get_err2(Eigen::VectorXd& a)
{
	int i;
	double ret;

	ret = 0.0;
	for(i = 0;i < a.size();i++)
	{
		ret += a[i] * a[i];
	}

	return ret;
}


double get_function3(point_data_list2& points, Eigen::VectorXd& beta, int i)
{
	double ret;
	double k;
	double tx;
	double ty;
	double tz;
	double rx;
	double ry;
	double rz;
	int i2;

	i2 = i / 2;

	k = beta[0];
	tx = beta[1];
	ty = beta[2];
	tz = beta[3];
	rx = beta[4] * DEG_TO_RAD; 
	ry = beta[5] * DEG_TO_RAD;
	rz = beta[6] * DEG_TO_RAD;

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


double get_r2(point_data_list2& points, int i, Eigen::VectorXd& vbeta)
{
	int i2;
	double ret;

	i2 = i / 2;
	if(i % 2 == 0)
	{
		ret = points[i2].xp - get_function3(points, vbeta, i);
	}
	else
	{
		ret = points[i2].yp - get_function3(points, vbeta, i);
	}

	return ret;
}


Eigen::VectorXd get_r_col2(point_data_list2& points2, Eigen::VectorXd& vbeta)
{
	int i;
	Eigen::VectorXd ret(NUM_CONDS);

	for(i = 0;i < NUM_CONDS;i++)
	{
		ret[i] = get_r2(points2, i, vbeta);
	}

	return ret;
}


void PointDataListPrint4(point_data_list2& a, Eigen::VectorXd& vbeta)
{
	int i;
	int n;

	std::cout << std::endl << "type: Point_Data_List2" << std::endl;
	n = a.size();
	std::cout << "num elements: " << n << std::endl;

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(6);
	for(i = 0;i < n;i++)
	{
		std::cout << std::setw(12) << a[i].xp << " " << \
			std::setw(12) << a[i].yp << "   (" << \
			std::setw(12) << get_function3(a, vbeta, 2 * i) << \
			" " << std::setw(12) << \
			get_function3(a, vbeta, 2 * i + 1) << ")    " << \
			std::setw(12) << a[i].x << " " << \
			std::setw(12) << a[i].y << " " << \
			std::setw(12) << a[i].z << std::endl;
	}

	return;
}


/* get partial derivative */
double get_partiald(point_data_list2& points, Eigen::VectorXd& vbeta, int i,
	int indx)
{
	double ret0;
	double ret_plus;
	double epsilon = 0.001;
	double partd;
	Eigen::VectorXd vbeta1;
	
	vbeta1 = vbeta;

	vbeta1[indx] += epsilon;

	ret0 = get_function3(points, vbeta, i);
	ret_plus = get_function3(points, vbeta1, i);

	partd = (ret_plus - ret0) / epsilon;

	return partd;
}


Eigen::MatrixXd jacobian(point_data_list2& points2, Eigen::VectorXd& vbeta)
{
	int i, j;
	Eigen::MatrixXd ret(NUM_CONDS, vbeta.size());

	for(i = 0;i < NUM_CONDS;i++)
	{
		for(j = 0;j < vbeta.size();j++)
		{
			ret(i, j) = get_partiald(points2, vbeta, i, j);
		}
	}

	return ret;
}


double calc_err2(point_data_list2& points, Eigen::VectorXd vbeta)
{
	double err2;
	Eigen::VectorXd r_vec;

	r_vec = get_r_col2(points, vbeta);
	err2 = get_err2(r_vec);

	return err2;
}


Eigen::VectorXd calc_diff(point_data_list2& points2,
	Eigen::VectorXd& vbeta, double *err2)
{
	Eigen::VectorXd vcol;
	Eigen::MatrixXd ejf;
	Eigen::MatrixXd ejft;
	Eigen::MatrixXd ejf2;
	Eigen::MatrixXd ejf_inv;
	Eigen::MatrixXd ejf3;
	Eigen::VectorXd ediff;

	ejf = jacobian(points2, vbeta);
	ejft = ejf.transpose();
	ejf2 = ejft * ejf;

	ejf_inv = ejf2.inverse();

	ejf3 = ejf_inv * ejft;

	vcol = get_r_col2(points2, vbeta);
	*err2 = get_err2(vcol);

	ediff = ejf3 * vcol;

	return ediff;
}


Eigen::VectorXd calc_beta(point_data_list2& points2, Eigen::VectorXd& vbeta)
{
	int j;
	int j1;
	int j2;
	double atten;
	double err2;
	double new_err2;
	Eigen::VectorXd ediff;
	Eigen::VectorXd ediff1;
	Eigen::VectorXd vbeta1;

	std::cout << "Here in calc_beta" << std::endl;

	for(j = 0;j < 100;j++)
	{
		std::cout << "iteration " << j << std::endl;

		std::cout << "VBETA" << std::endl;
		std::cout << vbeta << std::endl;
		std::cout << std::endl;

		ediff = calc_diff(points2, vbeta, &err2);

		std::cout << "err2 " << err2 << std::endl;

                ediff1 = ediff;

		atten = 1.0;
		for(j1 = 0;j1 < 10;j1++)
		{
			for(j2 = 0;j2 < ediff.rows();j2++)
			{
				ediff1(j2, 0) = atten * ediff(j2, 0);
			}

			vbeta1 = ediff1 + vbeta;

			new_err2 = calc_err2(points2, vbeta1);
			if(new_err2 < err2)
			{
				break;
			}
			else
			{
				atten *= 0.5;
			}
		}

                vbeta = vbeta1;

		if(new_err2 < 0.0000001)
		{
			std::cout << "solved" << std::endl;

			break;
		}

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


int main(int argc, char **argv)
{
	int i;
	point_data cur_point;
	point_data_list2 points2;

	double k = 1024.0;
	double tx = 10.0;
	double ty = 0.0;
	double tz = 10.0;
	double rx = 0.0;
	double ry = 45.0;
	double rz = 0.0;
	Eigen::VectorXd vbeta(7);
	Eigen::VectorXd vbeta1;

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

	vbeta1 = calc_beta(points2, vbeta);

	std::cout << "VBETA1" << std::endl;
	std::cout << vbeta1 << std::endl;
	std::cout << std::endl;

	PointDataListPrint4(points2, vbeta1);

	return 0;
}

