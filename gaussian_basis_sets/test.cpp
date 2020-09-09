#include <cmath>
//erf() needs C++11
#include <iostream>
#include <array>
#include <map>
#include <vector>
#include "simpson.hpp"
#include "high_order_partial.hpp"

#include "iter.hpp"

const double PI = 3.141592653589793;

struct Point3D
{
	double x, y, z;
	Point3D operator-(const Point3D& p) const
	{
		return Point3D{ x - p.x ,y - p.y, z - p.z };
	}

	Point3D operator+(const Point3D& p)const
	{
		return Point3D{ x + p.x ,y + p.y, z + p.z };
	}

	friend Point3D operator*(const double& a, const Point3D& b)
	{
		return Point3D{ a*b.x ,a*b.y, a*b.z };
	}

	friend Point3D operator*(const Point3D& a, const double& b)
	{
		return b * a;
	}

	friend bool operator*(const Point3D& a, const Point3D& b)
	{
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}

	friend Point3D operator/(const Point3D& a, const double& b)
	{
		return Point3D{ a.x / b ,a.y / b, a.z / b };
	}

	double& operator[](int a)
	{
		if (a == 0)
			return x;
		if (a == 1)
			return y;
		return z;
	}

	const double& operator[](int a) const
	{
		if (a == 0)
			return x;
		if (a == 1)
			return y;
		return z;
	}

	const double norm_sq() const
	{
		return x * x + y * y + z * z;
	}

	const double norm() const
	{
		return sqrt(norm_sq());
	}

};



struct GauBas
{//gaussian basis k*E^(-alpha*r^2)
	double k;
	double alpha;
	Point3D center;
	GauBas operator*(const GauBas& gb) const
	{
		return GauBas{
			k*gb.k*exp(-alpha * gb.alpha / (alpha + gb.alpha)*(center - gb.center).norm()),
			alpha + gb.alpha ,
			(alpha*center + gb.alpha*gb.center) / (alpha + gb.alpha)
		};
	}
	double operator()(const Point3D&p) const
	{
		return k * exp(-alpha * (center - p).norm_sq());
	}
};


struct TwoExGauBas
{// x^a*y^b*z^c*ker
	GauBas ker;
	Point3D center1, center2;
	int a1, b1, c1, a2, b2, c2;
	double operator()(const Point3D&p) const
	{
		auto R1 = p - center1;
		auto R2 = p - center2;
		return pow(R1.x, a1)*pow(R1.y, b1)*pow(R1.z, c1)*
			pow(R2.x, a2)*pow(R2.y, b2)*pow(R2.z, c2)*ker(p);
	}
};

struct ExGauBas
{// x^a*y^b*z^c*ker
	GauBas ker;
	int a, b, c;
	friend TwoExGauBas operator*(const ExGauBas& egb1, const ExGauBas& egb2)
	{
		return TwoExGauBas{
			egb1.ker*egb2.ker,
			egb1.ker.center,
			egb2.ker.center,
			egb1.a,egb1.b,egb1.c,
			egb2.a,egb2.b,egb2.c
		};
	}
	double operator()(const Point3D&p) const
	{
		auto R = p - ker.center;
		return pow(R.x, a)*pow(R.y, b)*pow(R.z, c)*ker(p);
	}
};

template <typename T, unsigned int N>
struct NDVect
{
	typedef  std::vector<typename NDVect<T, N - 1>::type> type;
};
template <typename T>
struct NDVect<T, 1>
{
	typedef std::vector<T> type;
};

class FourTermIntegral
{
public:
	FourTermIntegral(const ExGauBas& _A, const ExGauBas& _B,
		const ExGauBas& _C, const ExGauBas& _D) :m_A(_A), m_B(_B), m_C(_C), m_D(_D)
	{/*
		m_size_list.push_back(m_A.a);
		m_size_list.push_back(m_A.b);
		m_size_list.push_back(m_A.c);
		m_size_list.push_back(m_B.a);
		m_size_list.push_back(m_B.b);
		m_size_list.push_back(m_B.c);
		m_size_list.push_back(m_C.a);
		m_size_list.push_back(m_C.b);
		m_size_list.push_back(m_C.c);
		m_size_list.push_back(m_D.a);
		m_size_list.push_back(m_D.b);
		m_size_list.push_back(m_D.c);*/
	};

	double calc_integral_ker()
	{
		double a = term_T_ker();
		if (abs(a) < 0.0001)
		{
			return term_perfix_co_ker() * 2 / sqrt(0 + PI);
		}
		else
		{
			return term_perfix_co_ker() / sqrt(a)*erf(sqrt(a));
		}
	}

private:
	ExGauBas m_A, m_B, m_C, m_D;
	GauBas ker_list[4] = { m_A.ker,m_B.ker,m_C.ker,m_D.ker };
	GauBas gb1 = ker_list[0] * ker_list[1];
	GauBas gb2 = ker_list[2] * ker_list[3];
	/*
		std::vector<int> m_size_list;*/


	NDVect<double, 12>::type m_partial_derivate_val;
	NDVect<double, 12>::type m_ker_fn_val;


	double term_T_ker()
	{
		double term_var_v_sq = gb1.alpha*gb2.alpha / (gb1.alpha + gb2.alpha);
		return (gb1.center - gb2.center).norm_sq()*term_var_v_sq;
	}

	double term_perfix_co_ker()//ker means it is not consider x^a*y^b*z^c term
	{
		return gb1.k * gb2.k * pow(PI, 3) /
			(gb1.alpha*gb2.alpha*sqrt(gb1.alpha + gb2.alpha));
	}



	double calc_overlap_ker()
	{
		return gb1.k*gb2.k*pow(PI / (gb1.alpha + gb2.alpha), 1.5)*
			exp(-gb1.alpha*gb2.alpha*(gb1.center - gb2.center).norm_sq() /
			(gb1.alpha + gb2.alpha)
			);
	}


	//double partial_calc_integral(int &index_of_integrand, int &choose_of_xyz)
	//{
	//	const double delta = 0.00000001;
	//	double val1 = calc_integral_ker();
	//	GauBas &da = ker_list[index_of_integrand];
	//	da.center[choose_of_xyz] += delta;
	//	gb1 = ker_list[0] * ker_list[1];
	//	gb2 = ker_list[2] * ker_list[3];
	//	double val2 = calc_integral_ker();
	//	da.center[choose_of_xyz] -= delta;
	//	gb1 = ker_list[0] * ker_list[1];
	//	gb2 = ker_list[2] * ker_list[3];
	//	return (val2 - val1) / delta;

	//}

	double high_order_partial(int &index_of_integrand, int &choose_of_xyz, int order)
	{
		if (order == 0)
		{
			return calc_integral_ker();
		}
		else if (order == 1)
		{
			const double delta = 0.00000001;
			double val1 = calc_integral_ker();
			GauBas &da = ker_list[index_of_integrand];
			da.center[choose_of_xyz] += delta;
			gb1 = ker_list[0] * ker_list[1];
			gb2 = ker_list[2] * ker_list[3];
			double val2 = calc_integral_ker();
			da.center[choose_of_xyz] -= delta;
			gb1 = ker_list[0] * ker_list[1];
			gb2 = ker_list[2] * ker_list[3];
			return (val2 - val1) / delta;
		}
		else
		{
			const double delta = 0.00000001;
			double val1 = high_order_partial(index_of_integrand, choose_of_xyz, order - 1);
			GauBas &da = ker_list[index_of_integrand];
			da.center[choose_of_xyz] += delta;
			gb1 = ker_list[0] * ker_list[1];
			gb2 = ker_list[2] * ker_list[3];
			double val2 = high_order_partial(index_of_integrand, choose_of_xyz, order - 1);
			da.center[choose_of_xyz] -= delta;
			gb1 = ker_list[0] * ker_list[1];
			gb2 = ker_list[2] * ker_list[3];
			return (val2 - val1) / delta;
		}
	}

	std::map<int, double> G_map;// order of derivation | coef

	void G_helper(int m, int n, const double& alpha, double coef)
	{
		if (m == 0)
		{
			G_map[n] += coef;
		}
		else if (m == 1)
		{
			G_helper(m - 1, n + 1, alpha, -coef / 2 / alpha);
		}
		else
		{
			G_helper(m - 2, n, alpha, coef * m / 2 / alpha);
			G_helper(m - 1, n + 1, alpha, -coef / 2 / alpha);
		}
	}

	void G(int n, const double& alpha)
	{
		G_helper(n, 0, alpha, 1);
	}



};

//
//double partial_calc_integral(const std::array<GauBas, 4>& a, int index_of_integrand, int choose_of_xyz)
//{
//	const double delta = 0.00000001;
//	GauBas da = a[index_of_integrand];
//	da.center[choose_of_xyz] += delta;
//	double val = calc_integral(a[0] * a[1], a[2] * a[3]);
//	if (index_of_integrand == 0)
//		return (calc_integral(da*a[1], a[2] * a[3]) - val) / delta;
//	else if (index_of_integrand == 1)
//		return (calc_integral(a[0] * da, a[2] * a[3]) - val) / delta;
//	else if (index_of_integrand == 2)
//		return (calc_integral(a[0] * a[1], da * a[3]) - val) / delta;
//	else
//		return (calc_integral(a[0] * a[1], a[2] * da) - val) / delta;
//
//}
//int main()
//{


	//GauBas a{ 0.3696,0.4166,{-0.4,0,0} }, b{ 0.5881,0.7739,{0.4,0,0} };
	////GauBas c = a * b, dd = a * a, e = b * b;
	//ExGauBas p{ a,0,0,0 }, q{ b,0,0,0 };
	//FourTermIntegral fti(p, q, p, q);
	//std::cout << fti.calc_integral_ker() << std::endl;
	//std::cout << log(tgamma(10)) << std::endl;

	//int a, b, c;
	//MyRangeIter<3> iter
	//(
	//	{{&a,-5,4,1},{&b,11,16,1},{&c,12,13,1}},
	//	[&]() {std::cout << a << "," << b << "," << c << std::endl; }
	//);

	//NDVect<double, 12>::type m;

	//G(30, 0.5);
	//int k = count_micro_state(10, 10);

	//return 0;
//}


