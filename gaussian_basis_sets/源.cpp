#include <cmath>
//erf() needs C++11
#include <iostream>
#include <array>
#include <map>


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

struct ExGauBas
{// x^a*y^b*z^c*ker
	GauBas ker;
	int a, b, c;
	ExGauBas operator*(const ExGauBas& egb)
	{
		return ExGauBas{
			ker*egb.ker,
			a + egb.a,b + egb.b,c + egb.c
		};
	}
	double operator()(const Point3D&p) const
	{
		auto R = p - ker.center;
		return pow(R.x, a)*pow(R.y, b)*pow(R.z, c)*ker(p);
	}
};

double term_T(const GauBas& gb1, const GauBas& gb2)
{
	double term_var_v_sq = gb1.alpha*gb2.alpha / (gb1.alpha + gb2.alpha);
	return (gb1.center - gb2.center).norm_sq()*term_var_v_sq;
}

double term_perfix_co(const GauBas& gb1, const GauBas& gb2)
{
	return gb1.k * gb2.k * pow(PI, 3) / (gb1.alpha*gb2.alpha*sqrt(gb1.alpha + gb2.alpha));
}

double calc_integral(const GauBas& gb1, const GauBas& gb2)
{
	double a = term_T(gb1, gb2);
	if (abs(a) < 0.0001)
	{
		return term_perfix_co(gb1, gb2) * 2 / sqrt(0 + PI);
	}
	else
	{
		return term_perfix_co(gb1, gb2) / sqrt(a)*erf(sqrt(a));
	}
}

double calc_overlap(const GauBas& gb1, const GauBas& gb2)
{
	return gb1.k*gb2.k*pow(PI / (gb1.alpha + gb2.alpha), 1.5)*
		exp(-gb1.alpha*gb2.alpha*(gb1.center - gb2.center).norm_sq() /
		(gb1.alpha + gb2.alpha)
		);
}

template <unsigned int N>
struct HighDimSimpsonIntegration
{
	struct IntegrationDescriptor
	{
		int div;
		double min, max, val;
	};
	static double simpson_int(
		std::array<IntegrationDescriptor, N>& d,
		unsigned int i,
		double(*f)(std::array<IntegrationDescriptor, N>&)
	)
	{
		int &div = d[i].div;
		double &min = d[i].min, &max = d[i].max, &val = d[i].val;
		if (i != N - 1)
		{
			double sum = 0;
			val = min;
			sum += simpson_int(d, i + 1, f);
			val = max;
			sum += simpson_int(d, i + 1, f);
			double step = (max - min) / div;
			for (int j = 1; j < div; j += 2)
			{
				val = min + step * j;
				sum += 4 * simpson_int(d, i + 1, f);
				val = min + step * (j + 1);
				sum += 2 * simpson_int(d, i + 1, f);
				if (i == 0)
					std::cout << j << std::endl;
			}
			return (max - min) / 3 / div * sum;
		}
		else
		{
			double sum = 0;
			val = min;
			sum += f(d);
			val = max;
			sum += f(d);
			double step = (max - min) / div;
			for (int j = 1; j < div; j += 2)
			{
				val = min + step * j;
				sum += 4 * f(d);
				val = min + step * (j + 1);
				sum += 2 * f(d);
			}
			return (max - min) / 3 / div * sum;
		}
	}
};

double func(std::array<HighDimSimpsonIntegration<3>::IntegrationDescriptor, 3>& d)
{
	if (d[0].val*d[0].val + d[1].val*d[1].val + d[2].val*d[2].val > 1)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

GauBas a{ 0.3696,0.4166,{-0.4,0,0} }, b{ 0.5881,0.7739,{0.4,0,0} };
GauBas c = a * b, dd = a * a, e = b * b;

double func1(std::array<HighDimSimpsonIntegration<6>::IntegrationDescriptor, 6>& d)
{
	Point3D p{ d[0].val,d[1].val,d[2].val }, q{ d[3].val,d[4].val,d[5].val };
	double R = (p - q).norm();
	return c(p)*c(q) / R;
}

double func2(std::array<HighDimSimpsonIntegration<6>::IntegrationDescriptor, 6>& d)
{
	Point3D p{ d[0].val,d[1].val,d[2].val }, q{ d[3].val,d[4].val,d[5].val };
	double R = (p - q).norm();
	return dd(p)*c(q) / R;
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

double partial_calc_integral(const std::array<GauBas, 4>& a, int index_of_integrand, int choose_of_xyz)
{
	const double delta = 0.00000001;
	GauBas da = a[index_of_integrand];
	da.center[choose_of_xyz] += delta;
	double val = calc_integral(a[0] * a[1], a[2] * a[3]);
	if (index_of_integrand == 0)
		return (calc_integral(da*a[1], a[2] * a[3]) - val) / delta;
	else if (index_of_integrand == 1)
		return (calc_integral(a[0] * da, a[2] * a[3]) - val) / delta;
	else if (index_of_integrand == 2)
		return (calc_integral(a[0] * a[1], da * a[3]) - val) / delta;
	else
		return (calc_integral(a[0] * a[1], a[2] * da) - val) / delta;
}

int main()
{
	const int div = 40, dim = 6;
	const double range = 8;
	//std::array<HighDimSimpsonIntegration<dim>::IntegrationDescriptor, dim> d{
	//	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
	//{div,-range + range / div / 2,range + range / div / 2,0},
	//	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
	//{div,-range + range / div / 2,range + range / div / 2,0},
	//	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
	//{div,-range + range / div / 2,range + range / div / 2,0},
	//	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
	//{div,-range,range,0},
	//	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
	//{div,-range,range,0},
	//	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
	//{div,-range,range,0}
	//};
	//double k = HighDimSimpsonIntegration<6>::simpson_int(d, 0, func1);
	//std::cout << k << std::endl;
	//std::cout << calc_integral(a*b, a*b) << std::endl;
	//double k = HighDimSimpsonIntegration<6>::simpson_int(d, 0, func2);
	//std::cout << k << std::endl;
	//std::cout << calc_integral(a*a, a*b) << std::endl;


	G(30, 0.5);

	return 0;
}

