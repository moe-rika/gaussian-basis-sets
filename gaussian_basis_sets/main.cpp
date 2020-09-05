#include <cmath>
#include <vector>
#include <deque>


const int max_angular_momentum = 150;

using std::deque;
using std::vector;

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
		return Point3D{ a * b.x ,a * b.y, a * b.z };
	}

	friend Point3D operator*(const Point3D& a, const double& b)
	{
		return b * a;
	}

	friend bool operator*(const Point3D& a, const Point3D& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
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

	const double sum() const
	{
		return x + y + z;
	}
};


//unnormalized atom-centered GTOs

class GaussianTypeOrbital {
public:
	struct AngularMomentum
	{
		int a[3];//a[0] x;

		int& operator[](int _n)
		{
			return a[_n];
		}

		const int& operator[](int _n) const
		{
			return a[_n];
		}
	};
	AngularMomentum a;
	double alpha;
	Point3D center;
};//needs to be improved AoS and SoA


//J. Chem. Phys. 89, 5777 (1988); doi: 10.1063/1.455553
//Phys. Chem. Chem. Phys., 2006, 8, 3072¨C3077
class TwoElectronIntegral {
public:
	TwoElectronIntegral(const GaussianTypeOrbital gto[4])
	{
		const auto& ia = gto[0];
		const auto& ib = gto[1];
		const auto& ic = gto[2];
		const auto& id = gto[3];

		const auto& alpha = ia.alpha;
		const auto& beta = ib.alpha;
		const auto& gamma = ic.alpha;
		const auto& delta = id.alpha;


		a = ia.a; A = ia.center;
		b = ib.a; B = ib.center;
		c = ic.a; C = ic.center;
		d = id.a; D = id.center;

		zeta = alpha + beta;
		eta = gamma + delta;
		rho = zeta * eta / (zeta + eta);
		P = (alpha * A + beta * B) / zeta;
		Q = (gamma * C + delta * D) / eta;
		T = rho * (P - Q).norm_sq();
		S_ab = exp(-alpha * beta / zeta * (A - B).norm_sq());
		S_cd = exp(-gamma * delta / eta * (C - D).norm_sq());
		k = pow((PI / (zeta + eta)), 1.5);
	}
	GaussianTypeOrbital::AngularMomentum a, b, c, d;
	Point3D A, B, C, D;//used in next step,copy of 


	double zeta;
	double eta;
	///////////////////////
	double rho;
	///////////////////////
	Point3D P, Q;
	double T;
	double S_ab, S_cd;
	double k;
};


class TwoElectronIntegralHelper {
public:
	TwoElectronIntegralHelper(TwoElectronIntegral* _ei) :ei(_ei)
	{
		recurrence_queue.push_front({ 1,ei->a,ei->b,ei->c,ei->d,0 });
		delta_AB = ei->A - ei->B;
		delta_CD = ei->C - ei->D;
	};
	TwoElectronIntegral* ei;
	Point3D delta_AB, delta_CD;
	deque<TwoElectronIntegralHelperUnit> recurrence_queue, recurrence_queue_A;
	void horizontal_recurrence()
	{
		for (int i = 0; i < 3; i++)
		{
			while (!recurrence_queue.empty())
			{
				TwoElectronIntegralHelperUnit& a = recurrence_queue.front();
				if (a.b[i] == 0)
				{
					recurrence_queue_A.push_front(a);
					recurrence_queue.pop_front();
				}
				else
				{
					TwoElectronIntegralHelperUnit a1 = a;
					TwoElectronIntegralHelperUnit a2 = a;
					a1.a[i] += 1;
					a1.b[i] -= 1;

					a2.b[i] -= 1;
					a2.k *= delta_AB[0];

					recurrence_queue.push_back(a1);
					recurrence_queue.push_back(a2);

					recurrence_queue.pop_front();
				}
			}
			recurrence_queue.swap(recurrence_queue_A);
		}

	}

	vector<double> m{ max_angular_momentum };
};

struct TwoElectronIntegralHelperUnit
{
	double k = 1;
	GaussianTypeOrbital::AngularMomentum a, b, c, d;
	int m = 0;
};



int main()
{
	vector<int> a;
	a.resize(1);
	return 0;
}