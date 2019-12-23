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




//
//
//double func(std::array<HighDimSimpsonIntegration<3>::IntegrationDescriptor, 3>& d)
//{
//	if (d[0].val*d[0].val + d[1].val*d[1].val + d[2].val*d[2].val > 1)
//	{
//		return 0;
//	}
//	else
//	{
//		return 1;
//	}
//}
//
//GauBas a{ 0.3696,0.4166,{-0.4,0,0} }, b{ 0.5881,0.7739,{0.4,0,0} };
//GauBas c = a * b, dd = a * a, e = b * b;
/*
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
*/

//const int div = 40, dim = 6;
//const double range = 8;
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