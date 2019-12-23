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