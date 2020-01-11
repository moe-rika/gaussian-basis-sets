#include <initializer_list>
#include <assert.h>
#include <iostream>
#include <tuple>
#include <array>


//provide mathematica style iter
template <unsigned int N>
struct MyRangeIter{
	MyRangeIter(const std::initializer_list<std::tuple<int*, int, int>>& _init)
	{
		assert(N > 0);
		assert(_init.size() == N);
		int i = 0;
		for (auto & a : _init)
		{
			assert(std::get<1>(a) - std::get<2>(a) <= 0);//start must <= end
			ptr[i] = std::get<0>(a);
			*(ptr[i]) = current[i] = start[i] = std::get<1>(a);
			end[i] = std::get<2>(a);
			step[i] = 1;
			i++;
		}
	}
	MyRangeIter(const std::initializer_list<std::tuple<int*, int, int, int>>& _init)
	{
		assert(N > 0);
		assert(_init.size() == N);
		int i = 0;
		for (auto & a : _init)
		{
			assert(std::get<1>(a) - std::get<2>(a) <= 0);//start must <= end
			ptr[i] = std::get<0>(a);
			*(ptr[i]) = current[i] = start[i] = std::get<1>(a);
			end[i] = std::get<2>(a);
			step[i] = std::get<3>(a);
			i++;
		}
	}

	void next()
	{
		if (!finish)
		{
			one_index_next(0);
		}
	}

	const bool is_finished() const { return finish; };

private:
	std::array<int*, N> ptr;
	std::array<int, N> start;
	std::array<int, N> end;
	std::array<int, N> step;
	std::array<int, N> current;

	bool finish = false;

	void one_index_next(int n)
	{
		if (n == N)
		{
			finish = true;
		}
		else
		{
			int res = current[n] + step[n];
			if (res  > end[n])
			{
				current[n] = start[n];
				*(ptr[n]) = start[n];
				one_index_next(n + 1);
			}
			else
			{
				current[n] = res;
				*(ptr[n]) = res;
			}
		}
	}
};

