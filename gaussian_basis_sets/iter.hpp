#include <initializer_list>
#include <assert.h>
#include <iostream>
#include <tuple>
#include <array>


//provide mathematica style iter
template <unsigned int N>
struct MyRangeIter{
	template< typename _Func>
	MyRangeIter(const std::initializer_list<std::tuple<int*, int, int>>& _init, const _Func& _f)
	{
		assert(N > 0);
		assert(_init.size() == N);
		int i = 0;
		for (auto & a : _init)
		{
			assert(std::get<1>(a) - std::get<2>(a) <= 0);//start must <= end
			ptr[i] = std::get<0>(a);
			*(ptr[i]) = start[i] = std::get<1>(a);
			end[i] = std::get<2>(a);
			step[i] = 1;
			i++;
		}
		one_index_next(0,_f);
	}

	template< typename _Func>
	MyRangeIter(const std::initializer_list<std::tuple<int*, int, int, int>>& _init,const _Func& _f)
	{
		assert(N > 0);
		assert(_init.size() == N);
		int i = 0;
		for (auto & a : _init)
		{
			assert(std::get<1>(a) - std::get<2>(a) <= 0);//start must <= end
			ptr[i] = std::get<0>(a);
			*(ptr[i]) = start[i] = std::get<1>(a);
			end[i] = std::get<2>(a);
			step[i] = std::get<3>(a);
			i++;
		}
		one_index_next(0,_f);
	}


private:
	std::array<int*, N> ptr;
	std::array<int, N> start;
	std::array<int, N> end;
	std::array<int, N> step;

	template< typename _Func>
	void one_index_next(int n, const _Func& _f)
	{
		if (n != N)
		{
			int res = *(ptr[n]) + step[n];
			if (res  > end[n])
			{
				*(ptr[n]) = start[n];
				one_index_next(n + 1,_f);
			}
			else
			{
				*(ptr[n]) = res;
				_f();
				one_index_next(0, _f);
			}
		}
	}
};

