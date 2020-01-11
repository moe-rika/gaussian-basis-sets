#include <initializer_list>
#include <assert.h>
#include <iostream>
//#include <string>
#include <set>
#include <array>

//template <unsigned int N>
//class MyIter {
//public:
//	MyIter(const std::initializer_list<int>& _init)
//	{
//		assert(N > 0);
//		assert(_init.size() == N);
//		assert(*(_init.end() - 1) != 0);
//		int index = 0;
//		bool flag = true;
//		for (auto &a : _init)
//		{
//			if (a != 0)
//				flag = false;
//			assert(a >= 0);
//			cnt[index] = a - 1;
//			it[index] = a - 1;
//			index++;
//		}
//		m_is_zero = flag;
//	}
//	const bool is_zero() const
//	{
//		return m_is_zero;
//	}
//	const MyIter& operator --()
//	{
//		dec(0);
//		return *this;
//	}
//	const int& operator[](int i) const
//	{
//		return it[i];
//	}
//private:
//	int cnt[N];
//	int it[N];
//	bool m_is_zero;
//	void dec(int i)
//	{
//		if (!m_is_zero)
//		{
//			if (i == N - 1)
//			{
//				if (it[i] == 1)
//				{
//					it[i] = 0;
//					m_is_zero = true;
//				}
//				else
//				{
//					it[i]--;
//				}
//			}
//			else
//			{
//				if (it[i] != 0)
//				{
//					it[i]--;
//				}
//				else
//				{
//					it[i] = cnt[i];
//					dec(++i);
//				}
//			}
//		}
//	}
//};

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
			start[i] = std::get<1>(a);
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
			start[i] = std::get<1>(a);
			end[i] = std::get<2>(a);
			step[i] = std::get<3>(a);
		}
	}
	std::array<int*, N> ptr;
	std::array<int, N> start;
	std::array<int, N> end;
	std::array<int, N> step;
	std::array<int, N> current;
private:
	bool finish = false;
	void next()
	{
		if (!finish)
		{
			one_index_next(0);
		}
	}
	void one_index_next(int n)
	{
		if (n == N)
		{
			finish = true;
		}
		else
		{
			if (int res = current[n] + step[n] > end[n])
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

//int main()
//{
//	for (MyIter<3> iter = { 2,3,4 }; !iter.is_zero(); --iter)
//	{
//		std::cout << iter[0] << " " << iter[1] << " " << iter[2] << std::endl;
//	}
//	return 0;
//}

//
//class SimpleMultiThreadTter
//{
//public:
//	SimpleMultiThreadTter();
//	~SimpleMultiThreadTter();
//	class index
//private:
//
//};
//
//SimpleMultiThreadTter::SimpleMultiThreadTter()
//{
//}
//
//SimpleMultiThreadTter::~SimpleMultiThreadTter()
//{
//}