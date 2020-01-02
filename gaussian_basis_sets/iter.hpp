#include <initializer_list>
#include <assert.h>
#include <iostream>
#include <string>

template <unsigned int N>
class MyIter {
public:
	MyIter(const std::initializer_list<int>& _init)
	{
		assert(N > 0);
		assert(_init.size() == N);
		assert(*(_init.end() - 1) != 0);
		int index = 0;
		bool flag = true;
		for (auto &a : _init)
		{
			if (a != 0)
				flag = false;
			assert(a >= 0);
			cnt[index] = a - 1;
			it[index] = a - 1;
			index++;
		}
		m_is_zero = flag;
	}
	const bool is_zero() const
	{
		return m_is_zero;
	}
	const MyIter& operator --()
	{
		dec(0);
		return *this;
	}
	const int& operator[](int i) const
	{
		return it[i];
	}
private:
	int cnt[N];
	int it[N];
	bool m_is_zero;
	void dec(int i)
	{
		if (!m_is_zero)
		{
			if (i == N - 1)
			{
				if (it[i] == 1)
				{
					it[i] = 0;
					m_is_zero = true;
				}
				else
				{
					it[i]--;
				}
			}
			else
			{
				if (it[i] != 0)
				{
					it[i]--;
				}
				else
				{
					it[i] = cnt[i];
					dec(++i);
				}
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

