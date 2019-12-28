#pragma once
//C++ 98 required

const int max_order = 15;
int table[max_order][max_order + 1];

template <unsigned int N>
struct Factorial
{
	static const int value = N * Factorial<N - 1>::value;
};

template <>
struct Factorial<0>
{
	static const int value = 1;
};

template <unsigned int M, unsigned int N>
struct ValueTable
{
	ValueTable()
	{
		table[M][N] = Factorial<M>::value / Factorial<M - N>::value / Factorial<N>::value;
	}
	ValueTable<M, N - 1> next;
};

template <unsigned int M>
struct ValueTable<M, 0>
{
	ValueTable()
	{
		table[M][0] = Factorial<M>::value / Factorial<M - 0>::value / Factorial<0>::value;
	}
	ValueTable<M - 1, M - 1> next;
};

template <>
struct ValueTable<0, 0>
{
	ValueTable()
	{
		table[0][0] = Factorial<0>::value / Factorial<0 - 0>::value / Factorial<0>::value;
	}
	//static const int next = ValueTable<M - 1, M - 1>::value();
};

