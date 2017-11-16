#ifndef _MYMATRIX_H
#define _MYMATRIX_H
#include <algorithm>
template <const int N>
bool comp(float* a, float* b)
{
	return (a[N] < b[N]);
}

class matrix //矩阵类
{
protected:
	int m;
	int n;
	float **mat;
	float mdet;
	friend class charge_system_search;

public:
	//static matrix conv(const matrix &a, const matrix &b);
	matrix(const int &mm, const int &nn);//构造函数
	matrix();//构造函数
	matrix(const int &mm, const int &nn, const float* arr);//构造函数
	matrix(const int &mm, const int &nn, const float** arr);//构造函数
	matrix(const matrix &c);//构造函数
	~matrix();//析构函数

	void renew(const int mm, const int nn);

	void netsort(const int order[]);

	template<int N>
	void colsort()
	{
		sort(mat, mat + m, comp<N>);
	}

	void write(int i, int j, float val);//矩阵元素写入
	void swap(int i, int j, int p, int q);//矩阵元素交换
	int add(int i, int j, float val);//矩阵元素加等于
	void zero();//清零
	matrix combine(const matrix &b);
	float val(int i, int j) const;//外部调用元素值（旧）
	
	int length(int flag) const
	{
		if (flag == 0)
			return m;
		else if (flag == 1)
			return n;
		else
			return 0;

	}
	float sum() const;
	float max() const;
	float distance() const;
	void show() const;//打印矩阵
	
//重载运算符
	float operator () (int i, int j) const;//调用元素值（新）
	matrix operator + (const matrix &b) const;
	void operator += (matrix &b);
	matrix operator - (const matrix &b) const;
	void operator -= (matrix &b);
	void operator = (matrix &b);
	matrix operator * (const matrix &b) const;//矩阵乘法
	matrix operator * (const float &b);//数乘

	matrix operator ^(const char &flag) const;//求逆（I）和转置（T）

};

#endif
