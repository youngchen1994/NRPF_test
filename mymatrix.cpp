#include <iostream>
#include "math.h"
#include <algorithm>
#include <string>
#include <iomanip>
using namespace std;
#include "mymatrix.h"

matrix::matrix(const int &mm, const int &nn) :m(mm), n(nn)
{
	mat = new float*[m];
	for (int i = 0; i < m; i++)
		mat[i] = new float[n];
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mat[i][j] = 0.0F;
}
matrix::matrix() :m(0), n(0) {}
matrix::matrix(const int &mm, const int &nn, const float* arr) :m(mm), n(nn)
{
	mat = new float*[m];
	for (int i = 0; i < m; i++)
		mat[i] = new float[n];
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mat[i][j] = arr[i*n + j];
}
matrix::matrix(const int &mm, const int &nn, const float** arr) :m(mm), n(nn)
{
	mat = new float *[m];
	for (int i = 0; i < m; i++)
		mat[i] = new float[n];
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mat[i][j] = arr[i][j];
}

matrix::matrix(const matrix &c) :m(c.m),n(c.n)
{
	mat = new float *[m];
	for (int i = 0; i < m; i++)
		mat[i] = new float[n];
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			write(i, j, c.mat[i][j]);
		}
}

matrix::~matrix()
{
	for (int i = 0; i < m; i++)
		delete[] mat[i];
	delete []mat;

}

/*static matrix matrix::conv(const matrix &a, const matrix &b)
{
	int N = a.length(1), M = b.length(1);
	int NM = N + M;
	if (a.length(0) == 1)
	{
		matrix c(1, NM);
		for (int i = 0; i < NM; i++)
		{
			float val = 0;
			for (int j = 0; j < N; j++)
			{
				if (i - j < M)
					val += a.val(0, j)*b.val(0, i - j);
				else
					continue;
			}
			c.write(0, i, val);
		}
		return c;
	}
	else if (a.length(0) == 2)
	{
		matrix c(2, NM);
		for (int i = 0; i < NM; i++)
		{
			float val = 0;
			for (int j = 0; j < N; j++)
			{
				if (i - j < M)
					val += a.val(0, j)*b.val(0, i - j);
				else
					continue;
			}
			c.write(0, i, val);
			c.write(1, i, a.val(1, i) + b.val(1, i));
		}
		return c;
	}
	
}*/

void matrix::netsort(const int order[])
{

	int *Sorder = new int[m];
	for (int i = 0; i < m; i++)
	{
		Sorder[i] = order[i];
	}
	for (int i = 0; i < m; i++)
	{
		if (i + 1 != order[i])
		{
			for (int j = 0; j < m; j++)
				swap(Sorder[i] - 1, j, i, j);
			for (int j = 0; j < m; j++)
				swap(j, Sorder[i] - 1, j, i);
			for (int k = 0; k < m; k++)
			{
				if (Sorder[k] == i + 1)
				{
					Sorder[k] = Sorder[i];
				}
			}
			Sorder[i] = i + 1;
		}
	}
}



void matrix::write(int i, int j, float val)
{
	if (i == -1)
		i = m - 1;
	if (j == -1)
		j = n - 1;
	if (i*n + j < m*n)
		mat[i][j] = val;
}

void matrix::swap(int i, int j, int p, int q)
{
	float temp = mat[i][j];
	mat[i][j] = mat[p][q];
	mat[p][q] = temp;
}

int matrix::add(int i, int j, float val)
{
	if (i*n + j < m*n)
		mat[i][j] += val;
	else return -1;
	return 0;
}

float matrix::val(int i, int j) const
{
	if (i == -1)
		i = m - 1;
	if (j == -1)
		j = n - 1;
	if (i*n + j < m*n)
	{
		return mat[i][j];
	}
	else return -1;
}

void matrix::zero()
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mat[i][j] = 0.0F;
}

void matrix::renew(const int mm, const int nn)
{
	if (m > 0)
	{
		for (int i = 0; i < m; i++)
			delete[] mat[i];
		delete[]mat;

	}
	m = mm;
	n = nn;
	mat = new float *[m];
	for (int i = 0; i < m; i++)
		mat[i] = new float[n];
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mat[i][j] = 0.0F;
}

//inline int matrix::length(char flag) const
//{
	//if (flag == '0')
	//	return m;
	//if (flag == '1')
	//	return n;
//}

matrix matrix::combine(const matrix &b)
{
	matrix c(m, n + b.n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			c.write(i, j, val(i, j));
		}
		for (int j = n; j < n + b.n; j++) 
		{
			c.write(i, j, b.mat[i][j-n]);
		}
	}
	return c;
}


void matrix::show() const
{
	cout.unsetf(ios::fixed);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j <n; j++)
			std::cout<< setw(15)<<  setprecision(4) << val(i, j);
		std::cout << "\n";
	}
	std::cout << "\n";
}

float matrix::operator () (int i, int j) const
{
	if (i == -1)
		i = m - 1;
	if (j == -1)
		j = n - 1;
	if (i*n + j < m*n)
	{
		return mat[i][j];
	}
	
}

matrix matrix::operator + (const matrix &b) const
{
	matrix c(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c.mat[i][j] = mat[i][j] + b.mat[i][j];
	return c;
}

void matrix::operator += (matrix &b)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mat[i][j] = mat[i][j] + b.val(i,j);
}

matrix matrix::operator - (const matrix &b) const
{
	matrix c(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c.mat[i][j] = mat[i][j] - b.mat[i][j];
	return c;
}

void matrix::operator -= (matrix &b)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			mat[i][j] = mat[i][j] - b.val(i,j);
}

void matrix::operator = (matrix &b)
{
	if (m == b.m&&n == b.n)
	{
		for (int i = 0; i < b.m; i++)
			for (int j = 0; j < b.n; j++)
				mat[i][j] = b.val(i, j);
	}
	else
	{
		m = b.m;
		n = b.n;
		mat = new float*[b.m];
		for (int i = 0; i < b.m; i++)
			mat[i] = new float[b.n];
		for (int i = 0; i < b.m; i++)
			for (int j = 0; j < b.n; j++)
				mat[i][j] = b.val(i, j);
	}
	}

matrix matrix::operator * (const matrix &b) const
{
	if (n != b.m)
		return matrix();
	matrix c(m, b.n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < b.n; j++)
		{
			for (int k = 0; k < n; k++)
				c.add(i, j, val(i, k)*b.val(k, j));
		}
	}
	return c;
}

matrix matrix::operator * (const float &b)
{
	matrix c(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c.mat[i][j] = b*mat[i][j];
	return c;
}

float matrix::sum() const
{
	float sum = 0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			sum += mat[i][j];
	return sum;
}

float matrix::max() const
{
	float max = 0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (max < mat[i][j])
				max = mat[i][j];
	return fabs(max);
}

float matrix::distance() const
{
	float sum = 0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			sum += mat[i][j]*mat[i][j];
	return sqrt(sum);
}

matrix matrix::operator ^(const char &flag) const
{
	if (flag == 'I')
	{
		matrix c = *this;
		float fDet = 1.0F;
		int f = 1;
		int *iMax = new int[m];
		int *jMax = new int[m];
		for (int i = 0; i < m; i++)
		{
			jMax[i] = 0;
			iMax[i] = 0;
		}
		for (int k = 0; k < m; k++)
		{
			float fMax = 0.0F;
			for (int i = k; i < m; i++)
			{
				for (int j = k; j < m; j++)
				{
					float ff = abs(c.val(i, j));
					if (ff > fMax)
					{
						fMax = ff;
						iMax[k] = i;
						jMax[k] = j;
					}
				}
			}
			if (abs(fMax) < 0.00001F)
			{
				c.mdet = 0;
				cout << "No inverse matrix exists!!\n";
				return c;
			}

			if (iMax[k] != k)
			{
				f = -f;
				for (int count = 0; count < m; count++)
				{
					c.swap(k, count, iMax[k], count);
				}
			}
			if (jMax[k] != k)
			{
				f = -f;
				for (int count = 0; count < m; count++)
				{
					c.swap(count, k, count, jMax[k]);
				}
			}

			fDet *= c.val(k, k);

			c.write(k, k, 1.0F / c.val(k, k));

			for (int j = 0; j < m; j++)
			{
				if (j != k)
				{
					//MatA[k, j] *= MatA[k, k];
					c.write(k, j, c.val(k, k)*c.val(k, j));
				}
			}
			for (int i = 0; i < m; i++)
			{
				if (i != k)
				{
					for (int j = 0; j < m; j++)
					{
						if (j != k)
							//MatA[i, j] = MatA[i, j] - MatA[i, k] * MatA[k, j];
							c.write(i, j, c.val(i, j) - c.val(i, k)*c.val(k, j));
					}
				}
			}
			for (int i = 0; i < m; i++)
			{
				if (i != k)
				{
					//MatA[i, k] *= -MatA[k, k];
					c.write(i, k, -c.val(k, k)*c.val(i, k));
				}
			}
		}
		for (int k = m - 1; k >= 0; k--)
		{
			if (jMax[k] != k)
			{
				for (int count = 0; count < m; count++)
				{
					c.swap(k, count, jMax[k], count);
				}
			}
			if (iMax[k] != k)
			{
				for (int count = 0; count < m; count++)
				{
					c.swap(count, k, count, iMax[k]);
				}
			}
		}
		c.mdet = fDet * f;
		return c;
	}
	else if (flag == 'T')
	{
		matrix c(n, m);
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				c.write(j, i, val(i, j));
		return c;
	}
	matrix c(1,1);
	return c;
}


