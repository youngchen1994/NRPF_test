#ifndef _MYMATRIX_H
#define _MYMATRIX_H
#include <algorithm>
template <const int N>
bool comp(float* a, float* b)
{
	return (a[N] < b[N]);
}

class matrix //������
{
protected:
	int m;
	int n;
	float **mat;
	float mdet;
	friend class charge_system_search;

public:
	//static matrix conv(const matrix &a, const matrix &b);
	matrix(const int &mm, const int &nn);//���캯��
	matrix();//���캯��
	matrix(const int &mm, const int &nn, const float* arr);//���캯��
	matrix(const int &mm, const int &nn, const float** arr);//���캯��
	matrix(const matrix &c);//���캯��
	~matrix();//��������

	void renew(const int mm, const int nn);

	void netsort(const int order[]);

	template<int N>
	void colsort()
	{
		sort(mat, mat + m, comp<N>);
	}

	void write(int i, int j, float val);//����Ԫ��д��
	void swap(int i, int j, int p, int q);//����Ԫ�ؽ���
	int add(int i, int j, float val);//����Ԫ�ؼӵ���
	void zero();//����
	matrix combine(const matrix &b);
	float val(int i, int j) const;//�ⲿ����Ԫ��ֵ���ɣ�
	
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
	void show() const;//��ӡ����
	
//���������
	float operator () (int i, int j) const;//����Ԫ��ֵ���£�
	matrix operator + (const matrix &b) const;
	void operator += (matrix &b);
	matrix operator - (const matrix &b) const;
	void operator -= (matrix &b);
	void operator = (matrix &b);
	matrix operator * (const matrix &b) const;//����˷�
	matrix operator * (const float &b);//����

	matrix operator ^(const char &flag) const;//���棨I����ת�ã�T��

};

#endif
