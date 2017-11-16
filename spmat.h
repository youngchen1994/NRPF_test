#pragma once
#include "mymatrix.h"
#include <vector>
#include <iomanip>
using namespace std;

class spmat
{
	int m;
	int n;
	float Err = 1e-4;
	vector<float> U;//�д��������Ƿ���Ԫֵ
	vector<int> JU;//�д��������Ƿ���Ԫ�к�
	vector<int> IU;//�д���������ÿ�е�һ������Ԫ��U�е�λ��
	vector<int> Ulink;//��һ������Ԫ��U�е�λ��
	vector<int> Una;//ÿ�з���Ԫ����
	vector<float> D;
	vector<float> L;
	vector<int> JL;
	vector<int> IL;
	vector<int> Llink;
	vector<int> Lna;

public:
	//������ʾ
	void show()
	{
		int q = 0, ql=-1, qr=-1;
		for (int i = 0; i <m; i++)
		{
			for (int j = 0; j <n; j++)
			{
				if (i == j)
					cout << setprecision(4) << D[i] << "\t";
				else if (i < j)
				{
					if (find(i, j, q, ql, qr) > 0)
						cout << setprecision(4) << U[q] << "\t";
					else
						cout << 0 << "\t";
				}
				else
				{
					if (find(i, j, q, ql, qr) > 0)
						cout << setprecision(4) << L[q] << "\t";
					else
						cout << 0 << "\t";
				}
			}
			cout << "\n";
		}
	}

	//���Է��������
	matrix equsolve(const matrix &b)
	{
		matrix z(b);
		int q = 0, ql = -1, qr = -1;
		for (int i = 0; i < n - 1; i++)
		{
			if (z(i, 0) == 0)
				continue;
			for (int j = i + 1; j < n; j++)
			{
				if (find(j, i, q, ql, qr) > 0)
					z.add(j, 0, -L[q] * z(i, 0));
			}
		}
		matrix y(z);
		for (int i = 0; i < n; i++)
			y.write(i, 0, y(i, 0) / D[i]);
		matrix x(y);
		for (int j = n - 1; j > 0; j--)
			for (int i = j - 1; i >= 0; i--)
			{
				if (find(i, j, q, ql, qr) > 0)
					x.add(i, 0, -U[q] * x(j, 0));
			}
		return x;
	}

	//�����ж�λ����Ԫ������ֵ
	int find(int i, int j, int &q, int &ql, int &qr)
	{
		if (i < j)
		{
			if (Una[i] == 0)
				return 0;
			for (int k = IU[i]; k != -1; k = Ulink[k])
			{
				if (JU[k] == j)
				{
					q = k;
					return 1;
				}
				else if (JU[k] < j)
					ql = k;
				else
					qr = k;
			}
			return -1;
		}
		else
		{
			if (Lna[j] == 0)
				return 0;
			for (int k = JL[j]; k != -1; k = Llink[k])
			{
				if (IL[k] == i)
				{
					q = k;
					return 1;
				}
				else if (IL[k] < i)
					ql = k;
				else
					qr = k;
			}
			return -1;
		}
	}
	//���Ƿֽ�
	void decomposition()
	{
		int i;
		int j;
		for (int p = 0; p < n - 1; p++)//ѡ��n-1��
			for (int k = IU[p]; k != -1; k = Ulink[k])//ѡ��������Ԫ��
			{
				U[k] = U[k] / D[p];
				j = JU[k];
				for (int l = JL[p]; l != -1; l = Llink[l])
				{
					i = IL[l];
					if (i == j)
						D[i] -= U[k] * L[l];
					else if (i < j)
					{
						int q = 0, ql = -1, qr = -1;
						int jud = find(i, j, q, ql, qr);
						if (jud > 0)
							U[q] -= U[k] * L[l];
						else if (jud == 0)
						{
							Una[i]++;
							U.push_back(-U[k] * L[l]);
							JU.push_back(j);
							Ulink.push_back(-1);
							IU[i] = U.size() - 1;
						}
						else
						{
							Una[i]++;
							U.push_back(-U[k] * L[l]);
							JU.push_back(j);
							if (ql < 0)//��һ��Ԫ��
							{
								Ulink.push_back(IU[i]);
								IU[i] = U.size() - 1;
							}
							else if (qr < 0)//���һ��Ԫ��
							{
								Ulink.push_back(-1);
								Ulink[ql] = U.size() - 1;
							}
							else
							{
								Ulink[ql] = U.size() - 1;
								Ulink.push_back(qr);
							}
						}
					}
					else
					{
						int q = 0, ql = -1, qr = -1;
						int jud = find(i, j, q, ql, qr);
						if (jud > 0)
							L[q] -= U[k] * L[l];
						else if (jud == 0)
						{
							Lna[j]++;
							L.push_back(-U[k] * L[l]);
							IL.push_back(i);
							Llink.push_back(-1);
							JL[j] = L.size() - 1;
						}
						else
						{
							Lna[j]++;
							L.push_back(-U[k] * L[l]);
							IL.push_back(i);
							if (ql < 0)//��һ��Ԫ��
							{
								Llink.push_back(JL[j]);
								JL[j] = L.size() - 1;
							}
							else if (qr < 0)//���һ��Ԫ��
							{
								Llink.push_back(-1);
								Llink[ql] = L.size() - 1;
							}
							else
							{
								Llink[ql] = L.size() - 1;
								Llink.push_back(qr);
							}
						}
					}
				}
				show();
				cout << endl;
			}
		for (int p = 0; p < n - 1; p++)
			for (int l = JL[p]; l < JL[p] + Lna[p]; l++)
				L[l] /= D[p];
	}
	//���캯������ͨ����==��ϡ�����
	spmat(const matrix &a)
	{
		m = a.length(0), n = a.length(1);
		for (int i = 0; i < m - 1; i++)
		{
			Una.push_back(0);
			for (int j = i + 1; j < n; j++)
			{
				if (fabs(a(i, j)) < Err)
					continue;
				U.push_back(a(i, j));
				JU.push_back(j);
				if (Una[i] == 0)
					IU.push_back(U.size() - 1);
				
				Una[i]++;
				Ulink.push_back(U.size());
			}
			if (Una[i] != 0)
			{
				Ulink.back() = -1;
			}
			else
				IU.push_back(U.size());
		}
		for (int j = 0; j < n - 1; j++)
		{
			Lna.push_back(0);
			for (int i = j + 1; i < m; i++)
			{
				if (fabs(a(i, j)) < Err)
					continue;
				L.push_back(a(i, j));
				IL.push_back(i);
				if (Lna[j] == 0)
					JL.push_back(L.size() - 1);
				
				Lna[j]++;
				Llink.push_back(L.size());
			}
			if (Lna[j] != 0)
			{
				Llink.back() = -1;
			}
			else
				JL.push_back(L.size());
		}
		for (int i = 0; i < m; i++)
			D.push_back(a(i, i));
		IU.push_back(IU.back());
		JL.push_back(JL.back());


	}
};
