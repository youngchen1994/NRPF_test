#pragma once
#include "mymatrix.h"
using namespace std;
#include <vector>
#include <iomanip>
#include <iostream>
#include <typeinfo>

//ϡ����󡪡�ʮ�������ʽ

//����Ԫ�ṹ��
struct sp_element{
	double U = 0, L = 0;//ֵ
	int i;//�к�
	int j;//�к�
	int right = -1;//��Ԫ���Ҳ��һ������Ԫ��λ��
	int down = -1;//��Ԫ���²��һ������Ԫ��λ��
};

//��ͷ/��ͷ����ڣ��ṹ��
struct HEAD{
	int count = 0;//����Ԫ����
	int first = -1;//������Ԫ��ַ
	bool flag;//ѡ�б�־
};//�С��п�ͷ��ʽ

//ϡ�������
class sparse
{
	bool symmetric;
	int m;//����
	int n;//����
	vector<sp_element> element;//����Ԫ����
	vector<HEAD> RHEAD, CHEAD;//�������
	vector<int> order;
	


	//�ڵ��Ż��ں˺�����������ŷ����仯�ķ���Ԫ���ƶ�������λ�á�
	void swap_inner_update(int &p1, double &temp, int &ni, int &nj, int &posR, int &posC);
	void swap_inner(const int &p, const int &q);

	//��Ž�������
	void swap(const int &p, const int &q);

	//Ѱַ����
	//���룺����Ԫ�ص����к�i��j��λ��q
	//������ҵ�Ԫ�أ���������true��q����ֵΪ����element�е�λ�ã�����Ϊfalse
	bool find(const int& i, const int &j, int &p, int& q);

public: 
	int current_position;//��ǰָ�����Ԫ��λ��
	bool position(const int &i, const int &j, int &posR)
	{
		int posC;
		if (find(i, j, posR, posC))
		{
			return true;
		}
		else
			return false;
	}

	//��ȡ����ԪA��i��j��
	double operator() (const int &i, const int &j);

	bool chooserow(const int &row);//ѡ������Ԫ
	bool right(int &i, int &j, int &r);//�����Ҳ���һ������Ԫ

	//չʾ����
	//���ϵ������
	void show();
	void clear();
	//���Ƿֽ�
	//ʹ��Tinney-2�������нڵ��Ż����
	//�붯̬������ȥ���������ٵĽڵ�
	//��Ŀ�����ֽ��L��D��U����ϡ�����
	void T2decomposition();

	//���Է��������
	//Ax=b
	void solve(vector<double> &b);

	//���º���
	//���룺ע��Ԫ�����к�i��j��ע��Ԫ��ֵval
	//��ע��Ԫ��ӵ�ϡ������У����¾���
	bool update(const int &i, const int &j, const double &val);

	//���캯��
	//����1����֪����ά�ȣ�m*n������յ�ϡ�����
	//����2����֪����ά��m*n������Ԫ���кż���ֵ���б�����ϡ�����
	sparse(const int &mm, const int &nn) :RHEAD(mm), CHEAD(nn), m(mm), n(nn), order(mm) {};
	sparse(const int &mm, const int &nn, const matrix &list, const bool sym);
};

