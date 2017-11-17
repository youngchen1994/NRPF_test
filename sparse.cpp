#include "mymatrix.h"
using namespace std;
#include <vector>
#include <iomanip>
#include "sparse.h"

//�ڵ��Ż��ں˺�������������Ϊ�������仯����Ԫ�ƶ�������λ�á�
void sparse::swap_inner_update(int &p1, double &temp, int &ni, int &nj, int &posR, int &posC)
{
	if (find(ni, nj, posR, posC))//�ҵ��Ի�Ԫ
	{

		if (ni < nj)
		{
			if (p1 != posR)
			{
				temp = element[p1].U;
				element[p1].U = element[posR].U;
				element[posR].U = temp;
				temp = element[p1].L;
				element[p1].L = element[posR].L;
				element[posR].L = temp;
			}
		}
		else//��������ֵ����
		{
			if (p1 != posR)
			{
				temp = element[p1].U;
				element[p1].U = element[posR].L;
				element[posR].L = temp;
				temp = element[p1].L;
				element[p1].L = element[posR].U;
				element[posR].U = temp;
			}
			else
			{
				temp = element[p1].U;
				element[p1].U = element[posR].L;
				element[p1].L = temp;
			}
		}

	}
	else//������Ԫ�ƶ�����λ��
	{

		//ԭλ����Χ��ϵ������ɾ����������������
		int p2;
		if (p1 != posR)//ɾ����������
		{
			for (p2 = RHEAD[element[p1].i].first; element[element[p2].right].j < element[p1].j; p2 = element[p2].right);
			if (p1 == p2)
				RHEAD[element[p2].i].first = element[p2].right;
			else
				element[p2].right = element[p1].right;
		}
		if (p1 != posC)//ɾ����������
		{
			for (p2 = CHEAD[element[p1].j].first; element[element[p2].down].i < element[p1].i; p2 = element[p2].down);
			if (p1 == p2)
				CHEAD[element[p2].j].first = element[p2].down;
			else
				element[p2].down = element[p1].down;
		}
		//
		//��������ֵ����
		if (ni < nj)
		{
			element[p1].i = ni;
			element[p1].j = nj;
		}
		else
		{
			element[p1].i = nj;
			element[p1].j = ni;
			temp = element[p1].U;
			element[p1].U = element[p1].L;
			element[p1].L = temp;
		}


		if (p1 != posR)//���������λ��
		{
			element[p1].right = element[posR].right;
			element[posR].right = p1;
		}

		if (p1 != posC)//����������λ��
		{
			if (element[CHEAD[element[p1].j].first].i > element[p1].i)
			{
				element[p1].down = CHEAD[element[p1].j].first;
				CHEAD[element[p1].j].first = p1;

			}
			else
			{
				element[p1].down = element[posC].down;
				element[posC].down = p1;
			}
		}
	}
}
//�ڵ��Ż��ں˺�����������ŷ����仯�ķ���Ԫ���ƶ�������λ�á�
void sparse::swap_inner(const int &p, const int &q)
{
	int p1;
	double temp;
	int ni, nj, posR, posC;
	int next;
	//���������j<i
	for (p1 = CHEAD[p].first; element[p1].down >= 0; p1 = next)
	{
	
		next = element[p1].down;
		//���ѽ���������
		if (RHEAD[element[p1].i].flag)
			continue;
		
		//
		if (q == element[p1].i)
			ni = p;
		else
			ni = element[p1].i;
		nj = q;
		

		swap_inner_update(p1, temp, ni, nj, posR, posC);


		//��־Ԫ�ر��ƶ�����ֹ�ظ�
		if (ni <= nj)
			RHEAD[ni].flag = true;
		else
			CHEAD[ni].flag = true;
	}
	//�����Խ�Ԫ
	p1 = RHEAD[p].first;
	if (element[p1].i == element[p1].j&&RHEAD[element[p1].i].flag == false)
	{
		temp = element[p1].U;
		element[p1].U = element[RHEAD[q].first].U;
		element[RHEAD[q].first].U = temp;
		RHEAD[q].flag = true;
	}
	//���������j>i
	for (p1 = element[p1].right; p1 >= 0; p1 = next)
	{
	
		next = element[p1].right;
		if (CHEAD[element[p1].j].flag)
			continue;
		ni = q;
		if (element[p1].j == q)
			nj = p;
		else
			nj = element[p1].j;

		swap_inner_update(p1, temp, ni, nj, posR, posC);

		if (ni <= nj)
			CHEAD[nj].flag = true;
		else
			RHEAD[nj].flag = true;

	}
}

//��Ž�������
void sparse::swap(const int &p, const int &q)
{
	//�ڵ��Ž���
	swap_inner(p, q);
	swap_inner(q, p);
	//������
	for (int i = 0; i < m; i++)
		if (RHEAD[i].flag)
			RHEAD[i].flag = false;
	for (int i = 0; i < n; i++)
		if (CHEAD[i].flag)
			CHEAD[i].flag = false;
	//�Ի����߼�����
	int temp;
	temp = RHEAD[p].count;
	RHEAD[p].count = RHEAD[q].count;
	RHEAD[q].count = temp;

	CHEAD[p].count = RHEAD[p].count;
	CHEAD[q].count = RHEAD[q].count;

	temp = order[p];
	order[p] = order[q];
	order[q] = temp;
}

//Ѱַ����
//���룺����Ԫ�ص����к�i��j��λ��q
//������ҵ�Ԫ�أ���������true��q����ֵΪ����element�е�λ�ã�����Ϊfalse
bool sparse::find(const int& i, const int &j, int &p, int& q)
{
	int pi;
	int ui, uj;
	if (i < j)
	{
		ui = i;
		uj = j;
	}
	else if (i > j)
	{
		ui = j;
		uj = i;
	}
	else
	{
		p = RHEAD[i].first;
		return true;
	}

	for (pi = RHEAD[ui].first; element[pi].right >= 0 && element[element[pi].right].j < uj; pi = element[pi].right);
	if (element[pi].right >= 0 && element[element[pi].right].j == uj)
	{
		p = element[pi].right;
		return true;
	}
	else
	{
		p = pi;
		for (pi = CHEAD[uj].first; element[pi].down >= 0 && element[element[pi].down].i < ui; pi = element[pi].down);
		q = pi;
		return false;
	}
}

//չʾ����
//���ϵ������
void sparse::show()
{
	int posR, posC;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (find(i, j, posR, posC))
			{
				if (i <= j)
					cout << element[posR].U << " ";
				else
					cout << element[posR].L << " ";
			}
			else
				cout << "0 ";
		}
		cout << endl;
	}
	cout << endl;
}

void sparse::clear()
{
	vector <HEAD>(m).swap(CHEAD);
	vector <HEAD>(n).swap(RHEAD);
	vector <sp_element>().swap(element);


}


//���Ƿֽ�
void sparse::T2decomposition()
{
	for (int i = 0; i < m; i++)
	{
		order[i] = i;
	}
	int min;//��������Сֵ
	int q;//��������С�Ľڵ���
	for (int p = 0; p < m - 1; p++)
	{
		
		if (T2)
		{
			min = m;
			//ѡȡ�������ٵĽڵ�
			for (int i = p; i < m; i++)
			{
				if (RHEAD[i].count < min)
				{
					min = RHEAD[i].count;
					q = i;
				}
			}

			//
			if (p != q)
				swap(p, q);//ͨ����Ž�����ʹ�ó��������ٵĽڵ�Ϊ����ȥ�ڵ�

		}
		
		for (int k = element[RHEAD[p].first].right; k >= 0; k = element[k].right)
		{
			//��񻯲���
			//���ڲ��öԳ�ʽ�洢��ʡȥ�����Ǵ���ռ䣬��ʱ����ԪRHEADָ����еĶԽ�Ԫ

			element[k].U /= element[RHEAD[p].first].U;
			
			for (int l = element[RHEAD[p].first].right; l >= 0; l = element[l].right)
			{
				update(element[l].j, element[k].j, -element[k].U*element[l].L);//��ȥ
				RHEAD[element[l].i].count--;//���³�����
			}
		}
		
		
		
	}

	for (int p = 0; p < m - 1; p++)
	{
		for (int k = element[RHEAD[p].first].right; k >= 0; k = element[k].right)
		{
			element[k].L /= element[RHEAD[p].first].U;
		}
	}
}

//���º���
//���룺ע��Ԫ�����к�i��j��ע��Ԫ��ֵval
//��ע��Ԫ��ӵ�ϡ������У����¾���
//true��ʾ��ע��Ԫ
bool sparse::update(int i, int j, const double &val)
{
	int ei;
	sp_element add;
	int ui, uj;
	if (i <= j)
	{
		ui = i;
		uj = j;
		add.U = val;
	}
	else
	{
		ui = j;
		uj = i;
		add.L = val;
	}

	add.i = ui;
	add.j = uj;

	if (RHEAD[ui].first < 0 || element[RHEAD[ui].first].j > uj)//��������Ԫ
	{
		add.right = RHEAD[ui].first;
		element.push_back(add);
		RHEAD[ui].first = element.size() - 1;
	}
	else
	{
		for (ei = RHEAD[ui].first; element[ei].right >= 0 && element[element[ei].right].j <= uj; ei = element[ei].right);
		if (element[ei].j == uj)//��λ��Ϊ����Ԫ
		{
			if (i <= j)
				element[ei].U += val;
			else
				element[ei].L += val;
			current_position = ei;
			return false;
		}
		add.right = element[ei].right;
		element.push_back(add);
		element[ei].right = element.size() - 1;
	}
	if (CHEAD[uj].first < 0 || element[CHEAD[uj].first].i > ui)
	{
		(element.back()).down = CHEAD[uj].first;
		CHEAD[uj].first = element.size() - 1;
	}
	else
	{
		for (ei = CHEAD[uj].first; element[ei].down >= 0 && element[element[ei].down].i < ui; ei = element[ei].down);
		(element.back()).down = element[ei].down;
		element[ei].down = element.size() - 1;
	}

	if (i != j)//ͳ�Ƴ�����
	{
		RHEAD[i].count++;
		CHEAD[i].count++;
		RHEAD[j].count++;
		CHEAD[j].count++;
	}
	return true;

}

//���Է��������
//Ax=b
void sparse::solve(vector<double> &equa_b)
{
	T2decomposition();
	//ǰ�� 
	//Lz=b
	for (int p = 1; p < n; p++)
	{
		for (int k = CHEAD[p].first; element[k].down >= 0; k = element[k].down)
		{
			equa_b[order[p]] -= equa_b[order[element[k].i]] * element[k].L;
		}
	}
	//Dy=z;
	for (int p = 0; p < n; p++)
	{
		equa_b[order[p]] /= element[RHEAD[p].first].U;
	}
	//�ش�
	//Ux=y;
	for (int p = m - 2; p >= 0; p--)
	{
		for (int k = element[RHEAD[p].first].right; k >= 0; k = element[k].right)
		{
			equa_b[order[p]] -= equa_b[order[element[k].j]] * element[k].U;
		}
	}
}

double sparse::operator() (const int &i, const int &j)
{
	int k;
	int ui, uj;
	if (i <= j)
	{
		ui = i;
		uj = j;
	}
	else
	{
		ui = j;
		uj = i;
	}
	for (k = RHEAD[ui].first; element[k].right >= 0 && element[element[k].right].j <= uj; k = element[k].right);
	if (element[k].j == uj)
	{
		if (i <= j)
			return element[k].U;
		else
			return element[k].L;
	}
	else
		return 0;
		
}

//ѡ������Ԫ
bool sparse::chooserow(const int &row)
{
	if (row < m)
	{
		current_position = CHEAD[row].first;
		return true;
	}
	else
		return false;
}

bool sparse::right(int &i, int &j, int &r)
{
	
	if (current_position < 0)
		return false;
	else
	{
		if (element[current_position].i != r)
		{
			i = r;
			j = element[current_position].i;
		}
		else
		{
			i = r;
			j = element[current_position].j;
		}
		if (element[current_position].i != r)
			current_position = element[current_position].down;
		else
			current_position = element[current_position].right;
		return true;
	}
}

//���캯��
sparse::sparse(const int &mm, const int &nn, const matrix &list, const bool sym) :RHEAD(mm), CHEAD(nn), m(mm), n(nn), symmetric(sym), order(mm)
{
	T2 = true;
	int pi, pj;
	double pval;
	if (symmetric)
		for (int i = 0; i < list.length(0); i++)
		{
			pi = list(i, 0);
			pj = list(i, 1);
			pval = 1.0/list(i, 2);
			if (pi - 1 < m)
			{
				update(pi - 1, pi - 1, pval);
				if (pj - 1 < m)
				{
					update(pi - 1, pj - 1, -pval);
					update(pj - 1, pi - 1, -pval);
					update(pj - 1, pj - 1, pval);
				}
			}
			else if (pj - 1 < m)
				update(pj - 1, pj - 1, pval);
	
		}
	else
		for (int i = 0; i < list.length(0); i++)
		{
			pi = list(i, 0);
			pj = list(i, 1);
			pval = list(i, 2);
			update(pi - 1, pj - 1, pval);
		}
}