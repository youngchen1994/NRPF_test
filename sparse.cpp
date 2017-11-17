#include "mymatrix.h"
using namespace std;
#include <vector>
#include <iomanip>
#include "sparse.h"

//节点优化内核函数，函数功能为将发生变化非零元移动到其新位置。
void sparse::swap_inner_update(int &p1, double &temp, int &ni, int &nj, int &posR, int &posC)
{
	if (find(ni, nj, posR, posC))//找到对换元
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
		else//上下三角值交换
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
	else//将非零元移动到空位置
	{

		//原位置周围关系修正，删除与其相连的链表
		int p2;
		if (p1 != posR)//删除横向链接
		{
			for (p2 = RHEAD[element[p1].i].first; element[element[p2].right].j < element[p1].j; p2 = element[p2].right);
			if (p1 == p2)
				RHEAD[element[p2].i].first = element[p2].right;
			else
				element[p2].right = element[p1].right;
		}
		if (p1 != posC)//删除纵向链接
		{
			for (p2 = CHEAD[element[p1].j].first; element[element[p2].down].i < element[p1].i; p2 = element[p2].down);
			if (p1 == p2)
				CHEAD[element[p2].j].first = element[p2].down;
			else
				element[p2].down = element[p1].down;
		}
		//
		//上下三角值交换
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


		if (p1 != posR)//插入横向新位置
		{
			element[p1].right = element[posR].right;
			element[posR].right = p1;
		}

		if (p1 != posC)//插入纵向新位置
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
//节点优化内核函数，遍历编号发生变化的非零元，移动到其新位置。
void sparse::swap_inner(const int &p, const int &q)
{
	int p1;
	double temp;
	int ni, nj, posR, posC;
	int next;
	//纵向遍历，j<i
	for (p1 = CHEAD[p].first; element[p1].down >= 0; p1 = next)
	{
	
		next = element[p1].down;
		//若已交换，跳过
		if (RHEAD[element[p1].i].flag)
			continue;
		
		//
		if (q == element[p1].i)
			ni = p;
		else
			ni = element[p1].i;
		nj = q;
		

		swap_inner_update(p1, temp, ni, nj, posR, posC);


		//标志元素被移动，防止重复
		if (ni <= nj)
			RHEAD[ni].flag = true;
		else
			CHEAD[ni].flag = true;
	}
	//交换对角元
	p1 = RHEAD[p].first;
	if (element[p1].i == element[p1].j&&RHEAD[element[p1].i].flag == false)
	{
		temp = element[p1].U;
		element[p1].U = element[RHEAD[q].first].U;
		element[RHEAD[q].first].U = temp;
		RHEAD[q].flag = true;
	}
	//横向遍历，j>i
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

//编号交换函数
void sparse::swap(const int &p, const int &q)
{
	//节点编号交换
	swap_inner(p, q);
	swap_inner(q, p);
	//清除标记
	for (int i = 0; i < m; i++)
		if (RHEAD[i].flag)
			RHEAD[i].flag = false;
	for (int i = 0; i < n; i++)
		if (CHEAD[i].flag)
			CHEAD[i].flag = false;
	//对换出线计数器
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

//寻址函数
//输入：待查元素的行列号i，j，位置q
//输出：找到元素，函数返回true，q被赋值为其在element中的位置；否则为false
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

//展示函数
//输出系数矩阵
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
	int posR, posC;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (find(i, j, posR, posC))
			{
				if (i <= j)
					element[posR].U = 0;
				else
					element[posR].L = 0;
			}
		}
	}
}


//三角分解
void sparse::T2decomposition()
{
	for (int i = 0; i < m; i++)
	{
		order[i] = i;
	}
	int min;//出线数最小值
	int q;//出线数最小的节点编号
	for (int p = 0; p < m - 1; p++)
	{
		
		
		min = m;
		//选取出线最少的节点
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
			swap(p, q);//通过编号交换，使得出线数最少的节点为待消去节点
			
			
		
		for (int k = element[RHEAD[p].first].right; k >= 0; k = element[k].right)
		{
			//规格化操作
			//由于采用对称式存储，省去下三角储存空间，此时行首元RHEAD指向该行的对角元

			element[k].U /= element[RHEAD[p].first].U;
			
			for (int l = element[RHEAD[p].first].right; l >= 0; l = element[l].right)
			{
				update(element[l].j, element[k].j, -element[k].U*element[l].L);//消去
				RHEAD[element[l].i].count--;//更新出线数
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

//更新函数
//输入：注入元的行列号i和j，注入元的值val
//将注入元添加到稀疏矩阵中，更新矩阵
//true表示有注入元
bool sparse::update(const int &i, const int &j, const double &val)
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

	if (RHEAD[ui].first < 0 || element[RHEAD[ui].first].j > uj)//插入行首元
	{
		add.right = RHEAD[ui].first;
		element.push_back(add);
		RHEAD[ui].first = element.size() - 1;
	}
	else
	{
		for (ei = RHEAD[ui].first; element[ei].right >= 0 && element[element[ei].right].j <= uj; ei = element[ei].right);
		if (element[ei].j == uj)//该位置为非零元
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
	if (i != j)//统计出线数
	{
		RHEAD[i].count++;
		CHEAD[i].count++;
		RHEAD[j].count++;
		CHEAD[j].count++;
	}
	return true;

}

//线性方程组求解
//Ax=b
void sparse::solve(vector<double> &equa_b)
{
	T2decomposition();
	//前代 
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
	//回代
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

//选中行首元
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

//构造函数
sparse::sparse(const int &mm, const int &nn, const matrix &list, const bool sym) :RHEAD(mm), CHEAD(nn), m(mm), n(nn), symmetric(sym), order(mm)
{
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