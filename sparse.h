#pragma once
#include "mymatrix.h"
using namespace std;
#include <vector>
#include <iomanip>
#include <iostream>
#include <typeinfo>

//稀疏矩阵——十字链表格式

//非零元结构体
struct sp_element{
	double U = 0, L = 0;//值
	int i;//行号
	int j;//列号
	int right = -1;//该元素右侧第一个非零元的位置
	int down = -1;//该元素下侧第一个非零元的位置
};

//行头/列头（入口）结构体
struct HEAD{
	int count = 0;//非零元个数
	int first = -1;//行列首元地址
	bool flag;//选中标志
};//行、列开头格式

//稀疏矩阵类
class sparse
{
	bool symmetric;
	int m;//行数
	int n;//列数
	vector<sp_element> element;//非零元集合
	vector<HEAD> RHEAD, CHEAD;//行列入口
	vector<int> order;
	int current_position;//当前指向非零元的位置


	//节点优化内核函数，遍历编号发生变化的非零元，移动到其新位置。
	void swap_inner(const int &p, const int &q);

	//编号交换函数
	void swap(const int &p, const int &q);

	//寻址函数
	//输入：待查元素的行列号i，j，位置q
	//输出：找到元素，函数返回true，q被赋值为其在element中的位置；否则为false
	bool find(const int& i, const int &j, int &p, int& q);

public: 
	//读取矩阵元A（i，j）
	double operator() (const int &i, const int &j);

	bool chooserow(const int &row);//选中行首元
	bool right(int &i, int &j);//搜索右侧下一个非零元

	//展示函数
	//输出系数矩阵
	void show();

	//三角分解
	//使用Tinney-2方法进行节点优化编号
	//半动态优先消去出线数最少的节点
	//将目标矩阵分解成L、D、U三个稀疏矩阵
	void T2decomposition();

	//线性方程组求解
	//Ax=b
	void solve(vector<double> &b);

	//更新函数
	//输入：注入元的行列号i和j，注入元的值val
	//将注入元添加到稀疏矩阵中，更新矩阵
	void update(const int &i, const int &j, const double &val);

	//构造函数
	//方法1：已知矩阵维度，m*n，构造空的稀疏矩阵
	//方法2：已知矩阵维度m*n，非零元行列号及数值的列表，构造稀疏矩阵
	sparse(const int &mm, const int &nn) :RHEAD(mm), CHEAD(nn), m(mm), n(nn), order(mm) {};
	sparse(const int &mm, const int &nn, const matrix &list, const bool sym);
};

