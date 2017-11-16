#include <fstream>
#include <iostream>
#include "mymatrix.h"
#include "importion.h"
using namespace std;
matrix txt2mat(const string &addr)
{
	int m;
	int n;
	ifstream fin(addr);
	if (fin.is_open() == false)
	{
		cout << "Error! Cannot open the file!";
		matrix c(1, 1);
		return c;
	}
	while (1)
	{
		if (fin.get() == 'm')
		{
			fin >> m;
			continue;
		}
		if (fin.get() == 'n')
		{
			fin >> n;
			break;
		}
	}
	matrix c(m, n);
	float val;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			fin >> val;
			c.write(i, j, val);
		}
	fin.close();
	return c;

}