#pragma once
#include <string>
#include "importion.h"
#include "mymatrix.h"
#include "sparse.h"
class NRPowerFlow
{
	friend class sparse;
	int N;
	int genN;//PV节点数
	sparse G;
	sparse B;
	sparse Jac;
	vector<double> detPQ;
	vector<double> PQs;
	vector<double> ef;//电压恒分量纵分量
	
	

public:

	bool update()
	{
		
		//a=(G*e) b=(B*f) c=(G*f) d=(B*e) 
		double ai, bi, ci, di;
		int i, j;
		double gg, bb, g, b;
		double bii, cii;

		for (int r = 0; r < N - 1; r++)
		{
			Jac.update(2 * r, 2 * r, 0);
			Jac.update(2 * r + 1, 2 * r + 1, 0);
			G.chooserow(r);
			B.chooserow(r);
			//对角
			gg = G(r, r);
			bb = B(r, r);
			ai = 0;

			bi = 0;

			ci = 0;

			di = 0;

			
			while (G.right(i, j, r))
			{
				g = G(i, j);
				b = B(i, j);
				ai += g*ef[2 * j] - b*ef[2 * j + 1];//G*E-B*f
				bii = g*ef[2 * i] + b*ef[2 * i + 1];//G*e+B*f
				bi += g*ef[2 * j] + b*ef[2 * j + 1];
				cii = g*ef[2 * i + 1] - b*ef[2 * i];//G*f-B*e
				ci += g*ef[2 * j + 1] - b*ef[2 * j];
				di += g*ef[2 * j + 1] + b*ef[2 * j];//G*f+B*e
				if (i != j && j < N - 1)
				{
					//非对角
					// dP/de
					Jac.update(2 * i + 1, 2 * j, -bii);
					// dP/df
					Jac.update(2 * i + 1, 2 * j + 1, -cii);
					if (r < N - 1 - genN)
					{
						//dQ/de
						Jac.update(2 * i, 2 * j, -cii);
						//dQ/df
						Jac.update(2 * i, 2 * j + 1, bii);
					}
				}
					


			}
			
			//对角
			if (r < N - 1 - genN)
			{
				// dP/de
				Jac.update(2 * r + 1, 2 * r, -ai - gg*ef[2 * r] - bb*ef[2 * r + 1]);
				// dP/df
				Jac.update(2 * r + 1, 2 * r + 1, -di + bb*ef[2 * r] - gg*ef[2 * r + 1]);
				//dQ/de
				Jac.update(2 * r, 2 * r, di + bb*ef[2 * r] - gg*ef[2 * r + 1]);
				//dQ/df
				Jac.update(2 * r, 2 * r + 1, -ai + gg*ef[2 * r] + bb*ef[2 * r + 1]);

				//detPQ
				detPQ[2 * r + 1] = PQs[2 * r + 1] - ef[2 * r] * ai - ef[2 * r + 1] * di;
				detPQ[2 * r] = PQs[2 * r] - ef[2 * r + 1] * ai + ef[2 * r] * di;
		

			}
			else
			{
				// dP/de
				Jac.update(2 * r + 1, 2 * r, -ai - gg*ef[2 * r] - bb*ef[2 * r + 1]);
				// dP/df
				Jac.update(2 * r + 1, 2 * r + 1, -di + bb*ef[2 * r] - gg*ef[2 * r + 1]);
				//dV2/de
				Jac.update(2 * r, 2 * r, -2 * ef[2 * r]);
				//dV2/df
				Jac.update(2 * r, 2 * r + 1, -2 * ef[2 * r + 1]);

				//detPV
				detPQ[2 * r + 1] = PQs[2 * r + 1] - ef[2 * r] * ai - ef[2 * r + 1] * di;
				detPQ[2 * r] = PQs[2 * r] - ef[2 * r] * ef[2 * r] - ef[2 * r + 1] * ef[2 * r + 1];
			}


			
		}
		//Jac.show();
		Jac.solve(detPQ);
		Jac.clear();
	
		
		double val = 0;

		
		for (int k = 0; k < detPQ.size(); k++)
		{
			ef[k] -= detPQ[k];
			if (fabs(detPQ[k]) > val)
				val = fabs(detPQ[k]);
		}
		

		if (val > 1e-6)
		{
			cout << "err=" << val << endl;
			return true;
		}
		else
		{
			cout << endl;
			for (int k = 0; k < ef.size(); k++)
			{
				cout << ef[k] << endl;
			}
			return false;
		}
			
	
	}

	NRPowerFlow(int nn, string net, string power, string transformer) :N(nn), G(nn, nn), B(nn, nn), Jac(2 * (nn - 1), 2 * (nn - 1)), detPQ(2 * (nn - 1), 0), ef(2 * nn , 0), PQs(2 * (nn - 1), 0)
	{
		genN = 0;
		int from, to;
		double rx2, g, b, b0;
		vector<int> order(N, 0);
		for (int i = 0; i < N; i++)
		{
			order[i] = i;
		}
		for (int i = 0; i < N; i++)
		{
			
			ef[2 * i] = 1;
		}
		matrix data = txt2mat(net);
		matrix ppower = txt2mat(power);
		matrix trans = txt2mat(transformer);
		data.show();
		ppower.show();
		int temp = 0;
		for (int i = 0; i < N; i++)
		{
			if (ppower(i, 6) == 1)//genbus
			{
				temp=order[i];
				order[i] = order[N - 1 - genN - 1];
				order[N - 1 - genN - 1] = temp;
				ef[2 * order[i]] = ppower(i, 0);
				genN++;
			}
			else if (ppower(i, 6) == 2)
			{
				temp = order[i];
				order[i] = order[N - 1];
				order[N - 1] = temp;
				ef[2 * order[i]] = ppower(i, 0);
			}

		}

		int k;
		for (int i = 0; i < N; i++)
		{
			k = order[i];
			if (k < N - 1 - genN)
			{
				PQs[2 * k + 1] = ppower(i, 2) - ppower(i, 4);
				PQs[2 * k] = ppower(i, 3) - ppower(i, 5);
			}
			else if (k < N - 1)
			{
				PQs[2 * k + 1] = ppower(i, 2) - ppower(i, 4);
				PQs[2 * k] = ppower(i, 0)*ppower(i, 0);
				
			}
			

		}

		for (int i = 0; i < data.length(0); i++)
		{
			from = order[data(i, 0) - 1];
			to = order[data(i, 1) - 1];
			rx2 = data(i, 2)*data(i, 2) + data(i, 3)*data(i, 3);
			g = data(i, 2) / rx2;
			b = -data(i, 3) / rx2;
			b0 = data(i, 4);
			G.update(from, from, g);
			B.update(from, from, b + b0);
			G.update(from, to, -g);
			G.update(to, from, -g);
			B.update(from, to, -b);
			B.update(to, from, -b);
			G.update(to, to, g);
			B.update(to, to, b + b0);

		}
		
		double tap;
		for (int i = 0; i < trans.length(0); i++)
		{
			from = order[trans(i, 0) - 1];
			to = order[trans(i, 1) - 1];
			tap = trans(i, 2);
		
			g = G(from, to);
			b = B(from, to);
			G.update(from, to, g / tap - g);
			B.update(from, to, b / tap - b);
			G.update(to, from, g / tap - g);
			B.update(to, from, b / tap - b);
			G.update(to, to, g - g / (tap*tap));
			B.update(to, to, b - b / (tap*tap));

		}
			
		
		B.show();
		G.show();

	}

};