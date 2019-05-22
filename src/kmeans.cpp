# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

#include "kmeans.h"


void kmns(double a[], int m, int n, double c[], int k, int ic1[], int nc[],
	int iter, double wss[], int *ifault)

{
	double aa;
	double *an1;
	double *an2;
	double *d;
	double da;
	double db;
	double dc;
	double dt[2];
	int i;
	int *ic2;
	int ii;
	int ij;
	int il;
	int indx;
	int *itran;
	int j;
	int l;
	int *live;
	int *ncp;
	double temp;

	*ifault = 0;

	if (k <= 1 || m <= k)
	{
		*ifault = 3;
		return;
	}
	ic2 = new int[m];
	an1 = new double[k];
	an2 = new double[k];
	ncp = new int[k];
	d = new double[m];
	itran = new int[k];
	live = new int[k];
	for (i = 1; i <= m; i++)
	{
		ic1[i - 1] = 1;
		ic2[i - 1] = 2;

		for (il = 1; il <= 2; il++)
		{
			dt[il - 1] = 0.0;
			for (j = 1; j <= n; j++)
			{
				da = a[i - 1 + (j - 1)*m] - c[il - 1 + (j - 1)*k];
				dt[il - 1] = dt[il - 1] + da * da;
			}
		}

		if (dt[1] < dt[0])
		{
			ic1[i - 1] = 2;
			ic2[i - 1] = 1;
			temp = dt[0];
			dt[0] = dt[1];
			dt[1] = temp;
		}

		for (l = 3; l <= k; l++)
		{
			db = 0.0;
			for (j = 1; j <= n; j++)
			{
				dc = a[i - 1 + (j - 1)*m] - c[l - 1 + (j - 1)*k];
				db = db + dc * dc;
			}

			if (db < dt[1])
			{
				if (dt[0] <= db)
				{
					dt[1] = db;
					ic2[i - 1] = l;
				}
				else
				{
					dt[1] = dt[0];
					ic2[i - 1] = ic1[i - 1];
					dt[0] = db;
					ic1[i - 1] = l;
				}
			}
		}
	}

	for (l = 1; l <= k; l++)
	{
		nc[l - 1] = 0;
		for (j = 1; j <= n; j++)
		{
			c[l - 1 + (j - 1)*k] = 0.0;
		}
	}

	for (i = 1; i <= m; i++)
	{
		l = ic1[i - 1];
		nc[l - 1] = nc[l - 1] + 1;
		for (j = 1; j <= n; j++)
		{
			c[l - 1 + (j - 1)*k] = c[l - 1 + (j - 1)*k] + a[i - 1 + (j - 1)*m];
		}
	}
	*ifault = 1;

	for (l = 1; l <= k; l++)
	{
		if (nc[l - 1] == 0)
		{
			*ifault = 1;
			return;
		}

	}

	*ifault = 0;

	for (l = 1; l <= k; l++)
	{
		aa = (double)(nc[l - 1]);

		for (j = 1; j <= n; j++)
		{
			c[l - 1 + (j - 1)*k] = c[l - 1 + (j - 1)*k] / aa;
		}
		an2[l - 1] = aa / (aa + 1.0);

		if (1.0 < aa)
		{
			an1[l - 1] = aa / (aa - 1.0);
		}
		else
		{
			an1[l - 1] = r8_huge();
		}
		itran[l - 1] = 1;
		ncp[l - 1] = -1;
	}

	indx = 0;
	*ifault = 2;

	for (ij = 1; ij <= iter; ij++)
	{
		optra(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, &indx);
		if (indx == m)
		{
			*ifault = 0;
			break;
		}
		qtran(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &indx);
		if (k == 2)
		{
			*ifault = 0;
			break;
		}
		for (l = 1; l <= k; l++)
		{
			ncp[l - 1] = 0;
		}

	}
	if (*ifault == 2)
	{
		cout << "\n";
		cout << "KMNS - Warning!\n";
		cout << "  Maximum number of iterations reached\n";
		cout << "  without convergence.\n";
	}
	for (l = 1; l <= k; l++)
	{
		wss[l - 1] = 0.0;
		for (j = 1; j <= n; j++)
		{
			c[l - 1 + (j - 1)*k] = 0.0;
		}
	}

	for (i = 1; i <= m; i++)
	{
		ii = ic1[i - 1];
		for (j = 1; j <= n; j++)
		{
			c[ii - 1 + (j - 1)*k] = c[ii - 1 + (j - 1)*k] + a[i - 1 + (j - 1)*m];
		}
	}

	for (j = 1; j <= n; j++)
	{
		for (l = 1; l <= k; l++)
		{
			c[l - 1 + (j - 1)*k] = c[l - 1 + (j - 1)*k] / (double)(nc[l - 1]);
		}
		for (i = 1; i <= m; i++)
		{
			ii = ic1[i - 1];
			da = a[i - 1 + (j - 1)*m] - c[ii - 1 + (j - 1)*k];
			wss[ii - 1] = wss[ii - 1] + da * da;
		}
	}

	delete[] ic2;
	delete[] an1;
	delete[] an2;
	delete[] ncp;
	delete[] d;
	delete[] itran;
	delete[] live;

	return;
}


void optra(double a[], int m, int n, double c[], int k, int ic1[],
	int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[],
	int itran[], int live[], int *indx){
	double al1;
	double al2;
	double alt;
	double alw;
	double da;
	double db;
	double dc;
	double dd;
	double de;
	double df;
	int i;
	int j;
	int l;
	int l1;
	int l2;
	int ll;
	double r2;
	double rr;
	for (l = 1; l <= k; l++)
	{
		if (itran[l - 1] == 1)
		{
			live[l - 1] = m + 1;
		}
	}

	for (i = 1; i <= m; i++)
	{
		*indx = *indx + 1;
		l1 = ic1[i - 1];
		l2 = ic2[i - 1];
		ll = l2;
		if (1 < nc[l1 - 1])
		{
			if (ncp[l1 - 1] != 0)
			{
				de = 0.0;
				for (j = 1; j <= n; j++)
				{
					df = a[i - 1 + (j - 1)*m] - c[l1 - 1 + (j - 1)*k];
					de = de + df * df;
				}
				d[i - 1] = de * an1[l1 - 1];
			}
			da = 0.0;
			for (j = 1; j <= n; j++)
			{
				db = a[i - 1 + (j - 1)*m] - c[l2 - 1 + (j - 1)*k];
				da = da + db * db;
			}
			r2 = da * an2[l2 - 1];

			for (l = 1; l <= k; l++)
			{
				if ((i < live[l1 - 1] || i < live[l2 - 1]) && l != l1 && l != ll)
				{
					rr = r2 / an2[l - 1];

					dc = 0.0;
					for (j = 1; j <= n; j++)
					{
						dd = a[i - 1 + (j - 1)*m] - c[l - 1 + (j - 1)*k];
						dc = dc + dd * dd;
					}

					if (dc < rr)
					{
						r2 = dc * an2[l - 1];
						l2 = l;
					}
				}
			}
			if (d[i - 1] <= r2)
			{
				ic2[i - 1] = l2;
			}
			else
			{
				*indx = 0;
				live[l1 - 1] = m + i;
				live[l2 - 1] = m + i;
				ncp[l1 - 1] = i;
				ncp[l2 - 1] = i;
				al1 = (double)(nc[l1 - 1]);
				alw = al1 - 1.0;
				al2 = (double)(nc[l2 - 1]);
				alt = al2 + 1.0;
				for (j = 1; j <= n; j++)
				{
					c[l1 - 1 + (j - 1)*k] = (c[l1 - 1 + (j - 1)*k] * al1 - a[i - 1 + (j - 1)*m]) / alw;
					c[l2 - 1 + (j - 1)*k] = (c[l2 - 1 + (j - 1)*k] * al2 + a[i - 1 + (j - 1)*m]) / alt;
				}
				nc[l1 - 1] = nc[l1 - 1] - 1;
				nc[l2 - 1] = nc[l2 - 1] + 1;
				an2[l1 - 1] = alw / al1;
				if (1.0 < alw)
				{
					an1[l1 - 1] = alw / (alw - 1.0);
				}
				else
				{
					an1[l1 - 1] = r8_huge();
				}
				an1[l2 - 1] = alt / al2;
				an2[l2 - 1] = alt / (alt + 1.0);
				ic1[i - 1] = l2;
				ic2[i - 1] = l1;
			}
		}

		if (*indx == m)
		{
			return;
		}
	}
	for (l = 1; l <= k; l++)
	{
		itran[l - 1] = 0;
		live[l - 1] = live[l - 1] - m;
	}

	return;
}


void qtran(double a[], int m, int n, double c[], int k, int ic1[],
	int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[],
	int itran[], int *indx){
	double al1;
	double al2;
	double alt;
	double alw;
	double da;
	double db;
	double dd;
	double de;
	int i;
	int icoun;
	int istep;
	int j;
	int l1;
	int l2;
	double r2;
	icoun = 0;
	istep = 0;

	for (; ; )
	{
		for (i = 1; i <= m; i++)
		{
			icoun = icoun + 1;
			istep = istep + 1;
			l1 = ic1[i - 1];
			l2 = ic2[i - 1];
			if (1 < nc[l1 - 1])
			{
				if (istep <= ncp[l1 - 1])
				{
					da = 0.0;
					for (j = 1; j <= n; j++)
					{
						db = a[i - 1 + (j - 1)*m] - c[l1 - 1 + (j - 1)*k];
						da = da + db * db;
					}
					d[i - 1] = da * an1[l1 - 1];
				}				
				if (istep < ncp[l1 - 1] || istep < ncp[l2 - 1])
				{
					r2 = d[i - 1] / an2[l2 - 1];

					dd = 0.0;
					for (j = 1; j <= n; j++)
					{
						de = a[i - 1 + (j - 1)*m] - c[l2 - 1 + (j - 1)*k];
						dd = dd + de * de;
					}					
					if (dd < r2)
					{
						icoun = 0;
						*indx = 0;
						itran[l1 - 1] = 1;
						itran[l2 - 1] = 1;
						ncp[l1 - 1] = istep + m;
						ncp[l2 - 1] = istep + m;
						al1 = (double)(nc[l1 - 1]);
						alw = al1 - 1.0;
						al2 = (double)(nc[l2 - 1]);
						alt = al2 + 1.0;
						for (j = 1; j <= n; j++)
						{
							c[l1 - 1 + (j - 1)*k] = (c[l1 - 1 + (j - 1)*k] * al1 - a[i - 1 + (j - 1)*m]) / alw;
							c[l2 - 1 + (j - 1)*k] = (c[l2 - 1 + (j - 1)*k] * al2 + a[i - 1 + (j - 1)*m]) / alt;
						}
						nc[l1 - 1] = nc[l1 - 1] - 1;
						nc[l2 - 1] = nc[l2 - 1] + 1;
						an2[l1 - 1] = alw / al1;
						if (1.0 < alw)
						{
							an1[l1 - 1] = alw / (alw - 1.0);
						}
						else
						{
							an1[l1 - 1] = r8_huge();
						}
						an1[l2 - 1] = alt / al2;
						an2[l2 - 1] = alt / (alt + 1.0);
						ic1[i - 1] = l2;
						ic2[i - 1] = l1;
					}
				}
			}
		
			if (icoun == m)
			{
				return;
			}
		}
	}
}

double r8_huge()
{
	double value;

	value = 1.0E+30;

	return value;
}