#include <vector>
///////////////////////////////
#define DD 23
#define DD3 13
///////////////////////////////
double h2[] = { -0.002, -0.003, 0.006, 0.006, -0.013,
-0.012, 0.030, 0.023, -0.078, -0.035,
0.307, 0.542, 0.307, -0.035, -0.078,
0.023, 0.030, -0.012, -0.013, 0.006,
0.006, -0.003, -0.002 };  //Haar小波低通滤波系数

double g2[] = { 0.002, -0.003, -0.006, 0.006, 0.013, -0.012,
-0.030, 0.023, 0.078, -0.035, -0.307, 0.542,
-0.307, -0.035, 0.078, 0.023, -0.030, -0.012,
0.013, 0.006, -0.006, -0.003, 0.002 };   //Haar小波高通滤波系数
////////////////////////////
double h[] = { 0.3326705529500825, 0.8068915093110924, 0.4598775021184914,
-0.1350110200102546, -0.0854412738820267, 0.0352262918857095 };//Daubechies小波低通滤波系数
double g[] = { 0.0352262918857095, 0.0854412738820267, -0.1350110200102546,
-0.4598775021184914, 0.8068915093110924, -0.3326705529500825 };   //Daubechies小波高通滤波系数
////////////////////////////////

double h3[DD3] = { -0.00332761, 0.00569794, 0.0196637, -0.0482603, -0.0485391,
0.292562, 0.564406, 0.292562, -0.0485391, -0.0482602, -0.0196637, 0.00569794, -0.0033276 };
//Morlet小波低通滤波系数
double g3[DD3] = { 0.00332761, 0.00569794, -0.0196637, -0.0482603, 0.0485391,
0.292562, -0.564406, 0.292562, 0.0485391, -0.0482602, 0.0196637, 0.00569794, 0.0033276 };
//Morlet小波高通滤波系数
///////////////////////////////////////////
//小波长度(6),逼近(2倍原始长度),细节,尺度数组大小,尺度数组
void DWT(int wlen, double c[], double d[], int m, int sca[])  //Daubechies小波变换
{
	int i, j, k, mid, flag[20];
	double p, q;

	for (flag[0] = 0, i = 0; i < m; i++)
	{
		flag[i + 1] = flag[i] + sca[i];
	}

	for (j = 1; j <= m; j++)
	{
		for (i = 0; i < sca[j]; i++)
		{
			p = 0;
			q = 0;

			for (k = 0; k < wlen; k++)
			{
				mid = k + 2 * i;
				if (mid >= sca[j - 1])  mid = mid - sca[j - 1];

				p = p + h[k] * c[flag[j - 1] + mid];
				q = q + g[k] * c[flag[j - 1] + mid];
			}

			c[flag[j] + i] = p;
			d[flag[j] + i] = q;
		}
	}

}

void IDWT(double g[], double h[], int wlen, double c[], double d[], int m, int sca[])
{
	int i, j, k, mid, flag[20];
	double p, q;

	for (flag[0] = 0, i = 0; i < m; i++)
	{
		flag[i + 1] = flag[i] + sca[i];
	}

	for (k = m; k > 0; k--)
	{
		for (i = 0; i < sca[k]; i++)
		{
			p = 0;
			q = 0;

			for (j = 0; j < wlen / 2; j++)
			{
				mid = i - j;
				if (mid < 0)  mid = sca[k] + (i - j);

				p += h[2 * j] * c[flag[k] + mid] + g[2 * j] * d[flag[k] + mid];
				q += h[2 * j + 1] * c[flag[k] + mid] + g[2 * j + 1] * d[flag[k] + mid];
			}

			c[flag[k - 1] + 2 * i] = p;
			c[flag[k - 1] + 2 * i + 1] = q;
		}
	}

}
////////////////////////////
void Radom05(int n, double random[], double mean, double fangcha)
{

	int i, j, m;
	double s, w, v, r, t;
	s = 65536.0;
	w = 2053.0;
	v = 13849.0;
	r = 0.0;

	for (i = 0; i <= n - 1; i++)
	{
		t = 0.0;
		for (j = 0; j <= 11; j++)
		{
			r = w*r + v;
			m = r / s; r = r - m*s;
			t = t + r / s;
		}

		random[i] = mean + fangcha*(t - 6.0);
	}

	return;
}
/*
////////////////////////////////////////////////////////
int edge_function(int x, int xsize)
{
	if (x < 0)
		x = -x;
	if (x >= xsize)
		x = xsize * 2 - x - 2;
	return(x);
}
void DWT1(double c[], double d[], int m, int sca[])  //Haar小波变换
{
	int i, j, k, mid, flag[20];
	double p, q;

	for (flag[0] = 0, i = 0; i < m; i++)
	{
		flag[i + 1] = flag[i] + sca[i];
	}

	for (k = 1; k <= m; k++)
	{
		for (i = 0; i < sca[k - 1]; i += 2)
		{
			p = 0;
			q = 0;
			for (j = -(DD - 1) / 2; j <= (DD - 1) / 2; j++)
				p = p + h2[j + (DD - 1) / 2] * c[flag[k - 1] + edge_function(i + j, sca[k - 1])];
			for (j = -(DD - 1) / 2 + 1; j <= (DD - 1) / 2 + 1; j++)
				q = q + g2[j + (DD - 1) / 2 - 1] * c[flag[k - 1] + edge_function(i + j, sca[k - 1])];

			c[flag[k] + i / 2] = p;
			d[flag[k] + i / 2] = q;
		}

	}

}

void IDWT1(double c[], double d[], int m, int sca[])
{
	int i, j, k, n, mid, flag[20];
	double p, q;
	double temp1, temp2;
	vector<double> ctemp[SHUZU], dtemp[SHUZU];

	for (flag[0] = 0, i = 0; i < m; i++)
	{
		flag[i + 1] = flag[i] + sca[i];
	}

	for (k = m; k > 0; k--)
	{
		for (n = 0; n < sca[k]; n++)
		{
			ctemp[flag[k - 1] + n * 2] = c[flag[k] + n];
			ctemp[flag[k - 1] + n * 2 + 1] = 0;

			dtemp[flag[k - 1] + n * 2] = d[flag[k] + n];
			dtemp[flag[k - 1] + n * 2 + 1] = 0;
		}
		for (i = 0; i < sca[k - 1]; i++)
		{
			temp1 = 0;
			for (j = -(DD - 1) / 2; j <= (DD - 1) / 2; j++)
				temp1 = temp1 + h2[j + (DD - 1) / 2] * ctemp[flag[k - 1] + edge_function(i + j, sca[k - 1])];
			c[flag[k - 1] + i] = temp1;
		}
		for (i = 0; i < sca[k - 1]; i++)
		{
			temp1 = 0;
			for (j = -(DD - 1) / 2 - 1; j <= (DD - 1) / 2 - 1; j++)
				temp1 = temp1 + g2[j + (DD - 1) / 2 + 1] * dtemp[flag[k - 1] + edge_function(i + j, sca[k - 1])];
			temp2 = temp1 + c[flag[k - 1] + i];
			c[flag[k - 1] + i] = temp2;

			c[flag[k - 1] + i] = 2 * c[flag[k - 1] + i];
		}

	}

}
///////////////////////////////////
void DWT3(double c[], double d[], int m, int sca[])   //Morlet小波变换
{
	int i, j, k, mid, flag[20];
	double p, q;

	for (flag[0] = 0, i = 0; i < m; i++)
	{
		flag[i + 1] = flag[i] + sca[i];
	}

	for (k = 1; k <= m; k++)
	{
		for (i = 0; i < sca[k - 1]; i += 2)
		{
			p = 0;
			q = 0;
			for (j = -(DD3 - 1) / 2; j <= (DD3 - 1) / 2; j++)
				p = p + h3[j + (DD3 - 1) / 2] * c[flag[k - 1] + edge_function(i + j, sca[k - 1])];
			for (j = -(DD3 - 1) / 2 + 1; j <= (DD3 - 1) / 2 + 1; j++)
				q = q + g3[j + (DD3 - 1) / 2 - 1] * c[flag[k - 1] + edge_function(i + j, sca[k - 1])];

			c[flag[k] + i / 2] = p;
			d[flag[k] + i / 2] = q;
		}

	}

}

void IDWT3(double c[], double d[], int m, int sca[])
{
	int i, j, k, n, mid, flag[20];
	double p, q;
	double temp1, temp2;
	double ctemp[SHUZU], dtemp[SHUZU];

	for (flag[0] = 0, i = 0; i < m; i++)
	{
		flag[i + 1] = flag[i] + sca[i];
	}

	for (k = m; k > 0; k--)
	{
		for (n = 0; n < sca[k]; n++)
		{
			ctemp[flag[k - 1] + n * 2] = c[flag[k] + n];
			ctemp[flag[k - 1] + n * 2 + 1] = 0;

			dtemp[flag[k - 1] + n * 2] = d[flag[k] + n];
			dtemp[flag[k - 1] + n * 2 + 1] = 0;
		}
		for (i = 0; i < sca[k - 1]; i++)
		{
			temp1 = 0;
			for (j = -(DD3 - 1) / 2; j <= (DD3 - 1) / 2; j++)
				temp1 = temp1 + h3[j + (DD3 - 1) / 2] * ctemp[flag[k - 1] + edge_function(i + j, sca[k - 1])];
			c[flag[k - 1] + i] = temp1;
		}
		for (i = 0; i < sca[k - 1]; i++)
		{
			temp1 = 0;
			for (j = -(DD3 - 1) / 2 - 1; j <= (DD3 - 1) / 2 - 1; j++)
				temp1 = temp1 + g3[j + (DD3 - 1) / 2 + 1] * dtemp[flag[k - 1] + edge_function(i + j, sca[k - 1])];
			temp2 = temp1 + c[flag[k - 1] + i];
			c[flag[k - 1] + i] = temp2;

			c[flag[k - 1] + i] = 2 * c[flag[k - 1] + i];
		}
	}
}*/