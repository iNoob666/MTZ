#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace std;


double D[4][4] = {
	   {4, 2, 2, 1},
	   {2, 4, 1, 2},
	   {2, 1, 4, 2},
	   {1, 2, 2, 4},
};

double G1[4][4] = {
   {2, 1, -2, -1},
   {1, 2, -1, -2},
   {-2, -1, 2, 1},
   {-1, -2, 1, 2},
};

double G2[4][4] = {
   {2, -2, 1, -1},
   {-2, 2, -1, 1},
   {1, -1, 2, -2},
   {-1, 1, -2, 2},
};

double G3[4][4] = {
   {-2, 2, -1, 1},
   {-1, 1, -2, 2},
   {2, -2, 1, -1},
   {1, -1, 2, -2},
};

double G3T[4][4] = {
   {-2, -1, 2, 1},
   {2, 1, -2, -1},
   {-1, -2, 1, 2},
   {1, 2, -1, -2},
};

vector<double> operator/(const vector<double> &a, const vector<double> &b) {
	double coefB = b[0] * b[0] + b[1] * b[1];
	vector<double> oppositeB = { b[0] / coefB, -(b[1] / coefB) };
	return { a[0] * oppositeB[0] - a[1] * oppositeB[1], a[0] * oppositeB[1] + a[1] * oppositeB[0] };
}

vector<double> operator*(const vector<double>& a, const vector<double>& b) {
	return { a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] };
}

vector<double> operator+(const vector<double>& a, const vector<double>& b) {
	return { a[0] + b[0], a[1] + b[1] };
}

vector<double> operator-(const vector<double>& a, const vector<double>& b) {
	return { a[0] - b[0], a[1] - b[1] };
}


struct matrix // Структура для матрицы, которая хранится в
	// разряженном строчно-столбцовом формате относительно блочных элементов (размер блока 2x2) 
{
	vector<vector<double>> di;// диагональные элементы
	vector<vector<double>> al;// элементы нижнего тругольника
	vector<int> li;// профиль матрицы
	vector<int> lj;// портрет матрицы
};


vector<vector<double>> formxyz(const vector<double >& x, const vector<double >& y, const vector<double >& z) {
	ofstream fout("xyz.txt");
	vector<vector<double>> tmp;
	int count = 0;
	for (double iz : z) {
		for (double iy : y) {
			for (double ix : x) {
				vector<double> tmpxyz;
				tmpxyz.emplace_back(ix);
				tmpxyz.emplace_back(iy);
				tmpxyz.emplace_back(iz);
				tmp.emplace_back(tmpxyz);
				fout << count++ << ' ' << ix << ' ' << iy << ' ' << iz << std::endl;
			}
		}
	}
	fout.close();
	return tmp;
}

vector<vector<pair<int, pair<int, int>>>> formnvtr(const vector<double>& x, const vector<double>& y, const vector<double>& z) {
	vector<vector<pair<int, pair<int, int>>>> tmpEdges;
	vector<pair<int, pair<int, int>>> aEdges;

	ofstream fout("nvtr.txt");
	int xsize = x.size();
	int ysize = y.size();
	int zsize = z.size();
	int count = 0;
	int b = 0;
	for (int iz = 0; iz < zsize - 1; ++iz) {
		count = (int)iz * xsize * ysize;
		for (int iy = 0; iy < ysize - 1; ++iy) {
			count += (int)iy * xsize;
			for (int ix = 0; ix < xsize - 1; ++ix) {
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * iz + (xsize - 1) * iy,
					make_pair(ix + count, ix + count + 1)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * iz + (xsize - 1) * iy + (xsize - 1),
					make_pair(ix + xsize + count, ix + xsize + count + 1)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize + (xsize - 1) * ysize * iz + (xsize - 1) * iy,
					make_pair(ix + count + xsize * ysize, ix + count + xsize * ysize + 1)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize + (xsize - 1) * ysize * iz + (xsize - 1) * iy + (xsize - 1),
					make_pair(ix + xsize + xsize * ysize + count, ix + xsize + xsize * ysize + count + 1)));

				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + iz * xsize * (ysize - 1) + iy * xsize,
					make_pair(ix + count, ix + xsize + count)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + iz * xsize * (ysize - 1) + iy * xsize + 1,
					make_pair(ix + count + 1, ix + xsize + count + 1)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + iz * xsize * (ysize - 1) + iy * xsize + xsize,
					make_pair(ix + count + xsize * ysize, ix + xsize + xsize * ysize + count)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + iz * xsize * (ysize - 1) + iy * xsize + xsize + 1,
					make_pair(ix + count + xsize * ysize + 1, ix + xsize + xsize * ysize + count + 1)));

				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * (ysize - 1) * zsize + iz * xsize * ysize + iy * xsize,
					make_pair(ix + count, ix + count + xsize * ysize)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * (ysize - 1) * zsize + iz * xsize * ysize + iy * xsize + 1,
					make_pair(ix + count + 1, ix + count + xsize * ysize + 1)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * (ysize - 1) * zsize + iz * xsize * ysize + iy * xsize + xsize,
					make_pair(ix + xsize + count, ix + xsize + xsize * ysize + count)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * (ysize - 1) * zsize + iz * xsize * ysize + iy * xsize + xsize + 1,
					make_pair(ix + xsize + count + 1, ix + xsize + xsize * ysize + count + 1)));

				fout << aEdges[0].first << ' '
					<< aEdges[0].second.first << ' '
					<< aEdges[0].second.second << endl
					<< aEdges[1].first << ' '
					<< aEdges[1].second.first << ' '
					<< aEdges[1].second.second << endl
					<< aEdges[2].first << ' '
					<< aEdges[2].second.first << ' '
					<< aEdges[2].second.second << endl
					<< aEdges[3].first << ' '
					<< aEdges[3].second.first << ' '
					<< aEdges[3].second.second << endl
					<< aEdges[4].first << ' '
					<< aEdges[4].second.first << ' '
					<< aEdges[4].second.second << endl
					<< aEdges[5].first << ' '
					<< aEdges[5].second.first << ' '
					<< aEdges[5].second.second << endl
					<< aEdges[6].first << ' '
					<< aEdges[6].second.first << ' '
					<< aEdges[6].second.second << endl
					<< aEdges[7].first << ' '
					<< aEdges[7].second.first << ' '
					<< aEdges[7].second.second << endl
					<< aEdges[8].first << ' '
					<< aEdges[8].second.first << ' '
					<< aEdges[8].second.second << endl
					<< aEdges[9].first << ' '
					<< aEdges[9].second.first << ' '
					<< aEdges[9].second.second << endl
					<< aEdges[10].first << ' '
					<< aEdges[10].second.first << ' '
					<< aEdges[10].second.second << endl
					<< aEdges[11].first << ' '
					<< aEdges[11].second.first << ' '
					<< aEdges[11].second.second << endl << endl;
				tmpEdges.push_back(aEdges);
				aEdges.clear();
			}
			count = (int)iz * xsize * ysize;
		}
	}
	fout.close();
	return tmpEdges;
}

vector<vector<double>> locEdgesG(double coef, double hx, double hy, double hz)
{
	vector<vector<double>> tmp;
	tmp.resize(12);
	int size = tmp.size();
	for (int i = 0; i < size; ++i)
		tmp[i].resize(12);

	double constxy_z = hx * hy / (6 * hz * coef);
	double constxz_y = hx * hz / (6 * hy * coef);
	double constyz_x = hz * hy / (6 * hx * coef);
	double constx = -hx / (6 * coef);
	double consty = hy / (6 * coef);
	double constz = -hz / (6 * coef);
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			tmp[i][j] = constxy_z * G1[i][j] + constxz_y * G2[i][j];
			tmp[i + 4][j + 4] = constxy_z * G1[i][j] + constyz_x * G2[i][j];
			tmp[i + 8][j + 8] = constxz_y * G1[i][j] + constyz_x * G2[i][j];

			tmp[i][j + 4] = constz * G2[i][j];
			tmp[i + 4][j] = constz * G2[i][j];
			tmp[i][j + 8] = consty * G3[i][j];
			tmp[i + 8][j] = consty * G3T[i][j];
			tmp[i + 4][j + 8] = constx * G1[i][j];
			tmp[i + 8][j + 4] = constx * G1[i][j];
		}
	}
	return tmp;
}

vector<vector<double>> locEdgesM(double coef, double omega, double hx, double hy, double hz) {
	vector<vector<double>> locM;
	for (int i = 0; i < 4; ++i) {
		vector<double> tmp;
		for (int j = 0; j < 4; ++j) {
			tmp.push_back(D[i][j] * coef * omega * hx * hy * hz / 36);
		}
		locM.emplace_back(tmp);
	}
	return locM;
}

vector<vector<vector<double>>> createLocalMG(double mu, double sigma, double omega, double hx, double hy, double hz) { // Создание локальной матрицы масс
	vector<vector<vector<double>>> locMG;// общая матрица, которая хранит только нижний треугольник (в виде блочных элементов)
	//создаем локальную матрицу жесткости
	vector<vector<double>> locG = locEdgesG(mu, hx, hy, hz);

	// создаем локальную матрицу масс
	vector<vector<double>> locM = locEdgesM(sigma, omega, hx, hy, hz);
	
	//заполняем общую матрицу элементами матрицы жесткости
	for (int i = 0; i < 12; ++i) {
		vector<vector<double>> tmp;
		for(int j = 0; j <= i; ++j){
			vector<double> tmp2;
			tmp2.push_back(locG[i][j]);
			tmp2.push_back(0);
			tmp.emplace_back(tmp2);
		}
		locMG.emplace_back(tmp);
	}
	//заполняем общую матрицу элементами матрицы масс
	for (int k = 0; k < 3; ++k) {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j <= i; ++j) {
				locMG[i + 4 * k][j + 4 * k][1] = locM[i][j];
			}
		}
	}
	return locMG;
}

bool isZero(vector<double> &a) {
	for (int i = 0; i < 2; ++i) {
		if (a[i] != 0) {
			return false;
		}
	}
	return true;
}

bool compareLine(pair<vector<double>, pair<int, int>>& a, pair<vector<double>, pair<int, int>>& b) {
	return a.second.first < b.second.first;
}

bool compareColumn(pair<vector<double>, pair<int, int>>& a, pair<vector<double>, pair<int, int>>& b) {
	return a.second.second < b.second.second;
}

matrix createGlobalMG(vector<double> mu, vector<double> sigma, double omega, vector<vector<pair<int, pair<int, int>>>> nvtr, vector<vector<double>> xyz, int n) {
	matrix res;

	vector<pair<vector<double>, pair<int, int>>> tmp;

	vector<vector<double>> di;
	di.resize(n);
	for (int i = 0; i < n; ++i) {
		di[i].push_back(0);
		di[i].push_back(0);
	}
	
	vector<int> li;
	li.resize(n + 1);
	for (int i = 0; i < li.size(); ++i) {
		li[i] = 0;
	}
	vector<vector<double>> al;
	vector<int> lj;


	for (int k = 0; k < nvtr.size(); ++k) {
		vector<vector<vector<double>>> locMG = createLocalMG(mu[k], sigma[k], omega,
			abs(xyz[nvtr[k][0].second.first][0] - xyz[nvtr[k][0].second.second][0]),
			abs(xyz[nvtr[k][4].second.first][1] - xyz[nvtr[k][4].second.second][1]),
			abs(xyz[nvtr[k][8].second.first][2] - xyz[nvtr[k][8].second.second][2]));

		for (int i = 0; i < 12; ++i) {
			int count = 0;
			for (int j = 0; j <= i; ++j) {
				if (i == j) {
					di[nvtr[k][i].first][0] += locMG[i][j][0];
					di[nvtr[k][i].first][1] += locMG[i][j][1];
					continue;
				}
				if (isZero(locMG[i][j])) {
					continue;
				}
				auto res = find_if(tmp.begin(), tmp.end(), [b = make_pair(nvtr[k][i].first, nvtr[k][j].first)](pair<vector<double>, pair<int, int>> a) {
					if (a.second.first == b.first && a.second.second == b.second) {
						return true;
					}
					return false;
				});
				if (res != tmp.end()) {
					res->first = res->first + locMG[i][j];
				}
				else {
					tmp.emplace_back(make_pair(locMG[i][j], make_pair(nvtr[k][i].first, nvtr[k][j].first)));
					count++;
				}
			}
			for (int f = nvtr[k][i].first + 1; f < li.size(); ++f) {
				li[f] += count;
			}
		}
	}

	std::sort(tmp.begin(), tmp.end(), compareLine);
	int i = 1;
	while (i < tmp.size()) {
		if (tmp[i].second.first == tmp[0].second.first) {
			i++;
			continue;
		}
		else
		{
			std::sort(tmp.begin(), tmp.begin() + i - 1, compareColumn);
			for (int j = 0; j < i; ++j) {
				al.push_back(tmp[j].first);
				lj.push_back(tmp[j].second.second);
			}
			tmp.erase(tmp.begin(), tmp.begin() + i);
			i = 1;
		}
	}

	std::sort(tmp.begin(), tmp.end(), compareColumn);
	for (auto& elem : tmp) {
		al.push_back(elem.first);
		lj.push_back(elem.second.second);
	}
	tmp.clear();

	res.di = di;
	res.al = al;
	res.li = li;
	res.lj = lj;
	return res;
}

vector<double> MultMatrixToVector(vector<vector<double>> matrix, vector<double> vector1)
{
	vector<double> out;
	out.resize(12);
	for (int i = 0; i < 12; ++i) {
		out[i] = 0;
	}
	for (int ix = 0; ix < 4; ix++) {
		for (int jx = 0; jx < 4; jx++)
			out[ix] += matrix[ix][jx] * vector1[jx];
	}
	for (int ix = 4; ix < 8; ix++) {
		for (int jx = 0; jx < 4; jx++)
			out[ix] += matrix[ix - 4][jx] * vector1[jx + 4];
	}
	for (int ix = 8; ix < 12; ix++) {
		for (int jx = 0; jx < 4; jx++)
			out[ix] += matrix[ix - 8][jx] * vector1[jx + 8];
	}
	return out;
}

vector<double> funcSin(double x, double y, double z) {
	vector<double> xyz;
	xyz.resize(3);
	xyz[0] = 0;
	xyz[1] = x * x;
	xyz[2] = 0;
	return xyz;
}

vector<double> funcCos(double x, double y, double z) {
	vector<double> xyz;
	xyz.resize(3);
	xyz[0] = 0;
	xyz[1] = 0;
	xyz[2] = 0;
	return xyz;
}

vector<double> funcFSin(double x, double y, double z, double mu, double sigma, double omega) {
	vector<double> xyz;
	xyz.resize(3);
	xyz[0] = 0;
	xyz[1] = -2 / mu;
	xyz[2] = 0;
	return xyz;
}

vector<double> funcFCos(double x, double y, double z, double mu, double sigma, double omega) {
	vector<double> xyz;
	xyz.resize(3);
	xyz[0] = 0;
	xyz[1] = sigma * omega * x * x;
	xyz[2] = 0;
	return xyz;
}

vector<double> createFSin(vector<vector<pair<int, pair<int, int>>>> nvtr, vector<vector<double>> xyz, vector<double> mu, vector<double> sigma, double omega, int n) {
	vector<double> FSin;
	FSin.resize(n);
	vector<double> locF;
	vector<vector<double>> locM;
	locF.resize(12);
	for (int k = 0; k < nvtr.size(); ++k) {
		locM = locEdgesM(1, 1,
			abs(xyz[nvtr[k][0].second.first][0] - xyz[nvtr[k][0].second.second][0]),
			abs(xyz[nvtr[k][4].second.first][1] - xyz[nvtr[k][4].second.second][1]),
			abs(xyz[nvtr[k][8].second.first][2] - xyz[nvtr[k][8].second.second][2]));
		for (int i = 0; i < 12; ++i) {
			locF[i] = funcFSin((xyz[nvtr[k][i].second.second][0] + xyz[nvtr[k][i].second.first][0]) / 2.,
								(xyz[nvtr[k][i].second.second][1] + xyz[nvtr[k][i].second.first][1]) / 2.,
								(xyz[nvtr[k][i].second.second][2] + xyz[nvtr[k][i].second.first][2]) / 2., mu[k], sigma[k], omega)[i / 4];
		}
		locF = MultMatrixToVector(locM, locF);
		for (int i = 0; i < 12; ++i) {
			FSin[nvtr[k][i].first] += locF[i];
		}
	}

	return FSin;
}

vector<double> createFCos(vector<vector<pair<int, pair<int, int>>>> nvtr, vector<vector<double>> xyz, vector<double> mu, vector<double> sigma, double omega, int n) {
	vector<double> FCos;
	FCos.resize(n);
	vector<double> locF;
	vector<vector<double>> locM;
	locF.resize(12);
	for (int k = 0; k < nvtr.size(); ++k) {
		locM = locEdgesM(1, 1,
			abs(xyz[nvtr[k][0].second.first][0] - xyz[nvtr[k][0].second.second][0]),
			abs(xyz[nvtr[k][4].second.first][1] - xyz[nvtr[k][4].second.second][1]),
			abs(xyz[nvtr[k][8].second.first][2] - xyz[nvtr[k][8].second.second][2]));
		for (int i = 0; i < 12; ++i) {
			locF[i] = funcFCos((xyz[nvtr[k][i].second.second][0] + xyz[nvtr[k][i].second.first][0]) / 2.,
								(xyz[nvtr[k][i].second.second][1] + xyz[nvtr[k][i].second.first][1]) / 2.,
								(xyz[nvtr[k][i].second.second][2] + xyz[nvtr[k][i].second.first][2]) / 2., mu[k], sigma[k], omega)[i / 4];
		}
		locF = MultMatrixToVector(locM, locF);
		for (int i = 0; i < 12; ++i) {
			FCos[nvtr[k][i].first] += locF[i];
		}
	}

	return FCos;
}

vector<vector<double>> createF(vector<vector<pair<int, pair<int, int>>>> nvtr, vector<vector<double>> xyz, vector<double> mu, vector<double> sigma, double omega, int n) {
	vector<double> FCos = createFCos(nvtr, xyz, mu, sigma, omega, n);
	vector<double> FSin = createFSin(nvtr, xyz, mu, sigma, omega, n);
	vector<vector<double>> res;
	for (int i = 0; i < n; ++i) {
		res.push_back({ FSin[i], FCos[i] });
	}
	return res;
}

void conditions(matrix &matr, vector<vector<double>>& F, vector<vector<pair<int, pair<int, int>>>> nvtr, vector<vector<double>> xyz, const vector<double>& x, const vector<double>& y, const vector<double>& z) {

	for (int elemNum = 0; elemNum < nvtr.size(); ++elemNum) {
		for (int i = 0; i < 12; ++i) {
			if ((xyz[nvtr[elemNum][i].second.first][2] == z[0] &&
				xyz[nvtr[elemNum][i].second.second][2] == z[0]) ||

				(xyz[nvtr[elemNum][i].second.first][2] == z[z.size() - 1] &&
					xyz[nvtr[elemNum][i].second.second][2] == z[z.size() - 1]) ||

				(xyz[nvtr[elemNum][i].second.first][0] == x[0] &&
					xyz[nvtr[elemNum][i].second.second][0] == x[0]) ||

				(xyz[nvtr[elemNum][i].second.first][0] == x[x.size() - 1] &&
					xyz[nvtr[elemNum][i].second.second][0] == x[x.size() - 1]) ||

				(xyz[nvtr[elemNum][i].second.first][1] == y[0] &&
					xyz[nvtr[elemNum][i].second.second][1] == y[0]) ||

				(xyz[nvtr[elemNum][i].second.first][1] == y[y.size() - 1] &&
					xyz[nvtr[elemNum][i].second.second][1] == y[y.size() - 1])) {
				if (xyz[nvtr[elemNum][i].second.first][1] != xyz[nvtr[elemNum][i].second.second][1])
					F[nvtr[elemNum][i].first][0] = funcSin((xyz[nvtr[elemNum][i].second.second][0] + xyz[nvtr[elemNum][i].second.first][0]) / 2.,
						(xyz[nvtr[elemNum][i].second.second][1] + xyz[nvtr[elemNum][i].second.first][1]) / 2.,
						(xyz[nvtr[elemNum][i].second.second][2] + xyz[nvtr[elemNum][i].second.first][2]) / 2.)[1] * 1e10;
				else
					F[nvtr[elemNum][i].first][0] = 0;
				F[nvtr[elemNum][i].first][1] = 0;
				matr.di[nvtr[elemNum][i].first][0] = 1e10;
			}
		}
	}
}

double dotProduct(const vector<vector<double>>& a, const vector<vector<double>>& b) {
	if (a.size() != b.size()) {
		exit(1);
	}
	else {
		double res = 0;
		int n = a.size();
		for (int i = 0; i < n; ++i) {
			double tmp =  a[i][0] * b[i][0] + a[i][1] * b[i][1] ;
			res += tmp;
		}
		return res;
	}
}

vector<vector<double>> operator*(const matrix& A, const vector<vector<double>> &b) {
	if (A.di.size() != b.size()) {
		exit(1);
	}
	else {
		vector<vector<double>> res;
		int n = A.di.size();
		res.resize(n);
		for (auto& elem : res) {
			elem.push_back(0);
			elem.push_back(0);
		}

		for (int i = 0; i < n; ++i) {
			res[i] = res[i] + A.di[i] * b[i];
		}

		for (int i = 1; i < n; ++i) {
			int elemPos = A.li[i];
			for (; elemPos < A.li[i + 1]; ++elemPos) {
				res[i] = res[i] + A.al[elemPos] * b[A.lj[elemPos]];
				res[A.lj[elemPos]] = res[A.lj[elemPos]] + A.al[elemPos] * b[i];
			}
		}
		return res;
	}
}

vector<vector<double>> operator-(const vector<vector<double>> &a, const vector<vector<double>> &b) {
	if (a.size() != b.size()) {
		exit(1);
	}
	else {
		vector<vector<double>> res;
		int n = a.size();
		res.resize(n);
		for (int i = 0; i < n; ++i) {
			res[i].push_back(a[i][0] - b[i][0]);
			res[i].push_back(a[i][1] - b[i][1]);
		}

		return res;
	}
}

vector<vector<double>> operator*(const double& coef, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = b.size();
	for (int i = 0; i < n; ++i) {
		res.push_back({ coef * b[i][0], coef * b[i][1] });
	}
	return res;
}

vector<vector<double>> operator*(const vector<double>& coef, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = b.size();
	for (int i = 0; i < n; ++i) {
		res.push_back(coef * b[i]);
	}
	return res;
}

vector<vector<double>> operator+(const vector<vector<double>>& a, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = b.size();
	res.resize(n);
	for (int i = 0; i < n; ++i) {
		res[i] = a[i] + b[i];
	}
	return res;
}

void los(matrix &A, vector<vector<double>>& b, vector<vector<double>>& x0, double eps) {
	vector<vector<double>> r0, z0, p0;

	r0 = b - A * x0;
	p0 = r0;
	z0 = A * p0;

	double absB = dotProduct(b, b);
	double absR0 = dotProduct(p0, p0);
	double nevyazka = absR0 / absB;
	while (abs(nevyazka) > eps) {
		cout << abs(nevyazka) << endl;
		double scalarZ0 = dotProduct(z0, z0);

		double alphaK = dotProduct(z0, r0) / scalarZ0;

		absR0 = absR0 - alphaK * alphaK * scalarZ0;
		nevyazka = absR0 / absB;

		x0 = x0 + alphaK * p0;
		r0 = r0 - alphaK * z0;
		double betaK = dotProduct(z0, A * r0) / scalarZ0;
		p0 = r0 + betaK * p0;
		z0 = A * r0 + betaK * z0;
	}
	cout << abs(nevyazka) << endl;
}

int main() {
	ifstream fin("x.txt");
	vector<double> x, y, z, mu, omega, sigma;
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		x.emplace_back(tmp);
	}
	fin.close();
	fin.open("y.txt");
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		y.emplace_back(tmp);
	}
	fin.close();
	fin.open("z.txt");
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		z.emplace_back(tmp);
	}
	fin.close();
	fin.open("mu.txt");
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		mu.push_back(tmp);
	}
	fin.close();
	fin.open("sigma.txt");
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		sigma.push_back(tmp);
	}
	fin.close();
	fin.open("omega.txt");
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		omega.push_back(tmp);
	}
	fin.close();

	double eps = 1e-8;

	int n = (x.size() - 1) * y.size() * z.size() + x.size() * y.size() * (z.size() - 1) + x.size() * (y.size() - 1) * z.size();

	vector<vector<double>> xyz = formxyz(x, y, z);
	vector<vector<pair<int, pair<int, int>>>> nvtr = formnvtr(x, y, z);

	matrix A = createGlobalMG(mu, sigma, omega[0], nvtr, xyz, n);
	vector<vector<double>> F = createF(nvtr, xyz, mu, sigma, omega[0], n);

	conditions(A, F, nvtr, xyz, x, y, z);

	vector<vector<double>> x0;
	for (int i = 0; i < n; ++i) {
		x0.push_back({ 1., 1. });
	}

	los(A, F, x0, eps);

	for (int i = 0; i < x0.size(); ++i) {
		cout << i << ": " << x0[i][0] << "\t\t" << x0[i][1] << endl;
	}
	return 0;
}