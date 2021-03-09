#include <vector>
#include <fstream>
#include <algorithm>

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


using namespace std;


struct matrix // Структура для матрицы, которая хранится в
	// разряженном строчно-столбцовом формате относительно блочных элементов (размер блока 2x2) 
{
	vector<double[2][2]> di;// диагональные элементы
	vector<double[2][2]> al;// элементы нижнего тругольника
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
				tmpxyz.emplace_back(iz);
				tmpxyz.emplace_back(iy);
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
				aEdges.emplace_back(make_pair(iy + ix + (xsize - 1) * ysize * zsize + iz * xsize * ysize + iy * xsize,
					make_pair(ix + count, ix + count + xsize * ysize)));
				aEdges.emplace_back(make_pair(iy + ix + (xsize - 1) * ysize * zsize + iz * xsize * ysize + iy * xsize + 1,
					make_pair(ix + count + 1, ix + count + xsize * ysize + 1)));
				aEdges.emplace_back(make_pair(iy + ix + (xsize - 1) * ysize * zsize + iz * xsize * ysize + iy * xsize + xsize,
					make_pair(ix + xsize + count, ix + xsize + xsize * ysize + count)));
				aEdges.emplace_back(make_pair(iy + ix + (xsize - 1) * ysize * zsize + iz * xsize * ysize + iy * xsize + xsize + 1,
					make_pair(ix + xsize + count + 1, ix + xsize + xsize * ysize + count + 1)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * ysize * (zsize - 1) + iz * xsize * (ysize - 1) + iy * xsize,
					make_pair(ix + count, ix + xsize + count)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * ysize * (zsize - 1) + iz * xsize * (ysize - 1) + iy * xsize + 1,
					make_pair(ix + count + 1, ix + xsize + count + 1)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * ysize * (zsize - 1) + iz * xsize * (ysize - 1) + iy * xsize + xsize * (ysize - 1),
					make_pair(ix + count + xsize * ysize, ix + xsize + xsize * ysize + count)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + xsize * ysize * (zsize - 1) + iz * xsize * (ysize - 1) + iy * xsize + xsize * (ysize - 1) + 1,
					make_pair(ix + count + xsize * ysize + 1, ix + xsize + xsize * ysize + count + 1)));

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
			tmp[i][j + 4] = tmp[i + 4][j] = constz * G2[i][j];
			tmp[i][j + 8] = consty * G3[i][j];
			tmp[i + 8][j] = consty * G3T[i][j];
			tmp[i + 4][j + 8] = tmp[i + 8][j + 4] = constx * G1[i][j];
		}
	}
	return tmp;
}

vector<vector<double[2][2]>> createLocalMG(double mu, double omega, double hx, double hy, double hz) { // Создание локальной матрицы масс
	vector<vector<double[2][2]>> locMG;// общая матрица, которая хранит только нижний треугольник (в виде блочных элементов)
	//создаем локальную матрицу жесткости
	vector<vector<double>> locG = locEdgesG(mu, hx, hy, hz);

	vector<vector<double>> locM;
	// создаем локальную матрицу масс
	for (int i = 0; i < 4; ++i) {
		vector<double> tmp;
		for (int j = 0; j < 4; ++j) {
			tmp.push_back(D[i][j] * omega * hx * hy * hz / 36);
		}
		locM.emplace_back(tmp);
	}
	//заполняем общую матрицу элементами матрицы жесткости
	for (int i = 0; i < 12; ++i) {
		vector<double[2][2]> tmp;
		for(int j = 0; j < i; ++j){
			tmp.push_back({ { locG[i][i], 0 }, { 0, locG[i][j] } });
		}
		locMG.emplace_back(tmp);
	}
	//заполняем общую матрицу элементами матрицы масс
	for (int k = 0; k < 3; ++k) {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < i; ++j) {
				locMG[i + 4 * k][j + 4 * k][0][1] = -locM[i][j];
				locMG[i + 4 * k][j + 4 * k][1][0] = locM[i][j];
			}
		}
	}
	return locMG;
}

bool isZero(double a[2][2]) {
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			if (a[i][j] != 0) {
				return false;
			}
		}
	}
	return true;
}

bool compareLine(pair<double, pair<int, int>>& a, pair<double, pair<int, int>>& b) {
	return a.second.first < b.second.first;
}

bool compareColumn(pair<double, pair<int, int>>& a, pair<double, pair<int, int>>& b) {
	return a.second.second < b.second.second;
}

matrix createGlobalMG(vector<double> mu, vector<double> omega, vector<vector<pair<int, pair<int, int>>>> nvtr, vector<vector<double>> xyz, int n) {
	matrix res;

	vector<pair<double[2][2], pair<int, int>>> tmp;

	vector<double[2][2]> di;
	di.resize(n);
	
	vector<int> li;
	li.resize(n + 1);
	for (int i = 0; i < li.size(); ++i) {
		li[i] = 0;
	}
	vector<double[2][2]> al;
	vector<int> lj;


	for (int k = 0; k < nvtr.size(); ++k) {
		vector<vector<double[2][2]>> locMG = createLocalMG(mu[k], omega[k],
			abs(xyz[nvtr[k][0].second.first][0] - xyz[nvtr[k][0].second.second][0]),
			abs(xyz[nvtr[k][4].second.first][1] - xyz[nvtr[k][4].second.second][1]),
			abs(xyz[nvtr[k][8].second.first][2] - xyz[nvtr[k][8].second.second][2]));

		for (int i = 0; i < 12; ++i) {
			int count = 0;
			for (int j = 0; j <= i; ++j) {
				if (i == j) {
					for (int line = 0; line < 2; ++line) {
						for (int col = 0; col < 2; ++col) {
							di[nvtr[k][i].first][line][col] += locMG[i][j][line][col];
						}
					}
					continue;
				}
				if (isZero(locMG[i][j])) {
					continue;
				}

				tmp.emplace_back(make_pair(locMG[i][j], make_pair(nvtr[k][i].first, nvtr[k][j].first)));
				count++;
			}
			for (int f = nvtr[k][i].first + 1; f < li.size(); ++f) {
				li[f] += count;
			}
		}
	}

	sort(tmp.begin(), tmp.end(), compareLine);
	int i = 1;
	while (i < tmp.size()) {
		if (tmp[i].second.first == tmp[0].second.first) {
			i++;
			continue;
		}
		else
		{
			sort(tmp.begin(), tmp.begin() + i - 1, compareColumn);
			for (int j = 0; j < i; ++j) {
				al.push_back(tmp[j].first);
				lj.push_back(tmp[j].second.second);
			}
			tmp.erase(tmp.begin(), tmp.begin() + i);
			i = 1;
		}
	}

	sort(tmp.begin(), tmp.end(), compareColumn);
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

int main() {

	return 0;
}