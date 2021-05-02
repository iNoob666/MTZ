#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <ctime>

using namespace std;


unsigned int startTime;
unsigned int endTime;

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

struct matrix // Структура для матрицы, которая хранится в
	// разряженном строчно-столбцовом формате относительно блочных элементов (размер блока 2x2)
{
	vector<vector<double>> di;// диагональные элементы
	vector<vector<double>> al;// элементы нижнего тругольника
	vector<int> li;// профиль матрицы
	vector<int> lj;// портрет матрицы
};

struct LUmatrix // Структура для матрицы, которая хранится в
	// разряженном строчно-столбцовом формате относительно блочных элементов (размер блока 2x2)
{
	vector<double> di;// диагональные элементы
	vector<double> al;// элементы нижнего тругольника
	vector<double> au;// элементы верхнего тругольника
	vector<int> li;// профиль матрицы
	vector<int> lj;// портрет матрицы
};

struct material {
	double x0, x1;
	double y0, y1;
	double z0, z1;
	double mu, sigma;
};

vector<double> operator/(const vector<double>& a, const vector<double>& b) {
	double coefB = b[0] * b[0] + b[1] * b[1];
	vector<double> oppositeB = { b[0] / coefB, -(b[1] / coefB) };
	return { a[0] * oppositeB[0] - a[1] * oppositeB[1], a[0] * oppositeB[1] + a[1] * oppositeB[0] };
}

vector<double> operator*(const vector<double>& a, const vector<double>& b) {
	return { a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] };
}

vector<double> operator+(const vector<double>& a, const vector<double>& b) {
	vector<double> res;
	int size = a.size();
	res.resize(size);
	for (int i = 0; i < size; ++i) {
		res[i] = a[i] + b[i];
	}
	return res;
}

vector<double> operator-(const vector<double>& a, const vector<double>& b) {
	return { a[0] - b[0], a[1] - b[1] };
}

vector<vector<double>> operator*(const matrix& A, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = (int)A.di.size();
	res.resize(n);
	for (auto& elem : res) {
		elem.emplace_back(0);
		elem.emplace_back(0);
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

vector<vector<double>> operator-(const vector<vector<double>>& a, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = (int)a.size();
	res.resize(n);
	for (int i = 0; i < n; ++i) {
		res[i].emplace_back(a[i][0] - b[i][0]);
		res[i].emplace_back(a[i][1] - b[i][1]);
	}

	return res;
}

vector<vector<double>> operator*(const double& coef, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = (int)b.size();
	res.resize(n);
	for (int i = 0; i < n; ++i) {
		res[i] =  { coef * b[i][0], coef * b[i][1] };
	}
	return res;
}

vector<vector<double>> operator*(const vector<double>& coef, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = (int)b.size();
	res.reserve(n);
	for (int i = 0; i < n; ++i) {
		res.emplace_back(coef * b[i]);
	}
	return res;
}

vector<vector<double>> operator+(const vector<vector<double>>& a, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	int n = (int)b.size();
	res.resize(n);
	for (int i = 0; i < n; ++i) {
		res[i] = a[i] + b[i];
	}
	return res;
}

vector<double> operator*(const vector<double>& b, const double& coef) {
	vector<double> res = b;
	for (auto& elem : res) {
		elem *= coef;
	}
	return res;
}

void readObjects(vector<material>& materials, const char* fileName) {
	ifstream fin(fileName);
	int n;
	material tmp;

	fin >> n;

	for (int i = 0; i < n; ++i) {
		fin >> tmp.x0 >> tmp.x1 >> tmp.y0 >> tmp.y1 >> tmp.z0 >> tmp.z1 >> tmp.mu >> tmp.sigma;
		materials.emplace_back(tmp);
	}
}

void readNet(vector<double>& x, vector<double>& y, vector<double>& z,
	const char* fileName, const vector<material>& materials) {
	ifstream fin(fileName);
	//x
	int nx;

	fin >> nx;

	vector<double> vertexX(nx);
	vector<double> sectionsX(nx - 1);
	vector<double> ratiosX(nx - 1);
	for (int i = 0; i < nx; ++i)
		fin >> vertexX[i];
	for (int i = 0; i < nx - 1; ++i)
		fin >> sectionsX[i];
	for (int i = 0; i < nx - 1; ++i)
		fin >> ratiosX[i];

	double hx, tmpx;

	for (int j = 0; j < nx - 1; j++) {
		if (ratiosX[j] != 1)
			hx = (vertexX[j + 1] - vertexX[j]) * (1. - ratiosX[j]) / (1. - pow(ratiosX[j], sectionsX[j]));
		else hx = (vertexX[j + 1] - vertexX[j]) / sectionsX[j];

		tmpx = vertexX[j];

		for (int k = 0; k < sectionsX[j]; k++) {
			x.emplace_back(tmpx);
			tmpx += hx;
			hx *= ratiosX[j];
		}
	}
	x.emplace_back(tmpx);
	//y
	int ny;

	fin >> ny;

	vector<double> vertexY(ny);
	vector<double> sectionsY(ny - 1);
	vector<double> ratiosY(ny - 1);
	for (int i = 0; i < ny; ++i)
		fin >> vertexY[i];
	for (int i = 0; i < ny - 1; ++i)
		fin >> sectionsY[i];
	for (int i = 0; i < ny - 1; ++i)
		fin >> ratiosY[i];

	double hy, tmpy;

	for (int j = 0; j < ny - 1; j++) {
		if (ratiosY[j] != 1)
			hy = (vertexY[j + 1] - vertexY[j]) * (1. - ratiosY[j]) / (1. - pow(ratiosY[j], sectionsY[j]));
		else hy = (vertexY[j + 1] - vertexY[j]) / sectionsY[j];

		tmpy = vertexY[j];

		for (int k = 0; k < sectionsY[j]; k++) {
			y.emplace_back(tmpy);
			tmpy += hy;
			hy *= ratiosY[j];
		}
	}
	y.emplace_back(tmpy);
	//z
	int nz;

	fin >> nz;

	vector<double> vertexZ(nz);
	vector<double> sectionsZ(nz - 1);
	vector<double> ratiosZ(nz - 1);
	for (int i = 0; i < nz; ++i)
		fin >> vertexZ[i];
	for (int i = 0; i < nz - 1; ++i)
		fin >> sectionsZ[i];
	for (int i = 0; i < nz - 1; ++i)
		fin >> ratiosZ[i];

	double hz, tmpz;

	for (int j = 0; j < nz - 1; j++) {
		if (ratiosZ[j] != 1)
			hz = (vertexZ[j + 1] - vertexZ[j]) * (1. - ratiosZ[j]) / (1. - pow(ratiosZ[j], sectionsZ[j]));
		else hz = (vertexZ[j + 1] - vertexZ[j]) / sectionsZ[j];

		tmpz = vertexX[j];

		for (int k = 0; k < sectionsZ[j]; k++) {
			z.emplace_back(tmpz);
			tmpz += hz;
			hz *= ratiosZ[j];
		}
	}
	z.emplace_back(tmpz);
}

int checkWhichMaterial(double x, double y, double z, const vector<material>& materials) {
	int n = materials.size() - 1;
	for (; n >= 0; --n) {
		if (materials[n].x0 <= x && x <= materials[n].x1 &&
			materials[n].y0 <= y && y <= materials[n].y1 &&
			materials[n].z0 <= z && z <= materials[n].z1)
				return n;
	}
	return -1;
}

void fillMuSigma(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& mu, vector<double>& sigma, const vector<material> & materials) {
	int nx = x.size() - 1, ny = y.size() - 1, nz = z.size() - 1;
	for (int oz = 0; oz < nz; ++oz) {
		for (int oy = 0; oy < ny; ++oy) {
			for (int ox = 0; ox < nx; ++ox) {
				int midX = (x[ox] + x[ox + 1]) / 2;
				int midY = (x[oy] + x[oy + 1]) / 2;
				int midZ = (z[oz] + z[oz + 1]) / 2;
				int check = checkWhichMaterial(midX, midY, midZ, materials);
				if (check != -1) {
					mu.emplace_back(materials[check].mu);
					sigma.emplace_back(materials[check].sigma);
				}
				else {
					mu.emplace_back(0);
					sigma.emplace_back(0);
				}
			}
		}
	}
}

LUmatrix converter(const matrix& base) {
	LUmatrix res;

	size_t size = base.di.size();
	size_t newSize = size * 2;

	res.di.resize(newSize);
	res.li.resize(newSize + 1);

	res.li[0] = 0;
	res.li[1] = 0;

	vector<double> al1;
	vector<double> al2;
	vector<double> au1;
	vector<double> au2;
	vector<int> li;
	vector<int> lj1;
	vector<int> lj2;

	for (int line = 0; line < size; ++line) {
		res.di[line * 2] = base.di[line][0];
		res.di[line * 2 + 1] = base.di[line][0];

		int i0 = base.li[line];
		int i1 = base.li[line + 1];
		for (int elemN = i0; elemN < i1; ++elemN) {
			int column = base.lj[elemN];

			al1.emplace_back(base.al[elemN][0]);
			al1.emplace_back(-base.al[elemN][1]);
			al2.emplace_back(base.al[elemN][1]);
			al2.emplace_back(base.al[elemN][0]);

			au1.emplace_back(base.al[elemN][0]);
			au1.emplace_back(base.al[elemN][1]);
			au2.emplace_back(-base.al[elemN][1]);
			au2.emplace_back(base.al[elemN][0]);

			lj1.emplace_back(column * 2);
			lj1.emplace_back(column * 2 + 1);
			lj2.emplace_back(column * 2);
			lj2.emplace_back(column * 2 + 1);
		}
		//диагональные элементы
		al2.emplace_back(base.di[line][1]);
		au2.emplace_back(-base.di[line][1]);
		lj2.emplace_back(line * 2);
		//заполни массив li
		unsigned int count1 = al1.size();
		unsigned int count2 = al2.size();
		res.li[line * 2 + 1] += (int)count1;
		for (int i = line * 2 + 1; i < newSize; ++i) {
			res.li[i + 1] += (int)count1;
			res.li[i + 1] += (int)count2;
		}
		//перенос в результирующую структуру, и очищение массивов
		res.al.insert(res.al.end(), al1.begin(), al1.end());
		res.al.insert(res.al.end(), al2.begin(), al2.end());
		res.au.insert(res.au.end(), au1.begin(), au1.end());
		res.au.insert(res.au.end(), au2.begin(), au2.end());
		res.lj.insert(res.lj.end(), lj1.begin(), lj1.end());
		res.lj.insert(res.lj.end(), lj2.begin(), lj2.end());
		al1.clear();
		al2.clear();
		au1.clear();
		au2.clear();
		lj1.clear();
		lj2.clear();
	}

	return res;
}

vector<vector<double>> forward(const LUmatrix& A, const vector<vector<double>>& b) {
	vector<vector<double>> res;
	res.resize(b.size());
	for (auto& elem : res) {
		elem.resize(2);
	}

	int size = (int)A.di.size();
	for (int i = 0; i < size; ++i) {
		double sum = 0;
		int i0 = A.li[i], i1 = A.li[i + 1];
		for (int k = i0; k < i1; ++k) {
			int j = A.lj[k];
			sum += A.al[k] * res[j / 2][j % 2];
		}
		res[i / 2][i % 2] = (b[i / 2][i % 2] - sum);
	}
	return res;
}

vector<vector<double>> backward(const LUmatrix& A, const vector<vector<double>>& b) {
	vector<vector<double>> res = b;

	int size = (int)A.di.size();
	for (int i = size - 1; i >= 0; --i) {
		int i0 = A.li[i], i1 = A.li[i + 1] - 1;
		res[i / 2][i % 2] /= A.di[i];
		for (int k = i1; k >= i0; --k) {
			int j = A.lj[k];
			res[j / 2][j % 2] -= A.au[k] * res[i / 2][i % 2];
		}
	}

	return res;
}

LUmatrix luFactorization(const LUmatrix& base) {
	LUmatrix res;

	size_t size = base.di.size();

	res.di.resize(size);

	res.li = base.li;
	res.lj = base.lj;
	res.al.resize(base.al.size());
	res.au.resize(base.au.size());

	for (int line = 0; line < size; ++line) {
		double sumD = 0;

		int i0 = base.li[line];
		int i1 = base.li[line + 1];
		for (int elemN = i0; elemN < i1; ++elemN) {
			double sumL = 0, sumU = 0;

			int column = base.lj[elemN];

			int j0 = base.li[column];
			int j1 = base.li[column + 1];

			int kl = i0, ku = j0;

			for (; kl < i1 && ku < j1;) {
				int j_kl = base.lj[kl];
				int j_ku = base.lj[ku];

				if (j_kl == j_ku) {
					sumU += res.al[kl] * res.au[ku];
					sumL += res.au[kl] * res.al[ku];
					kl++;
					ku++;
				}
				if (j_kl > j_ku)
					ku++;
				if (j_kl < j_ku)
					kl++;
			}
			res.au[elemN] = base.au[elemN] - sumL;
			res.al[elemN] = base.al[elemN] - sumU;
			res.al[elemN] /= res.di[column];

			sumD += res.al[elemN] * res.au[elemN];
		}
		res.di[line] = base.di[line] - sumD;
	}

	return res;
}

double dotProduct(const vector<vector<double>>& a, const vector<vector<double>>& b) {
	double res = 0;
	int n = (int)a.size();
	for (int i = 0; i < n; ++i) {
		double tmp = a[i][0] * b[i][0] + a[i][1] * b[i][1];
		res += tmp;
	}
	return res;
}

void los(matrix& A, vector<vector<double>>& b, vector<vector<double>>& x0, double eps) {
	startTime = clock();
	LUmatrix temp = converter(A);
	endTime = clock();
	cout << "Converting time: " << endTime - startTime << endl;

	startTime = clock();
	LUmatrix LU = luFactorization(temp);
	endTime = clock();
	cout << "LU time: " << endTime - startTime << endl;

	int maxitter = 100000;
	int itterCount = 0;
	vector<vector<double>> r0 = b - A * x0;
	r0 = forward(LU, r0);
	vector<vector<double>> z0 = backward(LU, r0);
	vector<vector<double>> p0 = forward(LU, A * z0);

	double nevyazka = dotProduct(r0, r0);

	while (/*itterCount < maxitter && */nevyazka >= eps) {
		cout << nevyazka << endl;
		double dotP = dotProduct(p0, p0);
		double alphaK = dotProduct(p0, r0) / dotP;
		x0 = x0 + alphaK * z0;
		r0 = r0 - alphaK * p0;

		nevyazka = dotProduct(r0, r0);

		vector<vector<double>> LAU = backward(LU, r0);
		vector<vector<double>> Ur0 = LAU;
		LAU = A * LAU;
		LAU = forward(LU, LAU);
		double betaK = -dotProduct(p0, LAU) / dotP;

		z0 = Ur0 + betaK * z0;
		p0 = LAU + betaK * p0;

		itterCount++;
	}
	cout << nevyazka << endl;
	cout << "count of itterations: " << itterCount << endl;
}

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
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + iz * xsize * (ysize - 1) + iy * xsize + xsize * (ysize - 1),
					make_pair(ix + count + xsize * ysize, ix + xsize + xsize * ysize + count)));
				aEdges.emplace_back(make_pair(ix + (xsize - 1) * ysize * zsize + iz * xsize * (ysize - 1) + iy * xsize + xsize * (ysize - 1) + 1,
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
				tmpEdges.emplace_back(aEdges);
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
			tmp.emplace_back(D[i][j] * coef * omega * hx * hy * hz / 36);
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
			tmp2.emplace_back(locG[i][j]);
			tmp2.emplace_back(0);
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
		di[i].emplace_back(0);
		di[i].emplace_back(0);
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
	int i = 1, j = 0;
	while (i < tmp.size()) {
		if (tmp[i].second.first == tmp[j].second.first) {
			i++;
			continue;
		}
		else
		{
			std::sort(tmp.begin() + j, tmp.begin() + i, compareColumn);
			for (int k = j; k < i; ++k) {
				al.emplace_back(tmp[k].first);
				lj.emplace_back(tmp[k].second.second);
			}
			j = i;
		}
	}

	std::sort(tmp.begin() + j, tmp.end(), compareColumn);
	for (auto elem = tmp.begin() + j; elem < tmp.end(); ++elem) {
		al.emplace_back(elem->first);
		lj.emplace_back(elem->second.second);
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
	xyz[1] = 1. / x;
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
	xyz[1] = -2 / mu / ( x * x * x );
	xyz[2] = 0;
	return xyz;
}

vector<double> funcFCos(double x, double y, double z, double mu, double sigma, double omega) {
	vector<double> xyz;
	xyz.resize(3);
	xyz[0] = 0;
	xyz[1] = sigma * omega * ( 1. / x );
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
	res.resize(n);
	for (int i = 0; i < n; ++i) {
		res[i] = { FSin[i], FCos[i] };
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

vector<double> resInPoint(double x, double y, double z, vector<vector<pair<int, pair<int, int>>>> nvtr, vector<vector<double>> xyz, vector<vector<double>> F) {
	vector<double> res;

	for (auto& elem : nvtr) {
		if ((xyz[elem[0].second.first][0] <= x && x <= xyz[elem[0].second.second][0]) &&
			(xyz[elem[4].second.first][1] <= y && y <= xyz[elem[4].second.second][1]) &&
			(xyz[elem[8].second.first][2] <= z && z <= xyz[elem[8].second.second][2])) {

			double hx = xyz[elem[0].second.second][0] - xyz[elem[0].second.first][0];
			double hy = xyz[elem[4].second.second][1] - xyz[elem[4].second.first][1];
			double hz = xyz[elem[8].second.second][2] - xyz[elem[8].second.first][2];
			double Xm = (xyz[elem[0].second.second][0] - x) / hx;
			double Xp = (x - xyz[elem[0].second.first][0]) / hx;
			double Ym = (xyz[elem[4].second.second][1] - y) / hy;
			double Yp = (y - xyz[elem[4].second.first][1]) / hy;
			double Zm = (xyz[elem[8].second.second][2] - z) / hz;
			double Zp = (z - xyz[elem[8].second.first][2]) / hz;

			vector<double> phi0 = {Ym * Zm, 0, 0};
			vector<double> phi1 = {Yp * Zm, 0, 0};
			vector<double> phi2 = {Ym * Zp, 0, 0};
			vector<double> phi3 = {Yp * Zp, 0, 0};
			vector<double> phi4 = {0, Xm * Zm, 0};
			vector<double> phi5 = {0, Xp * Zm, 0};
			vector<double> phi6 = {0, Xm * Zp, 0};
			vector<double> phi7 = {0, Xp * Zp, 0};
			vector<double> phi8 = {0, 0, Xm * Ym};
			vector<double> phi9 = {0, 0, Xp * Ym };
			vector<double> phi10 = {0, 0, Xm * Yp };
			vector<double> phi11 = {0, 0, Xp * Yp };

			double xPoint = (xyz[elem[0].second.second][0] + xyz[elem[0].second.first][0]) / 2;
			double yPoint = (xyz[elem[4].second.second][1] + xyz[elem[4].second.first][1]) / 2;
			double zPoint = (xyz[elem[8].second.second][2] + xyz[elem[8].second.first][2]) / 2;

			res = phi0 * F[elem[0].first][0] + phi1 * F[elem[1].first][0] + phi2 * F[elem[2].first][0] + phi3 * F[elem[3].first][0] +
				phi4 * F[elem[4].first][0] + phi5 * F[elem[5].first][0] + phi6 * F[elem[6].first][0] + phi7 * F[elem[7].first][0] +
				phi8 * F[elem[8].first][0] + phi9 * F[elem[9].first][0] + phi10 * F[elem[10].first][0] + phi11 * F[elem[11].first][0];
			break;
		}
	}

	return res;
}

int main() {
	vector<material> materials;
	vector<double> x, y, z, mu, omega, sigma;
	readObjects(materials, "materials.txt");
	readNet(x, y, z, "net.txt", materials);
	fillMuSigma(x, y, z, mu, sigma, materials);
	ifstream fin("omega.txt");
	while (!fin.eof()) {
		double tmp;
		fin >> tmp;
		omega.emplace_back(tmp);
	}
	fin.close();

	double eps = 1e-16;

	int n = (x.size() - 1) * y.size() * z.size() + x.size() * y.size() * (z.size() - 1) + x.size() * (y.size() - 1) * z.size();

	vector<vector<double>> xyz = formxyz(x, y, z);
	vector<vector<pair<int, pair<int, int>>>> nvtr = formnvtr(x, y, z);

	startTime = clock();
	matrix A = createGlobalMG(mu, sigma, omega[0], nvtr, xyz, n);
	endTime = clock();
	cout << "Creating Global matrix time: " << endTime - startTime << endl;

	vector<vector<double>> F = createF(nvtr, xyz, mu, sigma, omega[0], n);

	conditions(A, F, nvtr, xyz, x, y, z);

	vector<vector<double>> x0;
	x0.resize(n);
	for (int i = 0; i < n; ++i) {
		x0[i] = { 0, 0 };
	}

	los(A, F, x0, eps);

	int numberOfNode = (x.size() - 1) * y.size() * z.size() + x.size() * (y.size() - 1) * (z.size() / 2) + x.size() / 2;

	cout << endl;
	for (int i = 0; i < y.size() - 1; ++i) {
		int div = i * x.size();
		cout << numberOfNode + i * x.size() << ": " << x0[numberOfNode + i * x.size()][0] << endl;
	}
	cout << endl;

	cout << "Result in point(11, 11, 11): " << resInPoint(11, 11, 11, nvtr, xyz, x0)[1] << endl;

	/*for (int i = 2; i < 20; ++i) {
		cout << "Result in point(" << i << ", 5, 5): " << resInPoint(i, 5, 5, nvtr, xyz, x0)[1] << endl;
	}*/

	cout << "Count of elements: " << (x.size() - 1) * (y.size() - 1) * (z.size() - 1) << endl;

	return 0;
}