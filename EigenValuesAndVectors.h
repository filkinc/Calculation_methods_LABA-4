#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"
#include "LinSolveAlgs.h"

#define MAX_ITER 10
#define EPS 1e-4

template<class T>
struct EigenValueSearchResults {
	vector<T> eigenValues;
	vector<vector<T>> eigenVector;
	int iterationCount = 0;
};

template <class T>
EigenValueSearchResults<T> qrSimpleSearchEigenvalue(const QuadMatrix<T>& A, size_t n) {
	auto [curQ, curR] = qrDecomposition(A);
	//size_t n = 4;
	size_t iter = 0;
	QuadMatrix<T> nextA = A;
	//vector<T>  vecEigenValue;
	if (n == 1) {
		EigenValueSearchResults<T> res;
		res.eigenValues = { nextA(0, 0) };

		res.iterationCount = iter;

		return res;
	}

		while (true) {
			nextA = curR * curQ;

			//Временно коментим
			if (iter < 5) {
				cout << endl << "Вращение №: " << iter;
				for (int i = 0; i < n; i++) {
					cout << endl;
					for (int j = 0; j < n; j++) {
						cout << nextA(i, j) << " ";
					}
				}
				cout << endl;
			}

			curQ = qrDecomposition(nextA).first;
			curR = qrDecomposition(nextA).second;

			++iter;

			if (fabs(nextA(n - 1,n - 2)) <= EPS) {
				break;
			}
		}
		
		/*QuadMatrix<T> smallA(n - 1);

		for (int i = 0; i < n - 1; ++i) {
			for (int j = 0; j < n - 1; ++j) {
				smallA(i, j) = nextA(i, j);
			}
		}*/

		for (int i = 0; i < n - 1; ++i) {
			for (int j = 0; j < n - 1; ++j) {
				smallA(i, j) = nextA(i, j);
			}
		}
			
		auto res = qrSimpleSearchEigenvalue(smallA, n - 1);
		res.eigenValues.push_back(nextA(n - 1, n - 1));
		res.iterationCount += iter;
		return res;
}

template <class T>
EigenValueSearchResults<T> qrPlusShiftSearchEigenvalue(QuadMatrix<T> A, int n) {
	size_t iter = 0;
	
	if (n == 1) {
		EigenValueSearchResults<T> res;
		res.eigenValues = { A(0, 0) };

		res.iterationCount = iter;

		return res;
	}
	
	T sigma = A(n - 1, n - 1);

	//Создание матрицы Е
	QuadMatrix<T> E(n);
	for (int i = 0; i < n; ++i) {
		E(i, i) = 1;
	}

	A = A - sigma * E;

	auto [curQ, curR] = qrDecomposition(A);
	QuadMatrix<T> nextA = A;

	while (true) {

		nextA = (curR * curQ) + (sigma * E);

		if (fabs(nextA(n - 1, n - 2)) <= EPS) {
			break;
		}

		//Временно коментим
		if (iter < 5) {
			cout << endl << "Вращение №: " << iter;
			for (int i = 0; i < n; i++) {
				cout << endl;
				for (int j = 0; j < n; j++) {
					cout << nextA(i, j) << " ";
				}
			}
			cout << endl;
		}

		sigma = nextA(n - 1, n - 1);
		nextA = nextA - (sigma * E);
		curQ = qrDecomposition(nextA).first;
		curR = qrDecomposition(nextA).second;

		++iter;
	}

	QuadMatrix<T> smallA(n - 1);

	for (int i = 0; i < n - 1; ++i) {
		for (int j = 0; j < n - 1; ++j) {
			smallA(i, j) = nextA(i, j);
		}
	}

	auto res = qrPlusShiftSearchEigenvalue(smallA, n - 1);
	res.eigenValues.push_back(nextA(n - 1, n - 1));
	res.iterationCount += iter;
	return res;
}

template <class T>
QuadMatrix<T> HesenbergMatrix(const QuadMatrix<T>& A, int k, int l) {
	size_t n = A.order();

	QuadMatrix<T> H(n);

	T alpha = (A(k - 1, k - 2)) / sqrt(A(k - 1, k - 2) * A(k - 1, k - 2) + A(l - 1, k - 2) * A(l - 1, k - 2));
	T beta = (A(l - 1, k - 2)) / sqrt(A(k - 1, k - 2) * A(k - 1, k - 2) + A(l - 1, k - 2) * A(l - 1, k - 2));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {	
			if (i == k - 1) {
				H(i, j) = alpha * A(k - 1, j) + beta * A(l - 1, j);
			}
			else if (i == l - 1) {
				H(i, j) = -beta * A(k - 1, j) + alpha * A(l - 1, j);
			}
			else {
				H(i, j) = A(i, j);
			}
		}
	}

	QuadMatrix<T> rH(n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j == k - 1) {
				rH(i, j) = alpha * H(i, k - 1) + beta * H(i, l - 1);
			}
			else if (j == l - 1){
				rH(i, j) = -beta * H(i, k - 1) + alpha * H(i, l - 1);
			} 
			else {
				rH(i, j) = H(i, j);
			}
		}
	}

	/*cout << endl << "Коэффициенты:" << endl;
	cout << "alpha: " << alpha << endl;
	cout << "beta: " << beta << endl; 

	for (int i = 0; i < n; ++i) {
		cout << endl;
		for (int j = 0; j < n; ++j) {
			cout << rH(i, j) << " ";
		}
	}
	cout << endl;*/
	return rH;
}

template <class T>
EigenValueSearchResults<T> qrHesenbergSearchEigenvalue(const QuadMatrix<T>& A, int n) {
	return qrSimpleSearchEigenvalue(A, n);
}

template <class T>
EigenValueSearchResults<T> qrPlusShiftHesenbergSearchEigenvalue(const QuadMatrix<T>& A, int n) {
	return qrPlusShiftSearchEigenvalue(A, n);
}

template <class T>
EigenValueSearchResults<T> eigenVectorReverseIteration(const QuadMatrix<T>& A, vector<T> lambda) {
	int n = A.order();
	int k = 0;
	vector<T> curVecEigen(n);
	vector<T> nextVecEigen(n);
	curVecEigen[0] = 1;

	EigenValueSearchResults<T> res;

	QuadMatrix<T> E(n);
	for (int i = 0; i < n; ++i) {
		E(i, i) = 1;
	}

	for (int i = 0; i < n; ++i) {
		int iter = 0;
		QuadMatrix<T> B = (A - lambda[i] * E);
		auto Binv = B.inv();
		while (iter < MAX_ITER) {
			nextVecEigen = Binv * curVecEigen;
			nextVecEigen = div(nextVecEigen, norm_2(nextVecEigen));
			if ( fabs(fabs(compose(nextVecEigen, curVecEigen)) - 1) <= EPS ) {
				break;
			}
			/*printVector(nextVecEigen);
			cout << endl;*/
			curVecEigen = nextVecEigen;
			++iter;
		}

		if (iter >= MAX_ITER) {
			cout << "больше " << MAX_ITER << " итераций" << endl;
		}
		res.eigenVector.push_back(curVecEigen);
	}
	return res;
}

//template<class T>
//EigenValueSearchResults<T> eigenVectorReley(const QuadMatrix<T>& A, const vector<T>& lambda) {
//
//}
