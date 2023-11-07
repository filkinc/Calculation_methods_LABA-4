#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"


template <class T>
QuadMatrix<T> qrSimpleSearchEigenvalue(QuadMatrix<T> A) {
	auto [curQ, curR] = qrDecomposition(A);
	size_t n = 4;
	size_t iter = 0;
	QuadMatrix<T> nextA = A;

	while (true) {
		nextA = curR * curQ;

		//Временно коментим
		/*if (iter < 5) {
			cout << endl << "Вращение №: " << iter;
			for (int i = 0; i < n; i++) {
				cout << endl;
				for (int j = 0; j < n; j++) {
					cout << nextA(i, j) << " ";
				}
			}
			cout << endl;
		}*/
		
		curQ = qrDecomposition(nextA).first;
		curR = qrDecomposition(nextA).second;
		
		++iter;

		if (fabs(nextA(3, 2)) <= 1e-14) {
			break;
		}

	}
	cout << endl << "Число итераций: " << iter << endl;

	cout << "Получившаяся матрица А:";
	for (int i = 0; i < n; i++) {
		cout << endl;
		for (int j = 0; j < n; j++) {
			cout << nextA(i, j) << " ";
		}
	}
	cout << endl << endl; 

	return nextA;
}

template <class T>
QuadMatrix<T> qrPlusShiftSearchEigenvalue(QuadMatrix<T> A) {
	T sigma = A(3,3);
	QuadMatrix<T> E({
		{1.00, 0.00, 0.00, 0.00},
		{0.00, 1.00, 0.00, 0.00},
		{0.00, 0.00, 1.00, 0.00},
		{0.00, 0.00, 0.00, 1.00}
		});

	A = A - sigma * E;

	auto [curQ, curR] = qrDecomposition(A);
	size_t n = 4;
	size_t iter = 0;
	QuadMatrix<T> nextA = A;

	while (true) {

		nextA = (curR * curQ) + (sigma * E);

		if (fabs(nextA(3, 2)) <= 1e-14) {
			break;
		}

		//Временно коментим
		/*if (iter < 5) {
			cout << endl << "Вращение №: " << iter;
			for (int i = 0; i < n; i++) {
				cout << endl;
				for (int j = 0; j < n; j++) {
					cout << nextA(i, j) << " ";
				}
			}
			cout << endl;
		}*/

		nextA = nextA - (sigma * E);
		//sigma = nextA(3, 3);
		curQ = qrDecomposition(nextA).first;
		curR = qrDecomposition(nextA).second;

		++iter;
	}
	cout << endl << "Число итераций: " << iter << endl;

	cout << "Получившаяся матрица А:";
	for (int i = 0; i < n; i++) {
		cout << endl;
		for (int j = 0; j < n; j++) {
			cout << nextA(i, j) << " ";
		}
	}
	cout << endl << endl;

	return nextA;
}