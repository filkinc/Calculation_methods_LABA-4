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
QuadMatrix<T> qrSearchEigenvalue(QuadMatrix<T> A) {
	auto [curQ, curR] = qrDecomposition(A);
	size_t n = A.order();
	//QuadMatrix<T> curA = A;
	QuadMatrix<T> nextA = A;
	while (fabs(nextA(4, 3)) <= 1e-7) {
		nextA = curQ * curR;

		/*for (size_t count = 1; count < 10; count++) {
			ansFile << "Приближение №" << count << "матрицы A" << endl;
			for (int i = 0; i < n; i++) {
				ansFile << endl;
				for (int j = 0; j < n; j++) {
					ansFile << nextA(i, j) << " ";
				}
			}
		}*/

		auto [nextQ, nextR] = qrDecomposition(nextA);
		curQ = nextQ;
		curR = nextR;
	}
	return nextA;
}