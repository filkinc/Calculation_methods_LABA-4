﻿#include <fstream>
#include <iostream>
#include <vector>
#include "QuadMatrix.h"
#include "LinSolveAlgs.h"
#include "Norma.h"
#include "EigenValuesAndVectors.h"

using namespace std;
ofstream ansFile;

template<class T>
void test(T eps) {
    /*vector<double> fileVector;
    ifstream matrixFile;
    matrixFile.open("D1.txt");

    double i;
    while (!matrixFile.eof()) {
        matrixFile >> i;
        fileVector.push_back(i);
    }

    int k = 0;

    QuadMatrix<double> matrix(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (k == fileVector.size()) {
                break;
            }
            else {
                matrix(i, j) = fileVector[k];
                k++;
            }
        }
    }
    matrixFile.close();*/

    QuadMatrix<T> matrix({
        {1.50, 0.00, -0.43, -0.75},
        {0.00, 3.00, 0.87, -0.50},
        {-0.43, 0.87, 2.90, -0.22},
        {-0.75, -0.50, -0.22, 2.60}
        });

    //ansFile.open("AnswerFileDATA.txt");

    size_t n = matrix.order();

    cout << "Исходная матрица (A):";
    for (int i = 0; i < n; i++) {
        cout << endl;
        for (int j = 0; j < n; j++) {
            cout << matrix(i, j) << " ";
        }
    }
    cout << endl << endl;

    cout << "Стандартный метод QR:";
    auto vecEigenValueQRSimple = qrSimpleSearchEigenvalue(matrix, 4, eps);
    cout << endl;

    cout << "Вектор собственных значений метод QR: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << vecEigenValueQRSimple.eigenValues[i] << " " << endl;
    }
    cout << "Количесво итераций QR:" << endl;
    cout << vecEigenValueQRSimple.iterationCount << endl << endl;

    cout << "Метод QR + сдвиг:";
    auto vecEigenValueQRPlusShift = qrPlusShiftSearchEigenvalue(matrix, 4, eps);
    cout << endl;

    cout << "Вектор собственных значений метод QR + сдвиг: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << vecEigenValueQRPlusShift.eigenValues[i] << " " << endl;
    }
    cout << "Количесво итераций QR + сдвиг:" << endl;
    cout << vecEigenValueQRPlusShift.iterationCount << endl << endl;

    cout << "Матрица Хесенберга: ";
    cout << endl << "Врашение № 1";
    auto matrixHesenbergA1 = HesenbergMatrix(matrix, 2, 3);
    for (int i = 0; i < n; ++i) {
        cout << endl;
        for (int j = 0; j < n; ++j) {
            cout << matrixHesenbergA1(i, j) << " ";
        }
    }
    cout << endl;

    cout << endl << "Врашение № 2";
    auto matrixHesenbergA2 = HesenbergMatrix(matrixHesenbergA1, 2, 4);
    for (int i = 0; i < n; ++i) {
        cout << endl;
        for (int j = 0; j < n; ++j) {
            cout << matrixHesenbergA2(i, j) << " ";
        }
    }
    cout << endl;

    cout << endl << "Врашение № 3";
    auto matrixHesenbergFinal = HesenbergMatrix(matrixHesenbergA2, 3, 4);
    for (int i = 0; i < n; ++i) {
        cout << endl;
        for (int j = 0; j < n; ++j) {
            cout << matrixHesenbergFinal(i, j) << " ";
        }
    }
    cout << endl;


    cout << endl << "Метод QR + матрица Хесенберга";
    auto vecEigenValueQRPlusHesenberg = qrHesenbergSearchEigenvalue(matrix, 4, eps);

    cout << endl << "Вектор собственных значений метод QR + матрица Хесенберга: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << vecEigenValueQRPlusHesenberg.eigenValues[i] << " " << endl;
    }

    cout << "Количесво итераций QR + матрица Хесенберга:" << endl;
    cout << vecEigenValueQRPlusHesenberg.iterationCount << endl;


    cout << endl << "Метод QR + сдвиг + матрица Хесенберга";
    auto vecEigenValueQRPlusShiftHesenberg = qrPlusShiftHesenbergSearchEigenvalue(matrix, 4, eps);

    cout << endl << "Вектор собственных значений метод QR + сдвиг + матрица Хесенберга: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << vecEigenValueQRPlusShiftHesenberg.eigenValues[i] << " " << endl;
    }

    cout << "Количесво итераций QR + сдвиг + матрица Хесенберга:" << endl;
    cout << vecEigenValueQRPlusShiftHesenberg.iterationCount << endl << endl;

    cout << "Нахождение собственных векторов методом обратной итерации: " << endl;
    auto vecVectorEigenReverseIteration = eigenVectorReverseIteration(matrix, vecEigenValueQRPlusShift.eigenValues, eps);

    for (int i = 0; i < n; ++i) {
        printVector(vecVectorEigenReverseIteration.eigenVector[i]);
        cout << endl;
    }

    cout << "Количество итераций для нахождения собственных векторов метод обратной итерации:" << endl;
    cout << vecVectorEigenReverseIteration.iterationCount << endl << endl;


    cout << endl << "Нахождение собственных векторов методом обратной итерации + Рэлея: " << endl;
    auto vecVectorEigenReley = eigenVectorReley(matrix, vecEigenValueQRPlusShift.eigenValues, eps);

    cout << endl;
    for (int i = 0; i < n; ++i) {
        printVector(vecVectorEigenReley.eigenVector[i]);
        cout << endl;
    }

    cout << "Количество итераций для нахождения собственных векторов метод обратной итарции + Рэлея:" << endl;
    cout << vecVectorEigenReley.iterationCount << endl;

    auto checkVector = checkEigenVector(matrix, vecVectorEigenReverseIteration.eigenVector, vecEigenValueQRPlusShift.eigenValues);
    cout << endl << "Проверка собственных векторов: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << checkVector[i] << endl;
    }
}


int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Точность double:" << endl;
    test<double>(1e-10);
    cout << endl << endl << "Точность float:" << endl;
    test<float>(1e-5);

    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.