
#include <iostream>
#include <fstream>
#include <vector>
#include "QuadMatrix.h"
#include "LinSolveAlgs.h"
#include "Norma.h"
#include "EigenValuesAndVectors.h"

using namespace std;
ofstream ansFile;
const double EPS = 1e-7;

template<class T>
void test() {

    vector<T> fileVector;


    ifstream matrixFile;
    matrixFile.open("DATA.txt");

    while (!matrixFile.eof()) {
        T i;
        matrixFile >> i;
        fileVector.push_back(i);
    }
    //fileVector.pop_back();
    int k = 0;
    //n* (n + 1) = fileVector.size();
    size_t n = (-1 + sqrt(1 + 4 * fileVector.size())) / 2;

    QuadMatrix<T> matrix(n);
    vector<T> b(n);

    if (n * (n + 1) != fileVector.size()) {
        ansFile << "Матрица системы не является квадратной!" << endl;
        system("pause");
        exit(1);
    }

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
        b[i] = fileVector[k];
        ++k;
    }
    matrixFile.close();


    ansFile << "Исходная матрица (A):";
    for (int i = 0; i < n; i++) {
        ansFile << endl;
        for (int j = 0; j < n; j++) {
            ansFile << matrix(i, j) << " ";
        }
    }
    ansFile << endl << endl;

    /*auto [res, C] = gaussLinSolve(matrix, b);

    auto [Q, R] = qrDecomposition(matrix);
    vector<T> resQr = qrLinSolve(Q, R, b);

    ansFile << endl;
    ansFile << "Матрица (Q):";
    for (int i = 0; i < n; i++) {
        ansFile << endl;
        for (int j = 0; j < n; j++) {
            ansFile << Q(i, j) << " ";
        }
    }
    ansFile << endl << endl;

    ansFile << "Матрица (R):";
    for (int i = 0; i < n; i++) {
        ansFile << endl;
        for (int j = 0; j < n; j++) {
            ansFile << R(i, j) << " ";
        }
    }
    ansFile << endl << endl;

    ansFile << "Результат полученный QR-разложением (x):" << endl;
    for (int i = 0; i < n; ++i) {
        ansFile << resQr[i] << endl;
    }*/

    QuadMatrix<T> MatrixEigen = qrSearchEigenvalue(matrix);

    ansFile << "Получившаяся матрица А:";
    for (int i = 0; i < n; i++) {
        ansFile << endl;
        for (int j = 0; j < n; j++) {
            ansFile << MatrixEigen(i, j) << " ";
        }
    }
    ansFile << endl << endl;

    /*ansFile << endl << "Норма вектора невязки в методе QR-разложения (||b - b1||): " << endl;
    ansFile << "При кубической норме: " << normDiscrepancyVectorQR(matrix, b, norm_inf) << endl;
    ansFile << "При октаэдральной норме: " << normDiscrepancyVectorQR(matrix, b, norm_1) << endl;

    ansFile << endl;
    ansFile << "Кубическая норма исходной матрицы: " << norm_inf(matrix) << endl;
    ansFile << "Октаэдральная норма исходной матрицы: " << norm_1(matrix) << endl;

    ansFile << endl;
    ansFile << "Обратная матрица: ";
    QuadMatrix<T> inversMatrix = matrix.inv();
    for (int i = 0; i < n; ++i) {
        ansFile << endl;
        for (int j = 0; j < n; ++j) {
            ansFile << inversMatrix(i, j) << ' ';
        }
    }

    ansFile << endl << endl << "Оценка числа обусловленности: " << endl;

    ansFile << "При кубической норме: " << condEstimate(matrix, norm_inf) << endl;
    ansFile << "При октаэдральной норме: " << condEstimate(matrix, norm_1) << endl;

    ansFile << endl;
    ansFile << "Число обусловленности: " << endl;
    T condMatrixInf = cond(matrix, norm_inf);
    T condMatrix_1 = cond(matrix, norm_1);
    ansFile << "При кубической норме: " << condMatrixInf << endl;
    ansFile << "При октаэдральной норме: " << condMatrix_1 << endl;
    ansFile << endl << endl << endl;*/
}


int main()
{
    setlocale(LC_ALL, "Russian");

    ansFile.open("AnswerFileDATA.txt");
    ansFile << "Точность double:" << endl;
    test<double>();
    ansFile << "Точность float:" << endl;
    test<float>();
    ansFile.close();

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
