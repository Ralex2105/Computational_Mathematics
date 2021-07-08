#include <iostream>
#include "MATRIX.H"
#include <windows.h>

#define n 8
#define N_TYPE double

int main() {
    SetConsoleOutputCP(CP_UTF8);
    //Массив из параметров
    N_TYPE pArray[5] = {0.0000010, 0.0000015, 0.0000020, 0.0000025, 0.0000030};

    //Цикл для разных параметров p
    for (double p : pArray) {

        //Параметр pS
        //Матрица A
        double A[8][8] = {{p - 3, -4,  -4,  7,   2,  3,  8,  7},
                          {0,     -15, -1,  5,   -3, 6,  6,  -6},
                          {-4,    2,   -16, 7,   0,  8,  -7, 6},
                          {0,     8,   -5,  -11, 1,  0,  4,  5},
                          {8,     6,   -8,  4,   27, -7, -1, 5},
                          {-4,    -2,  1,   2,   -8, 10, 7,  0},
                          {0,     -1,  5,   2,   -8, 2,  -2, 0},
                          {0,     -8,  -7,  3,   -7, -4, -8, 5}};

        //Переменные для работы с исходной системой
        N_TYPE *x1; //Результат вызова Solve
        N_TYPE cond = 0; //Обусловленность матрицы
        N_TYPE work[n];
        int ipvt[n];


        //Переменные для работы с новой системой;
        double A_LU[8][8]; //Для LU-разложения
        double A_ob[8][8]; //Обратная матрица A-LU
        double A_ob_tr[8][8]; //Транспонированная матрица A_ob
        double M[8][8]; // M = A*A^(-1)
        double R[8][8]{}; //Матрица невязки R = E - A*A^(-1)
        double E[8][8]; //Единичная матрица
        double b[8][8]=     {{1, 0, 0, 0, 0, 0, 0, 0},
                             {0, 1, 0, 0, 0, 0, 0, 0},
                             {0, 0, 1, 0, 0, 0, 0, 0},
                             {0, 0, 0, 1, 0, 0, 0, 0},
                             {0, 0, 0, 0, 1, 0, 0, 0},
                             {0, 0, 0, 0, 0, 1, 0, 0},
                             {0, 0, 0, 0, 0, 0, 1, 0},
                             {0, 0, 0, 0, 0, 0, 0, 1}};

        //Переменные для нахождения нормы матрицы
        double temp;
        double max = 0;


        //Создание A_LU копии исходной матрицы А
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A_LU[i][j] = A[i][j];
            }
        }

        //Создание единичной матрицы
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    E[i][j] = 1;
                    A_ob[i][j] = 1;
                }
                else {
                    E[i][j] = 0;
                    A_ob[i][j] = 0;
                }

            }
        }

        //LU-разложение матрицы
        for (int i = 0; i < n; i++) {
            work[i] = 0;
            ipvt[i] = 0;
        }
        decomp(n, A_LU, &cond, ipvt, work);

        //Вывод значения обусловленности в консоль
        std::cout << "Параметр p = " << p << std::endl;
        std::cout << "Обусловленность матрицы = " << cond << std::endl;

        //Поиск элементов обратной матрицы
        for (int i = 0; i < n; i++) {
            solve(n, A_LU, &A_ob[i][0], ipvt);
        }

        //Траспонируем матрицу
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A_ob_tr[i][j] = A_ob[j][i];
            }
        }

        // M = A*A^(-1)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++) {
                    R[i][j] += A_ob_tr[i][k] * A[k][j];
                }
        }

        //R = E - M
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                R[i][j] = E[i][j] - R[i][j];
        }

        //Вычисление нормы матрицы
        for (auto & i : R){
            temp = 0;
            for (double j : i){
                temp += fabs(j);
            }
            if (temp > max)
                max = temp;
        }

        //Вывод нормы матрицы в консоль
        std::cout << "Норма матрицы невязки R = " << max << std::endl;
        printf("\n");
    }
}