#include <ios>
#include <windows.h>
#include "rkf45.h"


//Параметры для решения
const int n = 2;
double T = 1; //начало промежутка
double Tout = 2; //конец промежутка
double X[] = {1, -1}; //массив из двух уравнений
double actualX[n]; //значение X

//Параметры для решения методом Эйлера-Коши
double h = 0.1; //Шаг

//Параметры для решения методом RKF45
double RELERR = 0.0001; //относительная погрешность
double ABSERR = 0.000001; //абсоляютная погрешность
int IFLAG = 1;
double WORK[15];
int IWORK[5];

//Реальное решение
void ActualSolution(double T, double *actualX) {
    actualX[0] =  1 / T; // y = 1 / t
    actualX[1] =  -1 * (1 / (T * T)); // y' = - (1 / T^2)
}

//Восстановление занчений
void reload() {
    T = 1; //начало промежутка
    Tout = 2; //конец промежутка
    X[0] = 1; // y(t=1) = 1
    X[1] = -1;//-1; // y'(t=1) = -1
    IFLAG = 1; //FLAG
    ActualSolution(T, actualX);
}

void Fun(double T, double *X, double *DX) {
    DX[0] = X[1];  // p = y'
    DX[1] = (-1 * (3 * T + 2) * DX[0] - X[0]) / (T * T + T); //p' =
}

//Решение + вывод результатов с помощью подпрограммы RKF45
void RKF()  {
    printf("\nРешение через программу RKF45\n");
    printf("RELERR = %.5f \nABSERR = %.5f\n", RELERR, ABSERR);
    double actualX[2];
    double rightBorder = Tout;
    Tout = T + h;
    while (T < rightBorder) {
        RKF45(Fun, n, X, &T, &Tout, &RELERR, &ABSERR, &IFLAG, WORK, IWORK);
        printf("Tout = %.1f\n", Tout);
        ActualSolution(T, actualX);
        for (int i = 0; i < n; i++) {
            printf("X[%d] = %.8f\t(%.16f)\n", i, X[i], X[i] - actualX[i]);
        }
        Tout += h;
    }
}

//Метод Эйлера-Коши
void CauchyEuler(bool RELERR_FLAG) {
    double curX[n];
    double intermediateX[n];
    double DX[n];
    double intermediateDX[n];

    for (; T < Tout; T += h) {
        for (int i = 0; i < n; i++) {
            curX[i] = X[i];
        }
        if (RELERR_FLAG) {
            ActualSolution(T, curX);
        }
        Fun(T, curX, DX);
        for (int i = 0; i < n; i++) {
            intermediateX[i] = 0;
            intermediateDX[i] += curX[i] + h * DX[i];
        }

        Fun(T + h, intermediateX, intermediateDX);
        for (int i = 0; i < n; i++) {
            X[i] = 0;
            X[i] += (curX[i] + h / 2 * (DX[i] + intermediateDX[i]));
        }

        if (RELERR_FLAG) {
            printf("Tn = %.4f\n", T);
            for (int i = 0; i < n; i++) {
                printf("Xn[%d] = %.8f\n", i, curX[i]);
            }
            printf("Tn+1 = %.4f\t   (полученное - настоящее)\n", T + h);
            ActualSolution(T + h, curX);
            for (int i = 0; i < n; i++) {
                printf("Xn+1[%d] = %.8f\t(%.8f)\n", i, X[i], X[i] - curX[i]);
            }
        }
    }
}


    int main() {
        SetConsoleOutputCP(CP_UTF8);

        //Решение с помощью подпрограммы RKF45
        RKF();

        //Решение RKF45
        reload();
        RKF45(Fun, n, X, &T, &Tout, &RELERR, &ABSERR, &IFLAG, WORK, IWORK);
        printf("\nРешение через программу RKF45\n");
        printf("Tout = %.1f\nRELERR = %.5f\nABSERR = %.5f\n", Tout, RELERR, ABSERR);
        ActualSolution(Tout, actualX);
        for (int i = 0; i < n; i++) {
            printf("X[%d] = %.8f\t(%.12f)\n", i, X[i], X[i] - actualX[i]);
        }

        //Решение методом Эйлера-Коши
        reload();
        CauchyEuler(false);
        printf("Решение методом Эйлера-Коши\n");
        printf("Tout = %.1f\nh = %.4f\n", Tout, h);
        for (int i = 0; i < n; i++) {
            printf("X[%d] = %.8f\t(%.8f)\n", i, X[i], X[i] - actualX[i]);
        }


        reload();
        printf("\nh=%.4f\t\n", h);
        CauchyEuler(true);

        return 0;
    }