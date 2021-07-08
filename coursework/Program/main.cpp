#include <iostream>
#include "quanc8.h"
#include "rkf45.h"
#include <cmath>


//Функиця e^2 для метода calcL
double funL(double x) {
    return exp(x);                //Подинтегральная функция
}

//Вычисление параметра L (начальной длины пружины)
double calcL() {
    double result;                //Переменная для результата
    double errest, flag;          //Параметры для quanc8
    int nofun;                    //Параметр для quanc8
    quanc8(funL, 0.0, 1.0, 1e-14, 1e-14, &result, &errest, &nofun, &flag);
                                  //Вычиление интеграла с помощью quanc8
    result *= 0.5836896;          //Вычисления результата
    return result;
}

//ПАРАМЕТРЫ:
const double g = 9.81;             //Гравитационная постоянная
const double M = 1;                //Масса маятника
double K;                          //Исследуемый параметр K
double L = calcL();                //Начальная длина пружины L

//Параметры для решения
const int n = 4;                   //Количество уравнений в системе
double X[] = {0.0, 0.0, 0.0, 4.0}; //массив из двух уравнений
double actualX[7] = {0.0, 0.303, -0.465, 0.592, -0.409, 0.164, 0.180};
                                   //Экспериментальные значения
double T = 0.0;                    //Начала интервала
double Tout = 2.4;                 //Конец интервала
double h = 0.4;                    //Шаг

//Параметры для решения методом RKF45
double RELERR = 0.000000001;             //относительная погрешность
double ABSERR = 0.000000001;             //абсоляютная погрешность
int IFLAG = 1;                     //
double WORK[60];                   // - Параметры для RKF45
int IWORK[20];                     //

//Система дифференциальных уравнений
void equations(double t, double *y, double *dy) {
    dy[0] = y[1];
    dy[1] = -K * y[0] / M - g * (1 - cos(y[2])) + (L + y[0]) * y[3] * y[3];
    dy[2] = y[3];
    dy[3] = -g * sin(y[2]) / (L + y[0]) - 2 * y[1] * y[3] / (L + y[0]);
}

//RKF45 для нахождения минимального значения K
double RKF(double Parameter)  {
    T = 0.0;
    Tout = 2.4;
    double sum = 0.0;
    K = Parameter;
    double X[] = {0, 0, 0, 4};
    double rightBorder = Tout;
    Tout = T;
    for (int i = 0; Tout <= rightBorder; i++) {
        RKF45(equations, n, X, &T, &Tout, &RELERR, &ABSERR, &IFLAG, WORK, IWORK);
        Tout += h;
        sum += (X[0] - actualX[i]) * (X[0] - actualX[i]);
    }
    IWORK[60];
    IFLAG = 1;
    WORK[20];
    return sum;
}

//RKF45 для вычисления точек и разницы между экспериментальными значениями
double RKF_end(double Parameter)  {
    T = 0.0;
    Tout = 2.4;
    double sum = 0.0;
    K = Parameter;
    double X[] = {0, 0, 0, 4};
    double rightBorder = Tout;
    Tout = T;
    printf("\nTout\tCalculated\t Real\t  EPS\n");
    for (int i = 0; Tout <= rightBorder; i++) {
        RKF45(equations, n, X, &T, &Tout, &RELERR, &ABSERR, &IFLAG, WORK, IWORK);
            printf("%.1f\t%.6f\t%.3f\t%.6f\n", Tout, X[0], actualX[i], X[0] - actualX[i]);
        Tout += h;
        sum += (X[0] - actualX[i]) * (X[0] - actualX[i]);
    }
    IWORK[60];
    IFLAG = 1;
    WORK[20];
    return sum;
}

//Поиск минимума функции
typedef double(*Func)(double x);
double Minimum(double Start, double End, Func f, double Step) {
    double xMin = Start;
    double fMin = f(Start);
    while(Start < End) {
        if(f(Start) < fMin) {
            xMin = Start;
            fMin = f(Start);
        }
        Start +=Step;
    }
    return xMin;
}

int main() {
    printf("Calculated parameter: L = %.4f\n", L);
    K = Minimum(36, 46, RKF, 0.001);
    printf("Calculated parameter K = %.3f\n", K);
    RKF_end(K);
    return 0;
}