#include <iostream>
#include "Function.h"
#include "iomanip"


int main() {
    Function spline = Function(); //Построение сплайн-функции

    double
    a = 0.0, //Начало промежутка
    b = 1.0, //Конец промежутка
    e = 0.000001; //Погрешность
    double bis = Function().bisectionSpline(spline, a, b, e); //Результат метода бисекции
    double bisl = Function().bisectionLagrange(spline, a, b, e); //Результат метода бисекции
    //Выведение результатов в консоль
    std::cout << "\nResult bisection = " << std::setprecision(15) << bis << std::endl;
    std::cout << "\nResult bisection = " << std::setprecision(15) << bisl << std::endl;
    return 0;

}




