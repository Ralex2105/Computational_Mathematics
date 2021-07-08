#ifndef SPLINE_FUNCTION_H
#define SPLINE_FUNCTION_H


#include <iomanip>

class Function {

    // Переменные для построения функции Spline
public:
    double bS[7], cS[7], dS[7];
    //Таблично заданная функция
    double
            y[7] = {0, 1.0000, 1.2214, 1.4918, 2.0138, 2.4596, 2.7183};
    double
            x[7] = {0, 0.0, 0.2, 0.4, 0.7, 0.9, 1.0};

    double
            y1[7] = {1.0000, 1.2214, 1.4918, 2.0138, 2.4596, 2.7183};
    double
            x1[7] = {0.0, 0.2, 0.4, 0.7, 0.9, 1.0};

    int n = 6; //степень полинома

    //Массив значений для графика
    double arrayG[41] = {0.0, 0.025, 0.05, 0.75,
                         0.1, 0.125, 0.15, 0.175,
                         0.2, 0.225, 0.25, 0.275,
                         0.3, 0.325, 0.35, 0.375,
                         0.4, 0.425, 0.45, 0.475,
                         0.5, 0.525, 0.55, 0.575,
                         0.6, 0.625, 0.65, 0.675,
                         0.7, 0.725, 0.75, 0.775,
                         0.8, 0.825, 0.85, 0.875,
                         0.9, 0.925, 0.95, 0.975, 1.0
    };

    //Сплайн
    //Построение функции, используя Spline
    Function() {
        spline(n, x, y, bS, cS, dS);
    };
    //Стандартная подпрограмма Spline
    void spline(int n,double * x,double * y,double * b,double * c,double * d)
    {
        int i,ib,nm1;
        double t ;
        nm1=n-1;
        if (n < 2)  return;
        if (n < 3)  goto l20;
        d[1]=x[2]-x[1];
        c[2]=(y[2]-y[1])/d[1];
        for (i=2 ;i<=nm1;i++)
        {
            d[i]=x[i+1]-x[i];
            b[i]=2*(d[i-1]+d[i]);
            c[i+1]=(y[i+1]-y[i])/d[i];
            c[i]=c[i+1]-c[i];
        }
        b[1]=-d[1];
        b[n]=-d[n-1];
        c[1]=0;
        c[n]=0;
        if (n == 3) goto l10;
        c[1]=c[3]/(x[4]-x[2])-c[2]/(x[3]-x[1]);
        c[n]=c[n-1]/(x[n]-x[n-2])-c[n-2]/(x[n-1]-x[n-3]);
        c[1]=c[1]*sqrt(d[1])/(x[4]-x[1]);
        c[n]=-c[n]*sqrt(d[n-1])/(x[n]-x[n-3]);
        l10:    for (i=2;i<=n;i++)
    {
        t=d[i-1]/b[i-1];
        b[i]=b[i]-t*d[i-1];
        c[i]=c[i]-t*c[i-1];
    }
        c[n]=c[n]/b[n];
        for (ib=1;ib<=nm1;ib++)
        {
            i=n-ib;
            c[i]=(c[i]-d[i]*c[i+1])/b[i];
        }
        b[n]=(y[n]-y[nm1])/d[nm1]+d[nm1]*(c[nm1]+2*c[n]);
        for (i=1;i<=nm1;i++)
        {
            b[i]=(y[i+1]-y[i])/d[i]-d[i]*(c[i+1]+2*c[i]);
            d[i]=(c[i+1]-c[i])/d[i];
            c[i]=3*c[i];
        }
        c[n]=3*c[n];
        d[n]=d[n-1];
        return;
        l20:  b[1]=(y[2]-y[1])/(x[2]-x[1]);
        c[1]=0;
        d[1]=0;
        b[2]=b[1];
        c[2]=0;
        d[2]=0;
        l30:
        return;
    }
    //Стандартная подпрограмма Seval
    static double seval(int n, double u, double * x, double * y, double * b, const double * c, double * d)
    {
        int i,j,k;
        double dx;
        i=1;
        if (i >= n) i=1;
        if (u < x[i]) goto l10;
        if (u <= x[i+1]) goto l30;
        l10: i=1;
        j=n+1;
        l20: k=(i+j) / 2;
        if (u < x[k]) j=k;
        if (u >= x[k]) i=k;
        if (j > (i+1)) goto l20;
        l30: dx=u-x[i];
        return y[i]+dx*(b[i]+dx*(c[i]+dx*d[i]));
    }
    //Получение значения функции при помощи Seval
    double getYSpline(double u) {
        return seval(n, u, x, y, bS, cS, dS) + u * u - 2;
    }
    //Метод бисекции для сплайн
    double bisectionSpline(Function function, double a, double d, double e) {
        double aX, bX, cX, // [A, B] - интервал, С - середнина
        aF, bF, cF; // Значение spline-функции в этих точках

        aX = 0.0; //Начало промежутка
        bX = 1.0; //Конец промежутка
        cX = 0.0; //Середина промежутка

        aF = Function().getYSpline(aX); //Значение spline-функции в точке a
        bF = Function().getYSpline(bX); //Значение spline-функции в точке b
        int i = 0;
        while (i < 300000) {
            i++;

            if (aF * bF > 0) {
                //Сообщение, если знаки на концах будут одинковыми
                printf("Value at the ends [%lf, %lf] have the same sign\n", aX, bX);
                break;
            }

            cX = (aX + bX) / 2; //Середина промежутка
            cF = Function().getYSpline(cX);
            //Присвоение точки середины одному из концов промежутка (середина становится концом)
            if (aF * cF < 0) {
                bF = cF;
                bX = cX;
            } else {
                aF = cF;
                aX = cX;
            }
            //Завершение работы, если данные удовлетворяют погрешности
            if (std::abs(aX - bX) < e) {
                break;
            } else {

            }
            /*
            //Вывод подробных результатов в терминал
            printf("Bisection Spline, interval [%lf, %lf]:\n", aX, bX);
            std::cout << "F(a)=" << std::setprecision(10) << aF << std::endl;
            std::cout << "F(b)=" << std::setprecision(10) << bF << std::endl;
            std::cout << "F(c)=" << std::setprecision(10) << cF << std::endl;
            printf("\n");
             */

        }
        /*
        // Подробные результаты решения
        printf("Need %i iterations with precision %lf\n", i, e);
        printf("Value of bisection %lf\n\n\n", cX);
         */


        for (int f = 0; f < 41; f++) {
            double cXG = Function().getYSpline(arrayG[f]);
            printf("%lf\n", cXG);
        }
        return cX;
    }



    //Лагранж
    //Стандартная подпрограмма Lagrange
    double lagrange(int n, double * x1, double * y1, double xc)
    {
        double Ch;
        double Zn;
        int k;
        double R=0;
        for (int i = 0; i < n; i++) {
            Ch = 1; Zn = 1;
            for (k = 1; k < n; k++ ) {
                if ( k == i ) continue;
                Ch *= xc - x[k];
            }
            for(k= 1; k < n;k++) {
                if (x[i] == x[k]) continue;
                Zn *= x[i] - x[k];
            }
            R += y[i]*Ch/Zn;
        }
        return R;
    }
    //Получение значения в точке
    double getYL(double u) {
           return lagrange(7, x1, y1, u) + u * u - 2;
       }
    //Метод бисекции для лагранж
    double bisectionLagrange(Function function, double a, double d, double e) {
        printf("\n\n\n");
        double aX, bX, cX, // [A, B] - интервал, С - середнина
        aF, bF, cF; // Значение spline-функции в этих точках

        aX = 0.0; //Начало промежутка
        bX = 1.0; //Конец промежутка
        cX = 0.0; //Середина промежутка

        aF = Function().getYL(aX); //Значение spline-функции в точке a
        bF = Function().getYL(bX); //Значение spline-функции в точке b
        int i = 0;
        while (i < 300000) {
            i++;
            cX = (aX + bX) / 2; //Середина промежутка
            cF = Function().getYL(cX);

            //Присвоение точки середины одному из концов промежутка (середина становится концом)
            if (aF * cF < 0) {
                bF = cF;
                bX = cX;
            } else {
                aF = cF;
                aX = cX;
            }
            //Завершение работы, если данные удовлетворяют погрешности
            if (std::abs(aX - bX) < e) {
                break;
            } else {

            }
            /*
            //Вывод подробных результатов в терминал
            printf("Bisection Lagrange, interval [%lf, %lf]:\n", aX, bX);
            std::cout << "F(a)=" << std::setprecision(10) << aF << std::endl;
            std::cout << "F(b)=" << std::setprecision(10) << bF << std::endl;
            std::cout << "F(c)=" << std::setprecision(10) << cF << std::endl;
            printf("\n");
             */

        }
        /*
        // Подробные результаты решения
        printf("Need %i iterations with precision %lf\n", i, e);
        printf("Value of bisection %lf\n\n\n", cX);
         */



        for (int f = 0; f < 41; f++) {

            double cXG = Function().lagrange(7, x1, y1, arrayG[f]) + arrayG[f] * arrayG[f] - 2;
            printf("%lf\n", cXG);
        }


        return cX;
    }
};
#endif