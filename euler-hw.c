#include "matrix.h"
#include "ode.h"
#include <stdio.h>
#include <math.h>

vector* exactsoln(double (*f)(double), double x0, double h, int tn)
{
    int i;
    vector *y;
    y = CreateVector(tn);
    for(i=0; i<tn; i++)
        setvalV(y, i, f(x0+i*h));
    return y;
}

double p1(double x, double y)
{
    return pow(x, 1./3.);
}

double p1exact(double x)
{
    return 3./4.*pow(x, 3./4.) + 2;
}

double p2(double x, double y)
{
    return y-x*sin(x);
}

double p2exact(double x)
{
    return .5*(exp(-x)+x*sin(x)-x*cos(x)+cos(x));
}

int main(int argc, char *argv)
{
    double h = .01;
    double x0;
    double x1;
    int xn;
    int i;
    vector *y1, *y1e, *y1rk;
    vector *y2, *y2e, *y2rk;
    matrix *output;

    /* Problem 1 */
    x0 = 0;
    x1 = 1;
    xn = (x1-x0)/h;
    y1 = exactsoln(&p1exact, x0, h, xn);
    y1e = euler(&p1, h, 2, x0, xn);
    y1rk = rungekutta4(&p1, h, 2, x0, xn);

    /* Problem 2 */
    x0 = 0;
    x1 = 5;
    xn = (x1-x0)/h;
    y2 = exactsoln(&p2exact, x0, h, xn);
    y2e = euler(&p2, h, 1, x0, xn);
    y2rk = rungekutta4(&p2, h, 1, x0, xn);

    /* Output */
    output = CatColVector(6, y1, y1e, y1rk, y2, y2e, y2rk);
    //mtxprnt(output);
    //PrintVector(y1);

    /* Cleanup */
    DestroyVector(y1);
    DestroyVector(y2);
    DestroyVector(y1e);
    DestroyVector(y2e);
    DestroyVector(y1rk);
    DestroyVector(y2rk);
    DestroyMatrix(output);

    return 0;
}


