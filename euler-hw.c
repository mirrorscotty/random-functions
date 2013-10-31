#include "matrix.h"
#include "ode.h"
#include <stdio.h>
#include <stdlib.h>
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
    return y-x*sin(x);
}

double p1exact(double x)
{
    return .5*(exp(-x)+x*sin(x)-x*cos(x)+cos(x));
}

double p2z1(double x, vector *z)
{
    return valV(z,1);
}

double p2z2(double x, vector *z)
{
    return 6*exp(x) - 2*valV(z, 0) - 3*valV(z, 1);
}

double p2exact(double x)
{
     return exp(-2*x)*(4*exp(x) + exp(3*x) - 3);
}

int main(int argc, char *argv)
{
    double h = .01;
    double x0;
    double x1;
    int xn;
    int i;

    vector *y1, *y1e, *y1rk;

    vector *y2;
    matrix *y2rk;
    vector *x02, *y02;

    double (**z)(double, vector*);
    z = (double (**)(double, vector*))
        calloc(sizeof(double (**)(double, vector*)), 2);


    /* Problem 1 */
    x0 = 0;
    x1 = 5;
    xn = (x1-x0)/h;
    y1 = exactsoln(&p1exact, x0, h, xn);
    y1e = euler(&p1, h, 2, x0, xn);
    y1rk = rungekutta4(&p1, h, 2, x0, xn);


    /* Problem 2 */
    x0 = 0;
    x1 = 1;
    xn = (x1-x0)/h;

    x02 = CreateVector(2);
    y02 = CreateVector(2);
    setvalV(y02, 0, 2);
    setvalV(y02, 1, 3);

    z[0] = &p2z1;
    z[1] = &p2z2;

    y2 = exactsoln(&p2exact, x0, h, xn);
    y2rk = rungekutta4s(z, h, y02, x02, xn);

    /* Output */
    //PrintVector(y1);
    //PrintVector(y1e);
    //PrintVector(y1rk));
    PrintVector(y2);
    mtxprnt(y2rk);

    /* Cleanup */
    DestroyVector(y1);
    DestroyVector(y2);
    DestroyVector(y1e);
    DestroyVector(y1rk);
    DestroyMatrix(y2rk);
    free(z);

    return 0;
}


