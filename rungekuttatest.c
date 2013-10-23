#include "matrix.h"
#include "ode.h"
#include <math.h>
#include <stdio.h>

/* This file is designed to test the rungekutta2 function against the example
 * provided on the wikipedia page.
 */

double f(double t, double y)
{
    return tan(y) + 1;
}

int main(int argc, char *argv[])
{
    vector *y;
    printf("Euler's Method:\n");
    y = euler(&f, 0.025, 1, 1, 4);
    PrintVector(y);
    DestroyVector(y);
    printf("Runge-Kutta Method:\n");
    y = rungekutta2(&f, .025, 2./3., 1, 1, 4);
    PrintVector(y);
    DestroyVector(y);
    return 0;
}

