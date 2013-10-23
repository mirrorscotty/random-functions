#include "matrix.h"
#include "ode.h"

/* Use Euler's Method to solve an ODE */
vector* euler(double (*f)(double, double), /* y' = f(t,y) */
              double h, /* Step size */
              double y0, /* Initial value */
              double t0, /* Initial time */
              double tn) /* Number of time steps */
{
    int i; /* Loop index */
    vector *y; /* Place to store all of the output variables */
    double yi, yi1; /* Current and next values for y */
    double t; /* Variable to store the time in */

    y = CreateVector(tn+1); /* Initialize the output vector */
    yi = y0; /* Set the current value to the initial value */

    setvalV(y, 0, y0); /* Save the initial value to the output vector */
    
    for(i=1; i<tn+1; i++) {
        t = t0 + h*i; /* Calculate the current time */
        yi1 = yi + h*f(t,yi);/*Calculate the next y value */
        setvalV(y, i, yi1); /* Save the value we just calculated */
        yi=yi1; /* Get ready for the next loop iteration. */
    }
    
    return y;
}

/* Integrate an ODE using the second order Runge-Kutta Method */
vector* rungekutta2(double (*f)(double, double), /* y' = f(t,y) */
                    double h, /* time step */
                    double a, /* alpha value */
                    double y0, /* Initial value */
                    double t0, /* Initial time */
                    double tn) /* Number of time steps */
{
    int i; /* Loop index */
    vector *y; /* Place to store all of the output variables */
    double yi, yi1; /* Current and next values for y */
    double k1, k2; /* Intermediate values */
    double t; /* Variable to store the time in */

    y = CreateVector(tn+1); /* Initialize the output vector */
    yi = y0; /* Set the current value to the initial value */

    setvalV(y, 0, y0); /* Save the initial value to the output vector */
    
    for(i=1; i<tn+1; i++) {
        t = t0 + h*i; /* Calculate the current time */
        k1 = f(t, yi); /* Set the first intermediate value */
        k2 = f(t+a*h, yi+a*h*k1); /* Set the second one */
        yi1 = yi + h*((1-1/(2*a))*k1+1/(2*a)*k2);/*Calculate the next y value */
        setvalV(y, i, yi1); /* Save the value we just calculated */
        yi=yi1; /* Get ready for the next loop iteration. */
    }
    
    return y;
}

