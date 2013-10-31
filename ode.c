#include "matrix.h"
#include "ode.h"

/* Use Euler's Method to solve an ODE */
vector* euler(double (*f)(double, double), /* y' = f(t,y) */
              double h, /* Step size */
              double y0, /* Initial value */
              double t0, /* Initial time */
              int tn) /* Number of time steps */
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
                    int tn) /* Number of time steps */
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

/* Solve a first order ODE using a fourth order Runge-Kutta method */
vector* rungekutta4(double (*f)(double, double), /* y' = f(t,y) */
              double h, /* Step size */
              double y0, /* Initial value */
              double t0, /* Initial time */
              int tn) /* Number of time steps */
{
    int i; /* Loop index */
    vector *y; /* Place to store all of the output variables */
    double yi, yi1; /* Current and next values for y */
    double k1, k2, k3, k4; /* Intermediate values */
    double t; /* Variable to store the time in */

    y = CreateVector(tn+1); /* Initialize the output vector */
    yi = y0; /* Set the current value to the initial value */

    setvalV(y, 0, y0); /* Save the initial value to the output vector */
    
    for(i=1; i<tn+1; i++) {
        t = t0 + h*i; /* Calculate the current time */
        k1 = h*f(t,yi);
        k2 = h*f(t+.5*h, yi+.5*k1);
        k3 = h*f(t+.5*h, yi+.5*k2);
        k4 = h*f(t+h, yi+k3);
        yi1 = yi + 1./6.*k1 + 1./3.*k2 + 1./3.*k3 + 1./6.*k4; /* + O(h^5) */
        setvalV(y, i, yi1); /* Save the value we just calculated */
        yi=yi1; /* Get ready for the next loop iteration. */
    }
    
    return y;
}

/* Solve a first order system of ODEs using a fourth order Runge-Kutta method */
matrix* rungekutta4s(double (**f)(double, vector*), /* y' = f(t,y) */
              double h, /* Step size */
              vector *y0, /* Initial value */
              vector *t0, /* Initial time */
              int tn) /* Number of time steps */
{
    int i = 0; /* Loop index */
    int j;
    matrix *y; /* Place to store all of the output variables */
    vector *yi; /* Vector containing all the current values */
    vector *tmp;
    double k1, k2, k3, k4; /* Intermediate values */
    vector *t; /* Variable to store the time in */
    int neq = len(y0); /* Number of equations to solve. */

    y = CreateMatrix(tn+1, neq); /* Initialize the output matrix */
    for(j=0; j<neq; j++)
        setval(y, valV(y0,j), i, j); /* Save the initial values to the output matrix */
    
    for(i=0; i<tn; i++) {
        yi = ExtractRowAsVector(y, i);

        for(j=0; j<neq; j++) {
            t = scalaraddV(h*(i+1), t0); /* Calculate the current time */

            k1 = h*(f[j])(valV(t, j),yi);

            tmp = scalaraddV(.5*k1, yi);
            k2 = h*(f[j])(valV(t, j)+.5*h, tmp);
            DestroyVector(tmp);

            tmp = scalaraddV(.5*k2, yi);
            k3 = h*(f[j])(valV(t, j)+.5*h, tmp);
            DestroyVector(tmp);

            tmp = scalaraddV(k3, yi);
            k4 = h*(f[j])(valV(t, j)+h, tmp);
            DestroyVector(tmp);

            setval(y,
                   valV(yi,j) + 1./6.*k1 + 1./3.*k2 + 1./3.*k3 + 1./6.*k4,
                   i+1, j); /* + O(h^5) */

            DestroyVector(t);
        }
        DestroyVector(yi);
    }
    
    return y;
}

