#ifndef RUNGEKUTTA2_H
#define RUNGEKUTTA2_H

vector* euler(double (*)(double, double),
        double, double, double, int);
vector* rungekutta2(double (*)(double, double),
        double, double, double, double, int);
vector* rungekutta4(double (*)(double, double),
        double, double, double, int);
matrix* rungekutta4s(double (**)(double, vector*),
        double, vector*, vector*, int);

#endif

