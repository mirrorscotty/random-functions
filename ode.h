#ifndef RUNGEKUTTA2_H
#define RUNGEKUTTA2_H

vector* euler(double (*)(double, double),
        double, double, double, double);
vector* rungekutta2(double (*)(double, double),
        double, double, double, double, double);
vector* rungekutta4(double (*)(double, double),
        double, double, double, double);

#endif

