//
//  Newton.cpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#include "Newton.hpp"
#include "Functor.hpp"

double newtonSolve(double tgt, double hint, Functor &f, Functor &d, int maxSteps, double tol, double epsilon)
{
    double x0 = hint;
    double x1 = hint;
    int deltaCount = 0;
    double lastDelta = 0.0;
    while(maxSteps--) {
        double yp = double(d(x0));
        if(fabs(yp) < epsilon) return x1;
        double newTry = double(f(x0));
        x1 = x0 - (newTry - tgt) / yp;
        if(fabs(x1-x0)/fabs(x1) < tol) return x1;
        double delta = fabs(x0 - x1);
        if(delta == lastDelta) {
            if(deltaCount > 10) return (x0 + x1) * 0.5;
            deltaCount++;
        }
        lastDelta = delta;
        x0 = x1;
    }
    return x1;
}
