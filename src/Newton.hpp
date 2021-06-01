//
//  Newton.hpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#ifndef Newton_hpp
#define Newton_hpp

#include <stdio.h>
#include <math.h>
#include "Functor.hpp"

double newtonSolve(double tgt, double hint, Functor &f, Functor &d, int maxSteps = 100, double tol = 0.0001, double epsilon = 0.0001);

#endif /* Newton_hpp */
