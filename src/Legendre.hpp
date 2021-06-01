//
//  Legendre.hpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#ifndef Legendre_hpp
#define Legendre_hpp

#include <stdio.h>
#include "Functor.hpp"

double legendreIntegrate(int order, double from, double to, Functor &f);

#endif /* Legendre_hpp */
