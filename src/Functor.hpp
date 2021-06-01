//
//  Functor.hpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#ifndef Functor_hpp
#define Functor_hpp

#include <stdio.h>

class Functor
{
    public:
        virtual float operator()(float t) = 0;
};

#endif /* Functor_hpp */
