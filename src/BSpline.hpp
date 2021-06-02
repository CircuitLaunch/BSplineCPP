//
//  BSpline.hpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#ifndef BSpline_hpp
#define BSpline_hpp

#include <stdio.h>

class BSpline
{
    public:
        float *cpBuffer;
        float *knots;
        int stride;
        int order;
        int cpCount;

    public:
        BSpline(float *iCPBuffer, float *iKnotBuffer, int iOrder = 4);

        virtual void init(int iStride, int iCPCount);

        float basis(int i, int k, float t);

        void eval(float t, float *oPoint);
        void deriv(float t, float *oPoint);
};

#endif /* BSpline_hpp */
