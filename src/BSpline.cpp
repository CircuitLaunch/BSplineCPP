//
//  BSpline.cpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#include "BSpline.hpp"

BSpline::BSpline(float *iCPBuffer, float *iKnotBuffer, int iOrder)
: cpBuffer(iCPBuffer), knots(iKnotBuffer), stride(0), order(iOrder), cpCount(0)
{ }

void BSpline::init(int iStride, int iCPCount)
{
    stride = iStride;
    cpCount = iCPCount;
    int knotCount = cpCount + order;

    float k = 1.0;
    int i = order;
    while(i--) knots[i] = 0.0;
    for(i = order; i < cpCount; i++, k++)
        knots[i] = k;
    for(i = cpCount; i < knotCount; i++)
        knots[i] = k;
}

float BSpline::basis(int i, int k, float t)
{
    if(!k) return ((knots[i] <= t) && (t <= knots[i+1])) ? 1.0 : 0.0;

    float n0 = t - knots[i];
    float d0 = knots[i + k] - knots[i];
    float b0 = basis(i, k - 1, t);

    float n1 = knots[i + k + 1] - t;
    float d1 = knots[i + k + 1] - knots[i + 1];
    float b1 = basis(i + 1, k - 1, t);

    float left = ((b0 != 0.0) && (d0 != 0.0)) ? n0 * b0 / d0 : 0.0;
    float right = ((b1 != 0.0) && (d1 != 0.0)) ? n1 * b1 / d1 : 0.0;
    return left + right;
}

void BSpline::eval(float t, float *oPoint)
{
    int knotCount = cpCount + order;
    if(t < 0.0) t = 0.0;

    if(t > knots[knotCount - 1])
        t = knots[knotCount - 1];

    for(int i = 0; i < stride; i++)
        oPoint[i] = 0.0;

    for(int j = 0; j < cpCount; j++) {
        float b = basis(j, order - 1, t);
        int offset = j * stride;
        for(int i = 0; i < stride; i++)
            oPoint[i] += cpBuffer[offset + i] * b;
    }
}

void BSpline::deriv(float t, float *oPoint)
{
    int knotCount = cpCount + order;
    if(t < 0.0) t = 0.0;

    if(t > knots[knotCount - 1])
        t = knots[knotCount - 1];

    int n = order - 1;

    for(int i = 0; i < stride; i++)
        oPoint[i] = 0.0;

    for(int j = 0; j < (cpCount - 1); j++) {
        float u0 = knots[j + n + 1];
        float u1 = knots[j + 1];
        float fn = (float(n) / (u0 - u1)) * basis(j + 1, n - 1, t);
        int offset = j * stride;
        for(int i = 0; i < stride; i++)
            oPoint[i] += (cpBuffer[offset + stride + i] - cpBuffer[offset + i]) * fn;
    }
}
