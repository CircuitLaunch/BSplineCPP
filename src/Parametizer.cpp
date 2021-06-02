//
//  Parametizer.cpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#include "Parametizer.hpp"

#include "Legendre.hpp"
#include "Newton.hpp"

void Parametizer::init()
{
    length = 0.0;
    int start = spline.order - 1;
    int end = spline.cpCount;

    spanLengths.clear();
    for(int i = start; i < end; i++) {
        float t0 = spline.knots[i];
        float t1 = spline.knots[i+1];
        MagDFunctor d(*this);
        float spanLength = float(legendreIntegrate(64, t0, t1, d));
        spanLengths.push_back(spanLength);
        length += spanLength;
    }
}

float Parametizer::arcLength(float t)
{
    int knotCount = spline.cpCount + spline.order;
    int start = spline.order - 1;
    float arcLen = 0.0;
    int i;
    if(t >= spline.knots[knotCount - 1]) return length;
    for(i = start; i < knotCount && t > spline.knots[i]; i++)
        arcLen += spanLengths[i - start];
    MagDFunctor d(*this);
    return arcLen + float(legendreIntegrate(64, spline.knots[i], t, d));
}

float Parametizer::timeForArc(float iArc)
{
    if(iArc >= length) return spline.knots[spline.cpCount];

    float len = 0.0;

    int i;
    for(i = 0; (len + spanLengths[i]) < iArc; i++)
        len += spanLengths[i];

    return timeForSegmentArc(i, iArc - len);
}

float Parametizer::segmentArc(int iSeg, float t)
{
    float t0 = spline.knots[iSeg + spline.order - 1];

    MagDFunctor d(*this);
    return legendreIntegrate(64, t0, t0 + t, d);
}

float Parametizer::segmentArcDeriv(int iSeg, float t)
{
    float t0 = spline.knots[iSeg + spline.order - 1];

    MagDFunctor d(*this);
    return d(t0 + t);
}

float Parametizer::timeForSegmentArc(int iSeg, float iArc)
{
    int start = spline.order - 1;
    float t0 = spline.knots[iSeg + start];
    float hint = 0.0;

    SegArcFunctor f(*this, iSeg);
    SegArcDFunctor d(*this, iSeg);

    return t0 + newtonSolve(iArc, hint, f, d);
}

vector<float> Parametizer::parametizeLinear(int iCount)
{
    vector<float> times;

    times.push_back(0.0);
    float step = float(length) / float(iCount);
    float arc = step;
    while(iCount--) {
        times.push_back(timeForArc(arc));
        arc += step;
    }

    return times;
}

vector<float> Parametizer::parametizeSigmoidal(int iCount)
{
    vector<float> times;

    times.push_back(0.0);
    float step = 1.0 / float(iCount);
    float x = step;
    while(iCount--) {
        float arc = length / (1 + exp(-16.0 * (x - 0.5)));
        times.push_back(timeForArc(arc));
        x += step;
    }

    return times;
}
