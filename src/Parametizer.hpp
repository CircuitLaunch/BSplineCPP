//
//  Parametizer.hpp
//  BSpline
//
//  Created by Edward Janne on 5/29/21.
//

#ifndef Parametizer_hpp
#define Parametizer_hpp

#include <stdio.h>
#include <vector>
#include <math.h>

#include "BSpline.hpp"
#include "Functor.hpp"

using namespace std;

class Parametizer
{
    public:
        Parametizer(BSpline &iSpline)
        : spline(iSpline), length(0), spanLengths()
        { }
        
        void init();
        
        float arcLength(float t);
        float timeForArc(float iArc);
        
        float segmentArc(int iSeg, float t);
        float segmentArcDeriv(int iSeg, float t);
        
        float timeForSegmentArc(int iSeg, float iArg);
        
        vector<float> parametizeLinear(int iCount);
        vector<float> parametizeSigmoidal(int iCount);
        
    public:
        class MagFunctor: public Functor
        {
            public:
                MagFunctor(Parametizer &iParametizer) : p(iParametizer) { }
                
                virtual float operator()(float t) { float buff[p.spline.stride]; p.spline.eval(t, buff); float mag = 0.0; int i = p.spline.stride; while(i--) mag += buff[i] * buff[i]; return sqrt(mag); }
                
            protected:
                Parametizer &p;
        };
        
        class MagDFunctor: public Functor
        {
            public:
                MagDFunctor(Parametizer &iParametizer) : p(iParametizer) { }
                
                virtual float operator()(float t) { float buff[p.spline.stride]; p.spline.deriv(t, buff); float mag = 0.0; int i = p.spline.stride; while(i--) mag += buff[i] * buff[i]; return sqrt(mag); }
                
            protected:
                Parametizer &p;
        };
        
        class SegArcFunctor: public Functor
        {
            public:
                SegArcFunctor(Parametizer &iParametizer, int iSeg) : p(iParametizer), seg(iSeg) { }
                
                virtual float operator()(float t) { return p.segmentArc(seg, t); }
                
            protected:
                Parametizer &p;
                int seg;
        };
        
        class SegArcDFunctor: public Functor
        {
            public:
                SegArcDFunctor(Parametizer &iParametizer, int iSeg) : p(iParametizer), seg(iSeg) { }
                
                virtual float operator()(float t) { return p.segmentArcDeriv(seg, t); }
                
            protected:
                Parametizer &p;
                int seg;
        };
    
    public:
        BSpline &spline;
        
        double length;
        vector<double> spanLengths;
};

#endif /* Parametizer_hpp */
