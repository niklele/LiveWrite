//
//  stroke.h
//  LiveWrite
//
//  Created by Ben Mildenhall on 11/8/13.
//  Copyright (c) 2013 Ben Mildenhall. All rights reserved.
//

#ifndef LiveWrite_stroke_h
#define LiveWrite_stroke_h

#include "utilities.h"

class Glyph {
public:
    void MouseDown() { ; }
    void MouseUp() { ; }
    void AddPoint(const Point &p, float t) {
        points.push_back(p);
        times.push_back(t);
    }
    void Draw() {
        for (Point &p : points) Circle(p, .001);
    }
    void Normalize() {
        width = NormalizePointSet(points, centroid);
        duration = NormalizeFloatSet(times);
    }
    void Empty() {
        points.clear();
        times.clear();
    }
    Glyph() { }
    Glyph(const Glyph &g) {
        *this = g;
    }
    Glyph &operator=(const Glyph &g) {
        Empty();
        width = g.width;
        duration = g.duration;
        centroid = g.centroid;
        for (int i = 0; i < g.times.size(); ++i)
            AddPoint(g.points[i], g.times[i]);
        return *this;
    }
    vector<Point> points;
    vector<float> times;
    float width, duration;
    Point centroid;
};


Glyph InterpGlyph(Glyph &g, int num) {
    Glyph ret;
    for (int n = 0; n <= num; ++n) {
        float t = (float)n / num;
        int i = 1;
        for (;i < g.times.size(); ++i) {
            if (g.times[i] >= t) break;
        }
        float interp = (t - g.times[i-1]) / (g.times[i] - g.times[i-1]);
        ret.points.push_back(g.points[i-1] * (1-interp) + g.points[i] * interp);
        ret.times.push_back(t);
    }
    return ret;
}

Glyph InterpSpaceGlyph(Glyph &g, int num) {
    Glyph ret;
    if (!g.points.size()) {
        cout << "No points!!!" << endl;
        return ret;
    }
    Point prev = g.points[0];
    float length = 0;
    for (int i = 1; i < g.points.size(); ++i) {
        Point next = g.points[i];
        length += Dist(prev, next);
        prev = next;
    }
    for (int n = 0; n <= num; ++n) {
        float dist = (float)n / num * length;
        prev = g.points[0];
        float acc = 0, accprev = 0;
        int i = 1;
        for (; i < g.times.size(); ++i) {
            Point next = g.points[i];
            acc += Dist(prev, next);
            if (acc >= dist) break;
            accprev = acc;
            prev = next;
        }
        float interp = (dist - accprev) / (acc - accprev);
        ret.points.push_back(g.points[i-1] * (1-interp) + g.points[i] * interp);
        ret.times.push_back((float)n / num);
    }
    return ret;
}

Glyph ReduceGlyph(Glyph &g, int num, float stdev) {
    Glyph ret;
    for (int n = 0; n <= num; ++n) {
        float t = (float)n / num;
        Point wavg;
        float w = 0;
        for (int i = 0; i < g.times.size(); ++i) {
            float weight = NormalPDF(g.times[i], t, stdev);
            wavg += g.points[i] * weight;
            w += weight;
        }
        ret.points.push_back(wavg / w);
        ret.times.push_back(t);
    }
    return ret;
}

Glyph SmoothGlyph(Glyph &g, int supern, int reducen, float std) {
    Glyph ret = InterpGlyph(g, supern);
    ret = ReduceGlyph(ret, reducen, std);
    ret.Normalize();
    return ret;
}


#endif
