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

class Stroke {
public:
    void AddPoint(const Point &p, float t) {
        pts.push_back(p);
        times.push_back(t);
    }
    void Normalize() {
        width = NormalizePointSet(pts, centroid);
        duration = NormalizeFloatSet(times);
    }
    void Draw() {
        SetBrightColor(.5);
        for (Point &p : pts) Circle(p, .001);
    }
    vector<Point> pts;
    vector<float> times;
    float width, duration;
    Point centroid;
};

class Glyph {
public:
    void MouseDown() { strokes.push_back(Stroke()); }
    void MouseUp() { strokes.back().Normalize(); }
    void AddPoint(const Point &p, float t) { strokes.back().AddPoint(p,t); }
    void Draw() { for (Stroke &s : strokes) s.Draw(); }
    vector<Stroke> strokes;
};


#endif
