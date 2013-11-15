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
    vector<Point> points;
    vector<float> times;
    float width, duration;
    Point centroid;
};


#endif
