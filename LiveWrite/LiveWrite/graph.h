//
//  graph.h
//  LiveWrite
//
//  Created by Ben Mildenhall on 12/1/13.
//  Copyright (c) 2013 Ben Mildenhall. All rights reserved.
//

#ifndef LiveWrite_graph_h
#define LiveWrite_graph_h

#include "utilities.h"
#include <vector>
#include <iostream>
using namespace std;

class Edge {
public:
    Edge() { }
    Edge(int _a, int _b, float _w)
    : a(_a), b(_b), w(_w) { }
    Edge(const Edge &e)
    : a(e.a), b(e.b), w(e.w) { }
    Edge &operator=(const Edge &e) {
        a = e.a;
        b = e.b;
        w = e.w;
        return *this;
    }
    bool operator<(const Edge &e) const { return w < e.w; }
    
    int a, b;
    float w;
};

float Kruskal(int n, vector<Edge> &edges, vector<Edge> &mst) {
    sort(edges.begin(), edges.end());
    vector<vector<int>> unions;
    for (int i = 0; i < n; ++i) {
        vector<int> v;
        v.push_back(i);
        unions.push_back(v);
    }
    int eind = 0;
    float maxw;
    while (unions.size() > 1 && eind < edges.size()) {
        Edge &e = edges[eind++];
        //cout << "Dealing with edge " << e.a << " to " << e.b << " weight " << e.w << endl;
        int au = -1, bu = -1;
        for (int ui = 0; ui < unions.size() && (au == -1 || bu == -1); ++ui) {
            vector<int> &u = unions[ui];
            if (au == -1 && find(u.begin(), u.end(), e.a) != u.end())
                au = ui;
            if (bu == -1 && find(u.begin(), u.end(), e.b) != u.end())
                bu = ui;
        }
        //cout << "Found in clusters " << au << " and " << bu << endl;
        if (au != bu) {
           // cout << "Merging: " << unions[au] << endl;
           // cout << " and: " << unions[bu] << endl;
            unions[au].insert(unions[au].end(), unions[bu].begin(), unions[bu].end());
            unions.erase(unions.begin() + bu);
            mst.push_back(e);
            maxw = e.w;
        }
    }
    cout << "Connected " << n << " vertices with ";
    cout << mst.size() << " edges" << endl;
    return maxw;
}

void DrawGraph(vector<Point> &pts, vector<Edge> &edges, vector<string> &labels, float cscale) {
    Color(1,1,1);
    for (int i = 0; i < pts.size(); ++i)
        WriteStr(labels[i], pts[i]);
    bool correct = true;
    for (int ei = 0; ei < edges.size(); ++ei) {
        Edge &e = edges[ei];
        if (correct && labels[e.a] != labels[e.b]) correct = false;
        float c = 1.f - e.w / cscale;
        if (correct) Color(1,c,c);
        else Color(0,1,0);
        LineCon(pts[e.a], pts[e.b], .01f);
    }
    Color(1,1,1);
}

void SpringIt(vector<Point> &pts, vector<Point> &vels, vector<Edge> &edges, float cscale) {
    const float k = .001f;
    const float drag = .995f;
    const float eqdist = .01f;
    for (Edge &e : edges) {
        float power = 1.f / (1.f - e.w / cscale * 3 / 4.f);
        Point disp = pts[e.a] - pts[e.b];
        float L = Length(disp) - eqdist * power;
        disp = Normalize(disp);
        vels[e.a] -= disp * L * k;
        vels[e.b] += disp * L * k;
    }
    const float G = .00000001f;
    const float epsilon = .001f;
    for (int i = 0; i < pts.size(); ++i) {
        for (int j = 0; j < pts.size(); ++j) {
            if (i == j) continue;
            Point disp = pts[i] - pts[j];
            float r2 = LengthSq(disp);
            if (r2 < epsilon) r2 = epsilon;
            vels[i] += disp / pow(r2, 1.5f) * G;
        }
    }
    for (int i = 0; i < pts.size(); ++i) {
        pts[i] += vels[i];
        DragOnScreen(pts[i], .05f);
        vels[i] *= drag;
    }
}

void DrawKruskal(int n, vector<Edge> &edges, vector<string> &labels) {
    vector<Edge> mst;
    float maxw = Kruskal(n, edges, mst);
    bool correct = true;
    cout << "Correct list:" << endl;
    for (Edge &e : mst) {
        if (correct && labels[e.a] != labels[e.b]) {
            cout << "Incorrect list:" << endl;
            correct = false;
        }
        cout << labels[e.a] << ' ' << labels[e.b] << ' ' << e.w << endl;
    }
    cout << "done" << endl;
    const float radius = .45f;
    const Point center(.5f,.5f);
    vector<Point> pts, vels;
    for (int i = 0; i < n; ++i) {
        Point p = center + DirPt(M_PI * .5f - 2.f * M_PI * i / n) * radius;
        pts.push_back(p);
        vels.push_back(Point(0,0));
    }
    while (true) {
        DrawClear();
        DrawGraph(pts, mst, labels, maxw);
        DrawSwap();
        usleep(1000.f);
        SpringIt(pts, vels, mst, maxw);
    }
    
    
}


#endif
