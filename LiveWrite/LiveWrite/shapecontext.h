//
//  shapecontext.h
//  LiveWrite
//
//  Created by Ben Mildenhall on 11/23/13.
//  Copyright (c) 2013 Ben Mildenhall. All rights reserved.
//

#ifndef LiveWrite_shapecontext_h
#define LiveWrite_shapecontext_h

#include "utilities.h"
#include "stroke.h"

using namespace std;

class ShapeContext {
public:
    ShapeContext() { }
    ShapeContext(const Glyph &g): glyph(g), npts(0) {
        if (g.points.empty()) return;
        glyph.Normalize();
        glyph = SmoothGlyph(glyph, 500, 50, .01f);
        npts = (int)glyph.points.size();
        rthetas.resize(npts);
        discrete.resize(npts);
        for (int i = 0; i < npts; ++i) {
            rthetas[i].resize(npts);
            for (int j = 0; j < npts; ++j) {
                rthetas[i][j] = RTheta(glyph.points[j] - glyph.points[i]);
            }
        }
        rnum = thetanum = -1;
    }
    void Discretize(int rn, int tn) {
        if (rnum == rn && thetanum == tn) return;
        rnum = rn;
        thetanum = tn;
        int k = rnum * thetanum;
        float npi = 1./npts;
        for (int i = 0; i < npts; ++i) {
            discrete[i] = vector<float>(k, 0);
            for (int j = 0; j < npts; ++j) {
                int logr = (logf(rthetas[i][j].first) - kLogLB) / (kLogUB - kLogLB) * rnum;
                int theta = rthetas[i][j].second / (2.f * M_PI) * thetanum;
                logr = min(rnum - 1, max(0, logr));
                theta = min(thetanum - 1, max(0, theta));
                discrete[i][logr * thetanum + theta] += npi;
            }
        }
    }
    void Discretize(int rn, int tn, vector<vector<float>> &v) {
        Discretize(rn, tn);
        v = discrete;
    }
    bool empty() { return npts == 0; }
    friend float compare(ShapeContext &s, ShapeContext &t, int rn, int tn, float dbias);

private:
    const float kLogLB = -3.5f;
    const float kLogUB = .5f * logf(2.f);
    int npts;
    Glyph glyph;
    vector<vector<pair<float,float>>> rthetas;
    vector<vector<float>> discrete;
    int rnum, thetanum;
};

float ChiSquaredTest(const vector<float>& H1, const vector<float>& H2) {
    if (H1.size() != H2.size()) return -1.f;
    int n = (int)H1.size();
    float result = 0;
    for (int i = 0; i < n; ++i)
        if (H1[i] || H2[i])
            result += .5f * (H1[i] - H2[i]) * (H1[i] - H2[i]) / (H1[i] + H2[i]);
    return result;
}

void GetSCHistograms(Glyph &g, vector<vector<float>> &hist) {
    ShapeContext sc(g);
    sc.Discretize(5, 12, hist);
}

float compare(ShapeContext &s, ShapeContext &t, int rn, int tn, float dbias = 1.f) {
    if (s.empty() || t.empty()) return 0;
    s.Discretize(rn, tn);
    t.Discretize(rn, tn);
    char dumb;
    for (int i = 0; i < 0; ++i) {
        DrawClear();
        DrawMatrix(s.discrete[i], tn, rn, vecmax(s.discrete[i]));
        DrawSwap();
        cin >> dumb;
        DrawClear();
        DrawMatrix(t.discrete[i], tn, rn, vecmax(t.discrete[i]));
        DrawSwap();
        cin >> dumb;
    }
    int n = (int)s.glyph.points.size();

    vector<float> costmat(n * n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            costmat[j * n + i] = ChiSquaredTest(s.discrete[i], t.discrete[j]);
        }
    }

    return LCP(costmat, n, dbias);
}


#endif
