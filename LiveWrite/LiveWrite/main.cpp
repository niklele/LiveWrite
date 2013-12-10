//
//  main.cpp
//  LiveWrite
//
//  Created by Ben Mildenhall on 11/8/13.
//  Copyright (c) 2013 Ben Mildenhall. All rights reserved.
//

#include <iostream>
#include "utilities.h"
#include "stroke.h"
#include <map>
#include <unistd.h>
#include "subcipher.h"
#include "shapecontext.h"
#include "graph.h"


static float mouseX, mouseY;
static Point mouse, mouseOld;
const static int NRES = 100;
const static float STDEV = .03f;
const static int UPSCALE = 500;

void DisplayCallback();


struct TrainEx {
    Glyph g;
    char label;
};
static char currlabel = 'a';
static map<char, vector<Glyph>> glyphs;
static vector<Glyph> glyphVec;
static Glyph currGlyph;
static Glyph compGlyph;
static vector<float> LogRegTheta;
const static bool gridOn = true;

static string wdir;
void Init() {
    SeedStartTime();
    wdir = string(__FILE__);
    wdir = wdir.substr(0, wdir.find("main.cpp"));
    cout << "Working directory for files is: " << endl << wdir << endl;
}

void UpdateMouse(int x, int y) {
    mouseOld = mouse;
    mouse[0] = mouseX = (float)x/width;
    mouse[1] = mouseY = 1. - (float)y/height;
}

Point MouseDelta(int x, int y) {
    float xn = (float)x/width;
    float yn = 1. - (float)y/height;
    return Point(xn - mouseX, yn - mouseY);
}

char Analyze(Glyph &g, bool verbose = false);

const static int gwidth = 8, gheight = 6;
static int gridentry = 0;
void CheckGrid();

static Timer timer;
static bool timerValid = false;
int lastbox = -1;
void MouseCallback(int button, int state, int x, int y) {
    //button = GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON
    //state = GLUT_UP or GLUT_DOWN
    UpdateMouse(x,y);
    if (state == GLUT_UP) {
        currGlyph.AddPoint(mouse, ElapsedMillis());
        Pause();
        currGlyph.MouseUp();
        currGlyph.penLift.push_back(true);
        timer.SetTimer(500.f);
        //timerValid = true;
    }
    if (state == GLUT_DOWN) {
        Unpause();
        if (gridOn) CheckGrid();
        currGlyph.AddPoint(mouse, ElapsedMillis());
        currGlyph.penLift.push_back(false);
        currGlyph.MouseDown();
        timerValid = false;
    }
    DisplayCallback();
}


void MouseMotionCallback(int x, int y) {
    UpdateMouse(x,y);
    currGlyph.AddPoint(mouse, ElapsedMillis());
    DisplayCallback();
}

void PassiveMotionCallback(int x, int y) {
    UpdateMouse(x,y);
}

// Adds glyph to glyph map and starts blank glyph
void NewGlyph() {
    if (currGlyph.points.size()) {
        cout << "Normalizing..." << endl;
        currGlyph.Normalize();
        DrawClear();
        currGlyph.Draw();
        DrawSwap();
        glyphs[currlabel].push_back(currGlyph);
        currGlyph.Empty();
        cout << "Recorded glyph with label " << currlabel;
        cout << " (example " << glyphs[currlabel].size() << ")" << endl;
    }
}

void WriteData(string filename, map<char, vector<Glyph>> &data) {
    string path = wdir + filename;
    ofstream outfile(path);
    int nex = 0;
    for (auto &p : data) nex += p.second.size();
    for (auto &p : data) {
        for (auto &g : p.second) {
            outfile << "Label: " << p.first << endl;
            size_t npoints = g.points.size();
            outfile << "Points: " << npoints << endl;
            outfile << "Width: " << g.width << "\nDuration: " << g.duration ;
            outfile << "\nCentroid: " << g.centroid[0] << " " << g.centroid[1] << endl;
            for (size_t i = 0; i < npoints; ++i)
                outfile << g.times[i] << " " << g.points[i][0] << " " << g.points[i][1] << endl;
            outfile << endl;
        }
    }
    outfile.close();
}

void WriteVectorData() {
    string path = wdir + "vecdata.dat";
    ofstream outfile(path);
    for (auto &g : glyphVec) {
        size_t npoints = g.points.size();
        outfile << "Points: " << npoints << endl;
        outfile << "Width: " << g.width << "\nDuration: " << g.duration ;
        outfile << "\nCentroid: " << g.centroid[0] << " " << g.centroid[1] << endl;
        for (size_t i = 0; i < npoints; ++i)
            outfile << g.times[i] << " " << g.points[i][0] << " " << g.points[i][1] << endl;
        outfile << endl;
    }
    outfile.close();
}


void PlotCurvature(Glyph &g, float time, bool integrate = true) {
    if (g.times.size() < 2) return;
    LineCon(0, .5, 1, .5);
    Color(1,0,0);
    Point prev = g.points[1] - g.points[0];
    Point graphpt(0,.5);
    Point graphint(0,.5);
    for (int i = 2; i < g.times.size(); ++i) {
        if (g.times[i] >= time) Color(0,1,0);
        else Color(1,0,0);
        Point next = g.points[i] - g.points[i-1];
        Point ngpt(g.times[i], .5f + Bearing(prev, next) / M_PI * .5f);
        Point ngipt(g.times[i], graphint[1] + (ngpt[1] - .5f) * .25f);
        Circle(ngpt, .001);
        
        if (g.times[i] >= time) Color(.3,.3,1);
        else Color(1,.5,0);
        if (integrate) Circle(ngipt, .001);
        
        graphpt = ngpt;
        graphint = ngipt;
        prev = next;
        
    }
    Color(1,1,1);
}

void PlotVector(vector<float> vals, float bound) {
    LineCon(0,.5f,1,.5f);
    for (int i = 0; i < vals.size(); ++i)
        Circle((float)i / vals.size(), .5f + .5f * vals[i] / bound, .003);
}

// Makes histogram of curvatures for Naive Bayes
void DumpCurvatures(Glyph &g, vector<int> &curvs, int nres, bool verbose = false) {
    if (g.times.size() < 2) return;
    Point prev = g.points[1] - g.points[0];
    if (verbose) cout << "\nBeginning..." << endl;
    for (int i = 2; i < g.times.size(); ++i) {
        Point next = g.points[i] - g.points[i-1];
        int xt = nres * g.times[i];
        int yc = nres * (.5f + Bearing(prev, next) / M_PI * .5f);
        if (verbose) cout << "Adding at timeval " << xt << " curval " << yc << endl;
        if (verbose) cout << "Bearing here is " << Bearing(prev, next) / M_PI << endl;
        curvs[xt * nres + yc]++;
        prev = next;
    }
}

// Plots a Naive Bayes curvature discretized histogram
void DrawCurvHist(vector<int> &curvs, int nres) {
    int max = 0;
    for (int i : curvs) max = std::max(i,max);
    for (int x = 0; x < nres * nres; ++x) {
        Color(0, 0, (float)(curvs[x]) / (max));
        DrawBox((float)(x / nres) / nres, (float)(x % nres) / nres, 1.f / nres, 1.f / nres);
    }
    Color(1,1,1);
}

void SimulateBigGlyph(Glyph &g, float totaltime, float extend,
                                    vector<int> &curvs, int nres) {
    float time;
    SeedStartTime();
    while ((time = ElapsedMillis()) <= totaltime * extend) {
        DrawClear();
        //DrawCurvHist(curvs, nres);
        PlotCurvature(g, time/totaltime);
        for (int i = 0; i < g.times.size(); ++i) {
            if (totaltime * g.times[i] < time) Circle(g.points[i], .001);
        }
        DrawSwap();
        //usleep(10000);
    }
}

void SimulateGlyph(Glyph &g, float totaltime = 1000.f, float extend = 1.5f) {
    cout << "running sim...";
    float time;
    SeedStartTime();
    while ((time = ElapsedMillis()) <= totaltime * extend) {
        DrawClear();
        PlotCurvature(g, time/totaltime);
        for (int i = 0; i < g.times.size(); ++i) {
            if (totaltime * g.times[i] < time) Circle(g.points[i], .001);
        }
        DrawSwap();
        usleep(10000);
    }
    cout << "over." << endl;
}

void PlotGlyph(Glyph &g, float t) {
    Color(1,1,1);
    for (int i = 0; i < g.times.size(); ++i)
        if (g.times[i] < t) Circle(g.points[i], .001);
}



float ShapeDistance(Glyph &a, Glyph &b, vector<float> &f) {
    if (a.points.empty() || b.points.empty()) return 0;
    ShapeContext s(a), t(b);
    f.push_back(1.f);
    f.push_back(compare(s, t, 5, 12));
    
    /*
    Glyph aa = SmoothGlyph(a, 500, 100, .01);
    Glyph bb = SmoothGlyph(b, 500, 100, .01);
    f.push_back(LeastCost<Point>(aa.points, bb.points,
                                 [] (const Point &a, const Point &b)
                                 { return LengthSq(a - b); }, 1.0f, false));
    
    aa = SmoothSpaceGlyph(a, 500, 100, .01);
    bb = SmoothSpaceGlyph(b, 500, 100, .01);
    f.push_back(LeastCost<Point>(aa.points, bb.points,
                                 [] (const Point &a, const Point &b)
                                 { return LengthSq(a - b); }, 1.0f, false));
    */
    //return f[1] + .44f * f[2] + .29f * f[3];
    return f[1];
}

float ShapeDistance(Glyph &a, Glyph &b) {
    vector<float> dummy;
    return ShapeDistance(a, b, dummy);
}

class Cluster {
public:
    float GetDistance(Glyph &g, bool verbose = false) {
        if (verbose) cout << GetLabels() << ": ";
        float score = 0;
        vector<float> curv;
        
        GetCurvature(g, curv);
        float curr = LeastCostFloat(curv, cf, 1.5f, false) * 20.f;
        if (verbose) cout << curr/20.f << " ";
        score += curr;
        curv.clear();
        
        GetCurvReverse(g, curv);
        curr = LeastCostFloat(curv, cr, 1.5f, false) * 20.f;
        if (verbose) cout << curr/20.f << " ";
        score += curr;
        curv.clear();
        
        GetCurvature(g, curv, true);
        curr = LeastCostFloat(curv, cif, 1.3f, false) * 4.f;
        if (verbose) cout << curr/4.f << " ";
        score += curr;
        curv.clear();
        
        GetCurvReverse(g, curv, true);
        curr = LeastCostFloat(curv, cir, 1.3f, false) * 4.f;
        if (verbose) cout << curr /4.f << " ";
        score += curr;
        
        curr = LeastCost<Point>(g.points, average.points,
                                [](const Point &a, const Point &b)
                                { return LengthSq(a - b); }, 1.2f, false) * 15.f;
        if (verbose) cout << curr /15.f << " ";
        score += curr;
        if (verbose) cout << "= " << score << endl;
        return score;
    }
    float GetDistanceThreaded(Glyph &g, bool verbose = false) {
        if (verbose) cout << GetLabels() << ": ";
        int numfeat = 5;
        float feats[5];
        float weights[] = {100.f, 100.f, 1.f, 1.f, 40.f};
        //float weights[] = {0,0,0,0,50.};
        thread threads[5];
        
        threads[0] = thread([this, &g] (float *dump) {
            vector<float> curv;
            GetCurvature(g, curv);
            *dump = LeastCostFloat(curv, cf, 1.0f, false);
        }, feats + 0);
        
        threads[1] = thread([this, &g] (float *dump) {
            vector<float> curv;
            GetCurvReverse(g, curv);
            *dump = LeastCostFloat(curv, cr, 1.0f, false);
        }, feats + 1);
        
        threads[2] = thread([this, &g] (float *dump) {
            vector<float> curv;
            GetCurvature(g, curv, true);
            *dump = LeastCostFloat(curv, cif, 1.0f, false);
        }, feats + 2);
        
        threads[3] = thread([this, &g] (float *dump) {
            vector<float> curv;
            GetCurvReverse(g, curv, true);
            *dump = LeastCostFloat(curv, cir, 1.0f, false);
        }, feats + 3);
        
        threads[4] = thread([this, &g] (float *dump) {
            *dump = LeastCost<Point>(g.points, average.points,
                                [] (const Point &a, const Point &b)
                                { return LengthSq(a - b); }, 1.0f, false);
        }, feats + 4);
        
        for (int i = 0; i < numfeat; ++i)
            threads[i].join();
        
        float score = 0;
        for (int i = 0; i < numfeat; ++i) {
            float curr = weights[i] * feats[i];
            score += curr;
            if (verbose) cout << curr << " ";
        }
        if (verbose) cout << "= " << score << endl;
        return score;
    }
    float GetShapeContextDiff(Glyph &g) {
        //return ShapeDistance(average, g);
        ShapeContext me(average);
        ShapeContext them(g);
        return compare(me, them, 5, 12);
    }
    void Recalculate() {
        average = ReduceGlyph(composite, NRES, STDEV);
        cf.clear(); cr.clear(); cif.clear(); cir.clear();
        GetCurvature(average, cf, false);
        GetCurvReverse(average, cr, false);
        GetCurvature(average, cif, true);
        GetCurvReverse(average, cir, true);
    }
    void Insert(Glyph &g, char label, bool recalc = true) {
        glyphs.push_back(g);
        labels += label;
        composite.points.insert(composite.points.end(), g.points.begin(), g.points.end());
        composite.times.insert(composite.times.end(), g.times.begin(), g.times.end());
        if (recalc) Recalculate();
    }
    void Clear() {
        average.Empty();
        cf.clear(); cr.clear(); cif.clear(); cir.clear();
        composite.Empty();
        glyphs.clear();
        labels.clear();

    }
    string GetLabels() {
        return labels;
    }
    float AverageScore() {
        float ret = 0;
        for (Glyph &g : glyphs) ret += GetDistance(g);
        return ret / glyphs.size();
    }
    float AverageScoreThreaded() {
        float ret = 0;
        mutex retLock;
        vector<thread> threads;
        for (Glyph &g : glyphs)
            threads.push_back(thread([&ret, &retLock, &g, this] () {
                float val = GetDistanceThreaded(g);
                lock_guard<mutex> lg(retLock);
                ret += val;
            }));
        for (thread &t : threads) t.join();
        return ret / glyphs.size();
    }
    char GetLabel() {
        map<char, int> counts;
        for (char c : labels) counts[c]++;
        map<int, char> inverse;
        for (auto &p : counts) inverse[p.second] = p.first;
        return inverse.begin()->second;
    }
    Glyph average;
    vector<float> cf, cr, cif, cir;
    Glyph composite;
    vector<Glyph> glyphs;
    string labels;
};

static vector<Cluster> clusters;

void ReadData(string filename, map<char, vector<Glyph>> &data) {
    string path = wdir + filename;
    ifstream infile(path);
    string dummy;
    int nexamples = 0;
    while (true) {
        TrainEx curr;
        int npoints;
        infile >> dummy >> curr.label >> dummy >> npoints;
        if (infile.eof()) break;
        cout << "Read in label " << curr.label << " with " << npoints << " points" << endl;
        infile >> dummy >> curr.g.width >> dummy >> curr.g.duration;
        infile >> dummy >> curr.g.centroid[0] >> curr.g.centroid[1];
        for (int j = 0; j < npoints; ++j) {
            float t;
            Point p;
            infile >> t >> p[0] >> p[1];
            curr.g.times.push_back(t);
            curr.g.points.push_back(p);
        }
        //curr.g.Normalize();
        data[curr.label].push_back(curr.g);
        ++nexamples;
    }
    infile.close();
    cout << "Read in " << nexamples << " total examples. " << endl;
}

void ReadDataSimple(string filename, map<char, vector<Glyph>> &data) {
    string path = wdir + filename;
    ifstream infile(path);
    string dummy;
    int nexamples = 0;
    while (true) {
        TrainEx curr;
        int npoints;
        infile >> dummy >> curr.label >> dummy >> npoints;
        if (infile.eof()) break;
        cout << "Read in label " << curr.label << " with " << npoints << " points" << endl;
        for (int j = 0; j < npoints; ++j) {
            float t;
            Point p;
            infile >> t >> p[0] >> p[1];
            curr.g.times.push_back(t);
            curr.g.points.push_back(p);
        }
        //curr.g.Normalize();
        data[curr.label].push_back(curr.g);
        ++nexamples;
    }
    infile.close();
    cout << "Read in " << nexamples << " total examples. " << endl;
}

void ReadVectorData(string filename = "vecdata.dat") {
    if (!glyphVec.empty()) glyphVec.push_back(Glyph());
    string path = wdir + filename;
    ifstream infile(path);
    string dummy;
    int nexamples = 0;
    while (true) {
        TrainEx curr;
        int npoints;
        infile >> dummy >> npoints;
        if (infile.eof()) break;
        cout << "Read in glyph with " << npoints << " points" << endl;
        infile >> dummy >> curr.g.width >> dummy >> curr.g.duration;
        infile >> dummy >> curr.g.centroid[0] >> curr.g.centroid[1];
        for (int j = 0; j < npoints; ++j) {
            float t;
            Point p;
            infile >> t >> p[0] >> p[1];
            curr.g.times.push_back(t);
            curr.g.points.push_back(p);
        }
        glyphVec.push_back(curr.g);
        ++nexamples;
    }
    infile.close();
    cout << "Read in " << nexamples << " total examples. " << endl;
}

void CreateClusters(bool threaded = false) {
    clusters.clear();
    for (auto &p : glyphs) {
        char label = p.first;
        vector<Glyph> &vec = p.second;
        for (Glyph &g : vec) {
            Glyph smoothed = SmoothGlyph(g, UPSCALE, NRES, STDEV);
            Cluster *best = nullptr;
            float score = INFINITY;
            for (Cluster &c : clusters) {
                if (c.GetLabels()[0] != label) continue;
                float sc = threaded ? c.GetDistanceThreaded(smoothed) : c.GetDistance(smoothed);
                if (sc < score) {
                    best = &c;
                    score = sc;
                }
            }
            if (score > 30) {
                if (best)
                    cout << "Best cluster was " << best->GetLabels() << ", score " << score << endl;
                clusters.push_back(Cluster());
                best = &clusters.back();
                cout << "New cluster. ";
            }
            cout << "Adding a \'" << label << "\' to cluster " << best->GetLabels() << endl;
            best->Insert(smoothed, label);
        }
    }
    /*
    auto iter = clusters.begin();
    for (int i = 0; i < clusters.size(); ++i) {
        if (clusters[i].glyphs.size() == 1) {
            clusters.erase(iter);
            i--;
        } else iter++;
    }
     */
    for (int i = 0; i < clusters.size(); ++i) {
        cout << "#" << i << " (" << (threaded ? clusters[i].AverageScoreThreaded() : clusters[i].AverageScore())
            << "): " << clusters[i].GetLabels() << endl;
    }
}

void CreateClusters2() {
    clusters.clear();
    for (auto &p : glyphs) {
        char label = p.first;
        vector<Glyph> &vec = p.second;
        for (Glyph &g : vec) {
            Glyph smoothed = SmoothGlyph(g, UPSCALE, NRES, STDEV);
            Cluster *best = nullptr;
            float score = INFINITY;
            for (Cluster &c : clusters) {
                if (c.GetLabels()[0] != label) continue;
                float sc = c.GetDistanceThreaded(smoothed);
                if (sc < score) {
                    best = &c;
                    score = sc;
                }
            }
            if (score > 30) {
                if (best)
                    cout << "Best cluster was " << best->GetLabels() << ", score " << score << endl;
                clusters.push_back(Cluster());
                best = &clusters.back();
                cout << "New cluster. ";
            }
            cout << "Adding a \'" << label << "\' to cluster " << best->GetLabels() << endl;
            best->Insert(g, label);
        }
    }
    /*
     auto iter = clusters.begin();
     for (int i = 0; i < clusters.size(); ++i) {
     if (clusters[i].glyphs.size() == 1) {
     clusters.erase(iter);
     i--;
     } else iter++;
     }
     */
    for (int i = 0; i < clusters.size(); ++i) {
        cout << "#" << i << " (" << clusters[i].AverageScoreThreaded()
            << "): " << clusters[i].GetLabels() << endl;
    }
}

void CreateClusters3() {
    clusters.clear();
    for (auto &p : glyphs) {
        char label = p.first;
        vector<Glyph> &vec = p.second;
        for (Glyph &g : vec) {
            Glyph smoothed = SmoothGlyph(g, UPSCALE, NRES, STDEV);
            Cluster *best = nullptr;
            float score = INFINITY;
            for (Cluster &c : clusters) {
                if (c.GetLabels()[0] != label) continue;
                float sc = c.GetShapeContextDiff(smoothed);
                if (sc < score) {
                    best = &c;
                    score = sc;
                }
            }
            if (score > 15) {
                if (best)
                    cout << "Best cluster was " << best->GetLabels() << ", score " << score << endl;
                clusters.push_back(Cluster());
                best = &clusters.back();
                cout << "New cluster. ";
            }
            cout << "Adding a \'" << label << "\' to cluster " << best->GetLabels();
            cout << " (score " << score << ")" << endl;
            best->Insert(g, label);
        }
    }
    auto iter = clusters.begin();
    for (int i = 0; i < clusters.size(); ++i) {
        if (clusters[i].glyphs.size() == 1) {
            clusters.erase(iter);
            i--;
        } else iter++;
    }
    for (int i = 0; i < clusters.size(); ++i) {
        cout << "#" << i << " (" << clusters[i].AverageScoreThreaded()
        << "): " << clusters[i].GetLabels() << endl;
    }
}


static Point dumbpoint;
char Analyze(Glyph &g, bool verbose) {
    if (g.points.size()) {
        g.Normalize();
        dumbpoint = g.points[10]; // Debugging crap
        Glyph currInt = SmoothGlyph(g, UPSCALE, NRES, STDEV);
        dumbpoint = currInt.points[10]; // Debugging crap
        cout << dumbpoint << endl;
        g.Empty();
        if (clusters.empty()) return 0;
        map<float, string> scores;
        for (Cluster &c : clusters) {
            //float curr = c.GetDistanceThreaded(currInt, false);
            float curr = c.GetShapeContextDiff(currInt);
            scores[curr] = c.GetLabels();
        }
        if (verbose) {
            auto iter = scores.begin();
            for (int i = 0; i < 5; ++i) {
                if (iter == scores.end()) break;
                cout << iter->second << ": " << iter->first << endl;
                iter++;
            }
            cout << endl;
        }
        cout << " score is " << scores.begin()->first << ", char " << scores.begin()->second << endl;
        return scores.begin()->second[0];
        
    } else return 0;
}

void CalculateError() {
    int classified(0);
    int wrong(0);
    for (auto &p : glyphs) {
        for (auto &g : p.second) {
            Glyph copy = g;
            char guess = Analyze(copy, false);
            if (guess != p.first) {
                copy = g;
                Glyph currInt = SmoothGlyph(copy, UPSCALE, NRES, STDEV);
                ++wrong;
                cout << "Misclassified " << p.first << " as " << guess << endl;
                for (Cluster &c : clusters) {
                    if (c.GetLabels()[0] == guess || c.GetLabels()[0] == p.first) {
                        c.GetDistanceThreaded(currInt, true);
                    }
                }
            }
            ++classified;
        }
    }
    cout << "Accuracy was " << 1. - (float)wrong / classified << endl;
    cout << " (ie " << wrong << " wrong out of " << classified << ")" << endl;
}


static char cOne = 's', cTwo = 's';
static int coneIndex = 11, ctwoIndex = 10;
void ShowMe() {
    /*
    DrawClear();
    vector<float> curv1, curv2;
    Glyph g1 = SmoothGlyph(glyphs[cOne][coneIndex], 800, 200, .01);
    Glyph g2 = SmoothGlyph(glyphs[cTwo][ctwoIndex], 800, 200, .01);
    GetCurvature(g1, curv1, true);
    GetCurvature(g2, curv2, true);
    LeastCostFloat(curv1, curv2, 1.f);
    PlotCurvature(g1, 1, true);
    PlotCurvature(g2, 0, true);
    cout << "Displaying for " << cOne << " (" << coneIndex << ")"
    << " vs " << cTwo << " (" << ctwoIndex << ")" << endl;
    DrawSwap();
     */
    DrawClear();
    const float velbound = 20.f;
    vector<float> vel1, vel2;
    Glyph g1 = SmoothGlyph(glyphs[cOne][coneIndex], 800, 200, .01);
    Glyph g2 = SmoothGlyph(glyphs[cTwo][ctwoIndex], 800, 200, .01);
    GetStrokeSpeeds(g1, vel1);
    GetStrokeSpeeds(g2, vel2);
    cout << vel1 << endl;
    cout << vel2 << endl;
    cout << *max_element(vel1.begin(), vel1.end()) << " for 1, ";
    cout << *max_element(vel2.begin(), vel2.end()) << " for 2" << endl;
    float lcp = LeastCostFloat(vel1, vel2, 1.f);
    float sd = SquareDif(vel1, vel2);
    cout << "Square dif is " << sd << ", ratio " << sd / lcp << endl;
    Color(1,0,0);
    PlotVector(vel1, velbound);
    Color(0,1,0);
    PlotVector(vel2, velbound);
    cout << "Displaying for " << cOne << " (" << coneIndex << ")"
    << " vs " << cTwo << " (" << ctwoIndex << ")" << endl;
    DrawSwap();
}

void ShowMeSpace() {
    DrawClear();
    LeastCost<Point>(glyphs[cOne][coneIndex].points, glyphs[cTwo][ctwoIndex].points,
              [](const Point &a, const Point &b) { return LengthSq(a - b); }, 1.3f);
    glyphs[cOne][coneIndex].Draw();
    glyphs[cTwo][ctwoIndex].Draw();
    cout << "Displaying for " << cOne << " (" << coneIndex << ")"
    << " vs " << cTwo << " (" << ctwoIndex << ")" << endl;
    DrawSwap();
}

void DrawCompare(vector<float> &one, vector<float> &two, vector<pair<int,int>> &ind) {
    float time, totaltime = 2000.f, extend = 2.f;
    SeedStartTime();
    while ((time = ElapsedMillis()) <= totaltime * extend) {
        int frac = two.size() * time / totaltime;
        DrawClear();
        Color(1,1,1);
        LineCon(0, .5, 1, .5);
        for (int i = 0; i < ind.size(); ++i) {
            if (ind[i].second > frac) break;
            float efftime = (float)ind[i].second / two.size();
            Color(1,0,0);
            Circle(efftime, .5f + two[ind[i].second], .003);
            Color(0,1,0);
            Circle(efftime, .5f + one[ind[i].first], .003);
        }
        for (int i = 0; i < one.size(); ++i) {
            if ((float)i / one.size() > time / totaltime) break;
            Color(0,0,1);
            Circle((float)i / one.size(), .5f + one[i], .003);
        }
        DrawSwap();
        usleep(10000);
    }
}

void HighlightMatches(string str) {
    Color(0,1,0);
    int n = (int)str.length();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && str[i] == str[j])
                Circle((i + .5f) / n, (j + .5f) / n, .015);
        }
    }
}

float GetGlyphDistance(Glyph &g, Glyph &h, vector<float> &feat, bool verbose = false) {
    int numfeat = 5;
    feat.clear();
    feat.push_back(1.f);
    float feats[5];
    //float weights[] = {1000.f, 1000.f, 1.f, 1.f, 40.f};
    //float weights[] = {0,0,0,0,50.};
    //float weights[] = {-78201.6 , 1030.86 , -10563,  -1284.93,  3180.93 , 17343.8};
    //float weights[] = {65.1378 , 0.220724 , 0.18566 , -32.5554  ,-115.064 , -6.89492};
    float weights [] = { 1 , -1030.86 , 10563 , 1284.93 , -3180.93 , -17343.8 };
    thread threads[5];
    
    threads[0] = thread([&g, &h] (float *dump) {
        vector<float> curvg, curvh;
        GetCurvature(g, curvg);
        GetCurvature(h, curvh);
        *dump = LeastCostFloat(curvg, curvh, 1.0f, false);
    }, feats + 0);
    
    threads[1] = thread([&g, &h] (float *dump) {
        vector<float> curvg, curvh;
        GetCurvReverse(g, curvg);
        GetCurvReverse(h, curvh);
        *dump = LeastCostFloat(curvg, curvh, 1.0f, false);
    }, feats + 1);
    
    threads[2] = thread([&g, &h] (float *dump) {
        vector<float> curvg, curvh;
        GetCurvature(g, curvg, true);
        GetCurvature(h, curvh, true);
        *dump = LeastCostFloat(curvg, curvh, 1.0f, false);
    }, feats + 2);
    
    threads[3] = thread([&g, &h] (float *dump) {
        vector<float> curvg, curvh;
        GetCurvReverse(g, curvg, true);
        GetCurvReverse(h, curvh, true);
        *dump = LeastCostFloat(curvg, curvh, 1.0f, false);
    }, feats + 3);
    
    threads[4] = thread([&g, &h] (float *dump) {
        *dump = LeastCost<Point>(g.points, h.points,
                                 [] (const Point &a, const Point &b)
                                 { return LengthSq(a - b); }, 1.0f, false);
    }, feats + 4);
    
    for (int i = 0; i < numfeat; ++i)
        threads[i].join();
    
    float score = weights[0];
    for (int i = 0; i < numfeat; ++i) {
        float curr = weights[i+1] * feats[i];
        score += curr;
        feat.push_back(feats[i]);
        if (verbose) cout << curr << " ";
    }
    if (verbose) cout << "= " << score << endl;
    return score;
}

inline float Sigmoid(float z) { return 1.f / (1.f + expf(-z)); }

vector<float> LogisticRegression(vector<vector<float>> &X, vector<float> &y, int iters = 10000) {
    float alpha = .001f;
    int m = (int)X.size();
    cout << X[100][0] << ' ';
    int n = (int)X[0].size();
    vector<float> theta(n, 0);
    cout << "Beginning logistic regression with n = " << n << ", m = " << m << endl;
    cout << "and " << iters << " iterations." << endl;
    for (int iter = 0; iter < iters; ++iter) {
        vector<float> newtheta = theta;
        for (int i = 0; i < m; ++i) {
            float hyp = 0;
            for (int j = 0; j < n; ++j)
                hyp += theta[j] * X[i][j];
            hyp = Sigmoid(hyp);
            
            if (y[i]) alpha *= 100.f;
            for (int j = 0; j < n; ++j)
                newtheta[j] += alpha * (y[i] - hyp) * X[i][j];
            if (y[i]) alpha /= 100.f;
            
        }
        theta = newtheta;
        if (!(iter % 100)) cout << ".";
    }
    cout << endl << theta << endl;
    for (int j = 1; j < n; ++j)
        theta[j] /= theta[0];
    theta[0] = 1.f;
    cout << "Logistic regression produced the following theta parameters" << endl;
    cout << "in " << iters << " iterations and using alpha value " << alpha << ": " << endl;
    cout << theta << endl;
    int fails = 0;
    int fails1 = 0, fails0 = 0, tot1 = 0, tot0 = 0;
    for (int i = 0; i < m; ++i) {
        float hyp = 0.f;
        for (int j = 0; j < n; ++j)
            hyp += theta[j] * X[i][j];
        float asgn = (Sigmoid(hyp) >= .5f ? 1.f : 0.f);
        if (y[i] == 0.f) tot0++;
        else tot1++;
        if (asgn != y[i]) {
            if (y[i] == 0.f) fails0++;
            else fails1++;
            fails++;
            //cout << "Failure on data point " << i << ", supposed to be " << y[i]
            //    << ", but hyp was " << hyp << ", assigned to " << asgn << endl;
        }
    }
    cout << "Fails on " << fails << " examples" << endl;
    cout << "0 fails : " << fails0 << " of " << tot0 << " (" << (float)fails0/tot0 << ")" << endl;
    cout << "1 fails : " << fails1 << " of " << tot1 << " (" << (float)fails1/tot1 << ")" << endl;
    return theta;
}

    
float ClusterMinDist(Cluster &c, Cluster &d) {
    float min = INFINITY;
    for (Glyph &g : c.glyphs) {
        for (Glyph &h : d.glyphs) {
            float curr = ShapeDistance(g, h);
            if (curr < min) min = curr;
        }
    }
    return min;
}

float ClusterVecMinDist(vector<Cluster> &cs, int i, int &j) {
    float glow = INFINITY;
    for (int k = 0; k < cs.size(); ++k) {
        if (k == i) continue;
        float curr = ClusterMinDist(cs[i], cs[k]);
        if (curr < glow) {
            glow = curr;
            j = k;
        }
    }
    return glow;
}

mutex cmLock;
static const int NTHREADS = 6;
static atomic<int> tofinish(NTHREADS);
template <class T>
void FillMatRange(vector<float> &cmat, vector<vector<T>> &timeseries,
                  function<float(const T&, const T&)> compareT,
                  int n, int start, int end) {
    
    for (int j = start; j < end; ++j) {
        for (int i = j + 1; i < n; ++i) {
            float val = LeastCost<T>(timeseries[i],
                                    timeseries[j], compareT, 1.f, false);
            lock_guard<mutex> lg(cmLock);
            cmat[j * n + i] = val;
        }
    }
    --tofinish;
}

void AnalyzeCostMat(vector<float> &costmat, string &labels) {
    int n = (int)labels.size();
    float smin = INFINITY, smax = -INFINITY;
    int dmin = INFINITY, dmax = -INFINITY;
    for (int ind = 0; ind < n * n; ++ind) {
        int i = ind % n, j = ind / n;
        if (i <= j) continue;
        float c = costmat[ind];
        if (labels[i] == labels[j]) {
            if (smin > c) smin = c;
            if (smax < c) smax = c;
        } else {
            if (dmin > c) dmin = c;
            if (dmax < c) dmax = c;
        }
    }
    cout << "Same range is " << smin << " to " << smax << endl;
    cout << "Diff range is " << dmin << " to " << dmax << endl;
    int sind = 0, dins = 0;
    int sames = 0, diffs = 0;
    for (int ind = 0; ind < n * n; ++ind) {
        int i = ind % n, j = ind / n;
        if (i <= j) continue;
        float c = costmat[ind];
        bool same = labels[i] == labels[j];
        if (same) ++sames;
        else ++diffs;
        if (c <= smax && !same) ++dins;
        if (c >= dmin && same) ++sind;
    }
    cout << "Found " << dins << " diffs in the \"same\" range, of " << diffs << " total diffs" << endl;
    cout << "Found " << sind << " sames in the \"diff\" range, of " << sames << " total sames" << endl;
    int miscat = 0;
    for (int j = 0; j < n; ++j) {
        float best = INFINITY;
        char match;
        for (int i = 0; i < n; ++i) {
            if (j == i) continue;
            if (costmat[j*n+i] < best) {
                best = costmat[j*n+i];
                match = labels[i];
            }
        }
        if (match != labels[j]) { ++miscat; cout << labels[j] << match << ' '; }
    }
    cout << endl;
    cout << "Miscat: " << miscat << ". As a fraction, " << (float)miscat / n << endl;
}

template <class T>
void EvaluateFeature(function<void(Glyph&, vector<T>&)> getSeries,
                     function<float(const T&, const T&)> compareT,
                     vector<vector<float>> &X) {
    vector<vector<T>> timeseries;
    string labels = "";
    cout << "Generating feature vectors" << endl;
    for (auto &p : glyphs) {
        for (Glyph &g : p.second) {
            labels += p.first;
            timeseries.push_back(vector<T>());
            getSeries(g, timeseries.back());
        }
    }
    
    cout << "Generating cost matrices" << endl;
    int n = (int)timeseries.size();
    int L = (int)timeseries.front().size();
    vector<float> costmatLCP(n * n);
    vector<float> costmatDirect(n * n);
    vector<thread> threads;
    tofinish = NTHREADS;
    for (int i = 0; i < NTHREADS; ++i) {
        int start = n * (1.f - sqrtf((float)(NTHREADS - i) / NTHREADS));
        int end = n * (1.f - sqrtf((float)(NTHREADS - i - 1) / NTHREADS));
        threads.push_back(thread(FillMatRange<T>, ref(costmatLCP), ref(timeseries),
                          compareT, n, start, end));
    }
    
    while (tofinish) {
        usleep(200 * 1000.f);
        lock_guard<mutex> lg(cmLock);
        DrawClear();
        DrawMatrix(costmatLCP, n, *max_element(costmatLCP.begin(), costmatLCP.end()));
        DrawSwap();
    }
    
    for (thread &t : threads)
        t.join();
    
    int fInd = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (i < j) {
                costmatLCP[j*n+i] = costmatLCP[i*n+j];
                costmatDirect[j*n+i] = costmatDirect[i*n+j];
            } else if (i == j) {
                costmatLCP[j*n+i] = 0;
                costmatDirect[j*n+i] = 0;
            } else {
                X[fInd++].push_back(costmatLCP[j*n+i]);
                //costmatLCP[j*n+i] = LeastCost<T>(timeseries[i], timeseries[j], compareT, 1.f, false);
                float dir = 0;
                for (int k = 0; k < L; ++k)
                    dir += compareT(timeseries[i][k], timeseries[j][k]);
                costmatDirect[j*n+i] = dir;
            }
        }
        DrawClear();
        DrawMatrix(costmatLCP, n, *max_element(costmatLCP.begin(), costmatLCP.end()));
        DrawSwap();
    }
    
    DrawClear();
    DrawMatrix(costmatDirect, n, *max_element(costmatDirect.begin(), costmatDirect.end()));
    DrawSwap();
    cout << "Analyzing direct comp" << endl;
    AnalyzeCostMat(costmatDirect, labels);
    cout << "Analyzing LCP comp" << endl;
    AnalyzeCostMat(costmatLCP, labels);
    cout << "Done evaluating!" << endl << endl;
    
}

string GetCipher(int n, vector<Edge>& mst, int numclusters) {
    cout << "Making cipher" << endl;
    int cutoff = (int)mst.size() - (numclusters - 1);
    char curr = 'a';
    string cipher(n, '.');
    for (int i = 0; i < glyphVec.size(); ++i)
        if (glyphVec[i].points.empty())
            cipher[i] = ' ';
    for (int i = 0; i < cutoff; ++i) {
        Edge &e = mst[i];
        if (cipher[e.a] != '.')
            cipher[e.b] = cipher[e.a];
        else if (cipher[e.b] != '.')
            cipher[e.a] = cipher[e.b];
        else cipher[e.a] = cipher[e.b] = curr++;
        for (int j = i + 1; j < cutoff; ++j) {
            Edge &d = mst[j];
            if (cipher[d.a] != '.')
                cipher[d.b] = cipher[d.a];
            else if (cipher[e.b] != '.')
                cipher[d.a] = cipher[d.b];
        }
        
    }
    for (int i = 0; i < cipher.size(); ++i)
        if (cipher[i] == '.')
            cipher[i] = curr++;
    return cipher;
}

void SolveCipher() {
    /*
    vector<Cluster> localClusters;
    char currlabel = 'a';
    vector<string> cipher;
    cipher.push_back(string());
    int n = (int)glyphVec.size();
    
    //string actual = "something else in this song cannot be remembered";
    //string actual = "the poison from a black widow spider is about fifteen times more potent";
    //string actual = "he found himself transformed in his bed into a horrible vermin";
    //actual = "the poison from a black widow spider is about fifteen times more potent he found himself transformed in his bed into a horrible vermin";
    //actual = "toaster and leader of a group of appliances consisting of a radio lamp electric blanket and vacuum cleaner";
    
    vector<string> nospace;
    for (char &c : actual)
        if (c != ' ') {
            string s = "";
            s += c;
            nospace.push_back(s);
        }
    vector<Edge> edges;
    
    vector<ShapeContext> contextVec;
    for (Glyph &g : glyphVec)
        contextVec.push_back(ShapeContext(g));
    
    vector<float> lcp(n * n, 0);
    float highest = 0.f;
    vector<vector<float>> features;
    for (int j = 0, jj = 0; j < n; ++j) {
        for (int i = 0, ii = 0; i < n; ++i) {
            vector<float> feat;
            lcp[j * n + i] = compare(contextVec[i], contextVec[j], 5, 12);
            //ShapeDistance(glyphVec[i], glyphVec[j], feat);
            //-GetGlyphDistance(smoothed[i], smoothed[j], feat);
            if (!feat.empty())
                features.push_back(feat);
            if (lcp[j * n + i] > highest) {
                highest = lcp[j * n + i];
                //cout << i << " = " << actual[i] << ", " << j
                  //  << " = " << actual[j] << ": " << lcp[j * n + i] << endl;
            }
            if ((lcp[j*n+i] && lcp[j*n+i] < 20) || actual[i] == actual[j])
                cout << actual[i] << actual[j] << ' ' << lcp[j*n+i] << "  ";
            
            if (lcp[j*n+i] > 0) {
                edges.push_back(Edge(ii, jj, lcp[j*n+i]));
            }
            if (actual[i] != ' ') ++ii;
        }
        cout << endl;
        DrawClear();
        DrawMatrix(lcp, n, highest);
        HighlightMatches(actual);
        DrawSwap();
        if (actual[j] != ' ') ++jj;
    }
    
    cout << "Edges has " << edges.size() << " edges. Tenth is " << edges[10].a;
    cout << " " << edges[10].b << " " << edges[10].w << endl;
    DrawClear();
    cout << "Starting..." << endl;
    DrawKruskal((int)nospace.size(), edges, nospace);
    cout << "Ended." << endl;
    DrawSwap();
    
    char d12;
    cin >> d12;
    
    DrawClear();
    DrawMatrix(lcp, n, highest);
    Color(0,1,0);
    int nstr = (int)actual.length();
    if (nstr != n) {
        cout << "Length mismatch, abort" << endl;
        return;
    }
    map<float, pair<char,char>> scorevalues;
    int needed = 0;
    vector<float> classes;
    for (int j = 0; j < nstr; ++j) {
        for (int i = 0; i < nstr; ++i) {
            classes.push_back(actual[i] == actual[j] ? 1.f : 0.f);
            if (actual[i] == actual[j]) {
                Circle((i + .5f) / nstr, (j + .5f) / nstr, .015);
                scorevalues[lcp[j * n + i]] = {actual[i], actual[j]};
                ++needed;
            }
        }
    }
    for (auto &p : scorevalues)
        cout << "(" << p.second.first << ", " << p.second.second << "): " << p.first << endl;
    
    int incorrect = 0, missing = 0;
    auto iter = scorevalues.end();
    iter--;
    float highestgood = iter->first;
    //highestgood = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (lcp[j * n + i] <= highestgood && actual[i] != actual[j]
                        && actual[i] != ' ' && actual[j] != ' ') {
                Color(1,0,0);
                Circle((i + .5f) / nstr, (j + .5f) / nstr, .015);
                ++incorrect;
                cout << "Wrong: " <<  i << " = " << actual[i] << ", " << j
                << " = " << actual[j] << ": " << lcp[j * n + i] << endl;
            }
            if (lcp[j * n + i] > highestgood && actual[i] == actual[j] && i != j) {
                ++missing;
                Color(0,1,1);
                Circle((i + .5f) / nstr, (j + .5f) / nstr, .015);
                cout << "Missing: " << i << " = " << actual[i] << ", " << j
                << " = " << actual[j] << ": " << lcp[j * n + i] << endl;
            }
        }
    }
    
    
    cout << "Highest good is " << highestgood << endl;
    cout << "Needed is " << needed / 2<< endl;
    cout << "Incorrect is " << incorrect / 2<< endl;
    cout << "Missing is " << missing / 2 << endl;
    cout << "Highest overall is " << highest << endl;
    DrawSwap();
    
    for (int i = 0; i < n; ++i) {
        if (actual[i] == ' ') continue;
        Glyph &g = glyphVec[i];
        bool inserted = false;
        for (Cluster &c : localClusters) {
            if (c.GetLabels()[0] == actual[i]) {
                inserted = true;
                c.Insert(g, actual[i]);
            }
        }
        if (!inserted) {
            localClusters.push_back(Cluster());
            localClusters.back().Insert(g, actual[i]);
        }
    }
    for (Cluster &c : localClusters) {
        float high = 0;
        for (int i = 0; i < c.glyphs.size(); ++i) {
            float glow = INFINITY;
            for (int j = 0; j < c.glyphs.size(); ++j) {
                if (i == j) continue;
                Glyph &g = c.glyphs[i], &h = c.glyphs[j];
                float curr = ShapeDistance(g, h);
                if (curr < glow) glow = curr;
            }
            if (glow > high) high = glow;
        }
        cout << c.GetLabels() << " has diameter ";
        cout << high << endl;
    }
    
    
    //LogisticRegression(features, classes);
    char dummy;
    cin >> dummy;
    localClusters.clear();
    
*/
    /*
    for (int i = 0; i < glyphVec.size(); ++i) {
        Glyph &g = glyphVec[i];
        if (g.points.empty()) cipher.push_back(string());
        else {
            cout << "Glyph number " << i << endl;
            //Glyph smoothed = SmoothGlyph(g, 1000, 100, .01);
            Cluster *best = nullptr;
            float score = INFINITY;
            for (int cc = 0; cc < localClusters.size(); ++cc) {
                Cluster &c = localClusters[cc];
                float sc = c.GetShapeContextDiff(g);
                if (sc < score) {
                    best = &c;
                    score = sc;
                }
            }
            if (score > 10) {
                if (best)
                    cout << "Best cluster was " << best->GetLabels() << ", score " << score << endl;
                localClusters.push_back(Cluster());
                best = &localClusters.back();
                best->labels = currlabel++;
                
                //best->labels = actual[i];
                
                cout << "New cluster. " << endl;
            }
            cout << "Labelled as " << best->labels[0] << " (score " << score << ")" << endl;
            best->Insert(glyphVec[i], best->labels[0]);
            cipher.back() += best->labels[0];
        }
    }
    for (int i = 0; i < localClusters.size(); ++i) {
        float best = 0;
        while (best < 12) {
            int j;
            best = ClusterVecMinDist(localClusters, i, j);
            if (best < 12) {
                for (Glyph &g : localClusters[j].glyphs)
                    localClusters[i].Insert(g, localClusters[i].GetLabels()[0]);
                cout << "Merged " << i << ", " << localClusters[i].GetLabels()[0];
                cout << " and " << j << ", " << localClusters[j].GetLabels()[0] << endl;
                localClusters.erase(localClusters.begin() + j);
            }
        }
        cout << "Moving on..." << endl;
    }
    cout << "Trying to solve: " << endl;
    for (string s : cipher)
        cout << s << " ";
    cout << endl;
    
    cipher.clear();
    cipher.push_back(string());
    cout << "Regenerating..." << endl;
    for (int gg = 0; gg < glyphVec.size(); ++gg) {
        Glyph &g = glyphVec[gg];
        if (g.points.empty()) { cipher.push_back(string()); continue; }
        char label;
        float score = INFINITY;
        for (int cc = 0; cc < localClusters.size(); ++cc) {
            Cluster &c = localClusters[cc];
            float sc = c.GetShapeContextDiff(g);
            if (sc < score) {
                label = c.GetLabels()[0];
                score = sc;
            }
        }
        cipher.back() += label;
    }
    
    cout << "Trying to solve: " << endl;
    for (string s : cipher)
        cout << s << " ";
    cout << endl;
    
    string charsUsed;
    for (string s : cipher)
        for (char ch : s)
            if (charsUsed.find(ch) == string::npos)
                charsUsed += ch;
    
    vector<string> cipher2 = cipher;
    for (string &s : cipher2)
        for (char &ch : s)
            ch = 'a' + charsUsed.find(ch);
    cipher = cipher2;
    cout << "Trying to solve: " << endl;
    for (string s : cipher)
        cout << s << " ";
    cout << endl;
    */
    
//    string sentence = "";
//    while (sentence.empty())
//        getline(cin, sentence);
//    cout << sentence << endl;
    
    int n = (int)glyphVec.size();
    
//    if (n != sentence.length()) {
//        cout << "Wrong sentence length" << endl;
//        return;
//    }
//    vector<bool> actualSame(n * n);
//    for (int i = 0; i < n * n; ++i)
//        actualSame[i] = (sentence[i % n] == sentence[i / n]);
    
    vector<ShapeContext> contextVec;
    for (Glyph &g : glyphVec)
        contextVec.push_back(ShapeContext(g));
    vector<Edge> edges, mst;
    cout << "Calculating compare...";
    for (int j = 0; j < n; ++j) {
        for (int i = j + 1; i < n; ++i) {
            float w = compare(contextVec[i], contextVec[j], 5, 12);
            if (w) edges.push_back(Edge(i, j, w));
        }
    }
    cout << "Getting MST" << endl;
    
    Kruskal(n, edges, mst);
    int clusters = 20;
    vector<Glyph> glyphVecCopy(glyphVec);
    float scale = .025f;
    auto drawthunk = [&glyphVecCopy, scale] (vector<Point>& pts) {
        for (int i = 0; i < pts.size(); ++i) {
            glyphVecCopy[i].LineDraw(pts[i], scale);
            Point bottom((float)i / glyphVecCopy.size(), 0);
            glyphVecCopy[i].LineDraw(bottom, 1.f / glyphVecCopy.size());
            Color(0,0,1);
            if (!glyphVecCopy[i].points.empty())
                LineCon(bottom, pts[i], .05f);
            Color(1,1,1);
        }
    };
    //DrawKruskal(n, edges, clusters, drawthunk);
    
    
    InitCipher();
    string result = "";
    int numclusters = 15;
    while (/*result.empty() && */numclusters < 27) {
        string ciphertext = GetCipher(n, mst, numclusters);
        cout << "Cipher, for numclusters " << numclusters << ": " << endl;
        
//        for (int j = 0; j < n; ++j)
//            for (int i = j + 1; i < n; ++i)
//                if (actualSame[j * n + i] && ciphertext[i] != ciphertext[j])
//                    cout << "Oversegmented: real sentence matches two " << sentence[i] << ", " << i << " " << j << endl;
//                else if (!actualSame[j * n + i] && ciphertext[i] == ciphertext[j])
//                    cout << "Wrongly grouped " << i << ": " << sentence[i] << " and "<< j << ": " << sentence[j] << endl;
        
        cout << ciphertext << endl;
        result = Infer(ciphertext);
        ++numclusters;
    }
    cout << "done!!! " << endl;
    
    glyphVec.clear();
    gridentry = 0;

}

void NewGlyphVec() {
    if (currGlyph.points.size()) {
        cout << "Normalizing..." << endl;
        currGlyph.Normalize();
        glyphVec.push_back(currGlyph);
        currGlyph.Empty();
        cout << "Recorded glyph, number " << glyphVec.size() << endl;
    } else {
        glyphVec.push_back(Glyph());
        cout << "Recorded blank glyph, number " << glyphVec.size() << endl;
    }
}

void CheckGrid() {
    int gx = mouse[0] * gwidth;
    int gy = (1 - mouse[1]) * gheight;
    int gbox = gy * gwidth + gx;
    if (gbox != gridentry) {
        NewGlyphVec();
        if (gbox > gridentry + 1)
            NewGlyphVec();
        gridentry = gbox % (gwidth * gheight);
    }
}

static bool WaitingForChar = false;
void KeyCallback(unsigned char key, int x, int y) {

    UpdateMouse(x, y);
    if (WaitingForChar) {
        currlabel = key;
        cout << "Changed label to " << currlabel << endl;
        WaitingForChar = false;
        return;
    }

    switch(key) {
        case 'k': {
            vector<Edge> edges, mst;
            int n = (int)currGlyph.points.size();
            for (int j = 0; j < n; ++j) {
                for (int i = j + 1; i < n; ++i) {
                    edges.push_back(Edge(i, j,
                            Length(currGlyph.points[i] - currGlyph.points[j])));
                }
            }
            Kruskal(n, edges, mst);
            DrawClear();
            Color(1,1,1);
            for (Point &p : currGlyph.points)
                Circle(p, .001f);
            int i;
            for (i = 0; i < mst.size() - 3; ++i)
                LineCon(currGlyph.points[mst[i].a], currGlyph.points[mst[i].b]);
            Color(0,1,0);
            LineCon(currGlyph.points[mst[i].a], currGlyph.points[mst[i].b]);
            ++i;
            Color(0,0,1);
            LineCon(currGlyph.points[mst[i].a], currGlyph.points[mst[i].b]);
            ++i;
            Color(1,0,0);
            LineCon(currGlyph.points[mst[i].a], currGlyph.points[mst[i].b]);
            ++i;
            DrawSwap();
            Color(1,1,1);
            currGlyph.Empty();
        }
            break;
        case '[':
            compGlyph = currGlyph;
            currGlyph.Empty();
            cout << "Saved it" << endl;
            break;
        case ']': {
            vector<float> x(4);
            x[0] = 1;
            vector<vector<float>> sch1, sch2;
            GetSCHistograms(compGlyph, sch1);
            GetSCHistograms(currGlyph, sch2);
            x[1] = LeastCost<vector<float>>(sch1, sch2, ChiSquaredTest, 1.f, false);
            Glyph g1 = SmoothGlyph(compGlyph, 500, 100, .01);
            Glyph g2 = SmoothGlyph(currGlyph, 500, 100, .01);
            x[2] = LeastCost<Point>(g1.points, g2.points, [] (const Point &p, const Point &q) {
                return LengthSq(p - q);
            }, 1.f, false);
            vector<float> c1, c2;
            GetCurvature(g1, c1, true);
            GetCurvature(g2, c2, true);
            x[3] = LeastCostFloat(c1, c2, 1.f, false);
            float result = 0;
            for (int i = 0; i < 4; ++i)
                result += LogRegTheta[i] * x[i];
            cout << "Predicted is " << result << endl;
            currGlyph.Empty();
        }
            break;
        case 'h': {
            vector<vector<float>> X;
            vector<float> y;
            int n = 0;
            string labels = "";
            for (auto &p : glyphs) {
                n += p.second.size();
                labels += string(p.second.size(), p.first);
            }
            for (int j = 0; j < n; ++j) {
                for (int i = j + 1; i < n; ++i) {
                    X.push_back(vector<float>());
                    X.back().push_back(1.f);
                    y.push_back(labels[i] == labels[j] ? 1.f : 0.f);
                }
            }
            cout << "Evaluating shape contexts" << endl;
            EvaluateFeature<vector<float>>(GetSCHistograms, ChiSquaredTest, X);
            LogisticRegression(X, y);
            cout << endl;
            cout << "Evaluating spacial distance" << endl;
            EvaluateFeature<Point>(
                [] (Glyph &g, vector<Point> &v) {
                    Glyph gg = SmoothGlyph(g, 500, 100, .01);
                    v = gg.points;
                }, [] (const Point &p, const Point &q) {
                   return LengthSq(p - q);
                }, X);
            LogisticRegression(X, y);
            cout << endl;
            cout << "Evaluating integrated curvature" << endl;
            EvaluateFeature<float>([](Glyph &g, vector<float> &v) {
                GetCurvature(g, v, true);
            }, [](const float &a, const float &b) {
                return (a-b)*(a-b);
            }, X);
            LogRegTheta = LogisticRegression(X, y);
            cout << endl;
        }
            break;
        case 'z':
            NewGlyphVec();
            break;
        case 'n':
            WriteVectorData();
            break;
        case 'b': {
            string file;
            cin >> file;
            ReadVectorData(file);
        }
            break;
        case 'm': {
  //          string actual;
//            getline(cin, actual);
            SolveCipher();
        }
            break;
        case '`':
            WriteData("data.dat", glyphs);
            break;
        case '.':
            WaitingForChar = true;
            break;
        case ' ':
            NewGlyph();
            break;
        case '\r':
            Analyze(currGlyph, true);
            break;
        case '/':
            cout << Analyze(currGlyph, false);
            break;
        case 'f': {
            Glyph g1 = SmoothGlyph(glyphs[cOne][coneIndex], 800, 200, .01);
            Glyph g2 = SmoothGlyph(glyphs[cTwo][ctwoIndex], 800, 200, .01);
            SimulateGlyph(g1, 2000, 1.4);
            SimulateGlyph(g2, 2000, 1.4);
        }
            break;
        case 'd':
            coneIndex = rng.GetInt((int)glyphs[cOne].size());
            ctwoIndex = rng.GetInt((int)glyphs[cTwo].size());
            ShowMeSpace();
            break;
        case 'q':
            exit(0);
            break;
        case 'c':
            SeedStartTime();
            CreateClusters2();
            cout << "\nAlternative Clustered in " << ElapsedMillis() * .001f << " sec\n\n\n\n\n";
            break;
        case 'v':
            SeedStartTime();
            CreateClusters3();
            cout << "\nClustered in " << ElapsedMillis() * .001f << " sec\n\n\n\n\n";
            break;
        case 'w': {
            string filename = wdir + "screen.ppm";
            WriteScreen(filename);
        }
            break;
        case 'e':
            CalculateError();
            break;
        case '5': {
            string inputfile;
            cout << "Which file to read? ";
            cin >> inputfile;
            ReadData(inputfile, glyphs);
            cout << "\n\n\n\n\n\n"; }
            break;
        case '6': {
            string inputfile;
            cout << "Which file to read? ";
            cin >> inputfile;
            ReadDataSimple(inputfile, glyphs);
            cout << "\n\n\n\n\n\n"; }
            break;
            
        case 'x':
            currGlyph.Empty();
            cout << "Cleared current glyph" << endl;
            DisplayCallback();
            break;
        case 'g':
            glyphs.clear();
            cout << "Cleared the map" << endl;
            break;
        case 'r': {
            cout << "Which chars? " << endl;
            cin >> cOne >> cTwo;
            cout << "Ok." << endl;
            
        }
            break;
        case 's': {
            coneIndex = rng.GetInt((int)glyphs[cOne].size());
            ctwoIndex = rng.GetInt((int)glyphs[cTwo].size());
            ShowMe();
        }
            break;
        case 'a':
            ShowMe();
            break;
        case '1':
            ReadData("data.dat", glyphs);
            break;
        case '3':    // DISPLAYS ALL CHARACTERS YOU'vE WRITTEN
            for (auto &p : glyphs) {
                for (int i = 0; i < p.second.size(); ++i) {
                    cout << "Simulating " << p.first << "...";
                    Glyph g(p.second[i]);
                    g.Normalize();
                    cout << g.points[2] << endl;
                    SimulateGlyph(g, 1000, 1.5);
                    cout << "done" << endl;
                }
            }
            break;
    }
}

void KeyUpCallback(unsigned char key, int x, int y) {
    UpdateMouse(x, y);
    switch(key) {
            
    }
}

void SpecialKeyCallback(int key, int x, int y) {
    UpdateMouse(x, y);
    //arrows are GLUT_KEY_LEFT, etc.
    switch(key) {
            
    }
}

void SpecialKeyUpCallback(int key, int x, int y) {
    UpdateMouse(x, y);
    //arrows are GLUT_KEY_LEFT, etc.
    switch(key) {
            
    }
}

// For live recognition, turned off.
void IdleCallback() {
    if (timerValid && timer.Check()) {
        cout << Analyze(currGlyph, false);
        DrawClear();
        DrawSwap();
    }
}

void ReshapeCallback(int x, int y) {
    width = x;
    height = y;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, 1, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void LineGrid(int m, int n) {
    Color(1,1,1);
    for (int i = 1; i < m ; ++i)
        LineCon(0, (float)i / m, 1, (float)i / m);
    for (int j = 1; j < n ; ++j)
        LineCon((float)j / n, 0, (float)j / n, 1);
}

void DisplayCallback() {
    glClear(GL_COLOR_BUFFER_BIT);
    
    currGlyph.Draw();
    
    if (gridOn)
        LineGrid(gheight, gwidth);
    
    glutSwapBuffers();
}


int main(int argc, const char * argv[])
{
    Init();
    glutInit(&argc, const_cast<char **>(argv));
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowPosition(50, 50);
    glutInitWindowSize(width, height);
    glutCreateWindow("LiveWrite");
    //glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_MULTISAMPLE);
    
    
    glutMouseFunc(MouseCallback);
    glutMotionFunc(MouseMotionCallback);
    glutPassiveMotionFunc(PassiveMotionCallback);
    glutKeyboardFunc(KeyCallback);
    glutKeyboardUpFunc(KeyUpCallback);
    glutSpecialFunc(SpecialKeyCallback);
    glutSpecialUpFunc(SpecialKeyUpCallback);
    //glutIdleFunc(IdleCallback);
    glutReshapeFunc(ReshapeCallback);
    glutDisplayFunc(DisplayCallback);
    
    glutMainLoop();
    
    return 0;
}

