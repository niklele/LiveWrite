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


static float mouseX, mouseY;
static Point mouse, mouseOld;
const static int NRES = 100;
const static float STDEV = .03f;
const static int UPSCALE = 500;

void DisplayCallback();

float checker = 1.f;
struct TrainEx {
    Glyph g;
    char label;
};
static char currlabel = 'A';
static map<char, vector<Glyph>> glyphs;
static vector<Glyph> glyphVec;
static Glyph currGlyph;
static map<char, vector<float>> phivals;
static map<char, vector<vector<float>>> curvint;

static int mousecount = 0;
void Init() {
    mousecount = 0;
    SeedStartTime();
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


static unsigned char micePix[700 * 700 * 3];
char Analyze(Glyph &g, bool verbose = false);
static Timer timer;
static bool timerValid = false;
void MouseCallback(int button, int state, int x, int y) {
    //button = GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON
    //state = GLUT_UP or GLUT_DOWN
    UpdateMouse(x,y);
    if (state == GLUT_UP) {
        currGlyph.MouseUp();
        timer.SetTimer(500.f);
        //timerValid = true;
    }
    if (state == GLUT_DOWN) {
        currGlyph.AddPoint(mouse, ElapsedMillis());
        currGlyph.MouseDown();
        timerValid = false;
    }
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

const static string wdir = "/Users/ben/Documents/projects/LiveWrite/LiveWrite/LiveWrite/";
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
    glyphVec.clear();
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
            if (score > 100) {
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

/*
void KMeans(vector<Glyph> &glyphs, vector<Cluster> &clusters, int iters, int numClusters) {
    if (glyphs.size() < numClusters) {
        cout << "Can't have more clusters than data points " << endl;
        return;
    }
    clusters.clear();
    for (int i = 0; i < numClusters; ++i) {
        clusters.push_back(Cluster());
        clusters.back().Insert(glyphs[i], i);
    }
    for (int iter = 0; iter < iters; ++iter) {
        vector<int> clusterInds(glyphs.size());
        for (int gi = 0; gi < glyphs.size(); ++gi) {
            Glyph &g = glyphs[gi];
            Glyph smoothed = SmoothGlyph(g, UPSCALE, NRES, STDEV);
            float score = INFINITY;
            for (int ci = 0; ci < clusters.size(); ++ci) {
                float sc = clusters[ci].GetDistanceThreaded(smoothed);
                if (sc < score) {
                    score = sc;
                    clusterInds[gi] = ci;
                }
            }
        }
        for (int ci = 0; ci < clusters.size(); ++ci) {
            clusters[ci].Clear();
            for (int gi = 0; gi < glyphs.size(); ++gi) {
                if (clusterInds[gi] == ci)
                    clusters[ci].Insert(glyphs[gi], 0, false);
            }
        }
    }
}
*/

static Point dumbpoint;

char Analyze(Glyph &g, bool verbose) {
    if (g.points.size()) {
        g.Normalize();
        dumbpoint = g.points[10]; // Debugging crap
        Glyph currInt = SmoothGlyph(g, UPSCALE, NRES, STDEV);
        dumbpoint = currInt.points[10]; // Debugging crap
        cout << dumbpoint;
        g.Empty();
        if (clusters.empty()) return 0;
        map<float, char> scores;
        for (Cluster &c : clusters) {
            float curr = c.GetDistanceThreaded(currInt, false);
            scores[curr] = c.GetLabels()[0];
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
        return scores.begin()->second;
        
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

void DrawMatrix(vector<float> &vals, int n, float max) {
    for (int x = 0; x < n * n; ++x) {
        float col = vals[x] / max;
        Color(col,col,col);
        DrawBox((float)(x % n) / n, (float)(x / n) / n, 1.f / n, 1.f / n);
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

void LogisticRegression(vector<vector<float>> &X, vector<float> &y, int iters = 10000) {
    float alpha = .01f;
    int m = (int)X.size();
    cout << X[100][3] << endl;
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
    }
    for (int j = 0; j < n; ++j)
        theta[j] /= theta[0];
    cout << "Logistic regression produced the following theta parameters" << endl;
    cout << "in " << iters << " iterations and using alpha value " << alpha << ": " << endl;
    cout << "(";
    for (int i = 0; i < n; ++i)
        cout << " " << theta[i] << " ";
    cout << ")" << endl;
    int fails = 0;
    for (int i = 0; i < m; ++i) {
        float hyp = 0.f;
        for (int j = 0; j < n; ++j)
            hyp += theta[j] * X[i][j];
        hyp = (Sigmoid(hyp) >= .5f ? 1.f : 0.f);
        if (hyp != y[i])
            fails++;
            //cout << "Failure on data point " << i << ", supposed to be " << y[i] << endl;
    }
    cout << "Fails on " << fails << " examples" << endl;
}

void SolveCipher() {
    
    vector<Cluster> localClusters;
    char currlabel = 'a';
    vector<string> cipher;
    cipher.push_back(string());
    vector<Glyph> smoothed;
    for (int i = 0; i < glyphVec.size(); ++i)
        if (!glyphVec[i].points.empty())
            smoothed.push_back(SmoothGlyph(glyphVec[i], 500, 100, .01));
    
    int n = (int)smoothed.size();
    vector<float> lcp(n * n, 0);
    string actual = "thepoisonfromablackwidowspiderisaboutfifteentimesmorepotent";
    float highest = 0.f;
    vector<vector<float>> features;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            vector<float> feat;
            lcp[j * n + i] = -GetGlyphDistance(smoothed[i], smoothed[j], feat);
            features.push_back(feat);
            if (lcp[j * n + i] > highest) {
                highest = lcp[j * n + i];
                cout << i << " = " << actual[i] << ", " << j
                    << " = " << actual[j] << ": " << lcp[j * n + i] << endl;
            }
        }
        DrawClear();
        DrawMatrix(lcp, n, highest);
        HighlightMatches(actual);
        DrawSwap();
    }
    DrawClear();
    DrawMatrix(lcp, n, highest);
    Color(0,1,0);
    int nstr = (int)actual.length();
    if (nstr != n) {
        cout << "Length mismatch, abort" << endl;
        return;
    }
    map<float, pair<int,int>> scorevalues;
    int needed = 0;
    vector<float> classes;
    for (int j = 0; j < nstr; ++j) {
        for (int i = 0; i < nstr; ++i) {
            classes.push_back(actual[i] == actual[j] ? 1.f : 0.f);
            if (actual[i] == actual[j]) {
                Circle((i + .5f) / nstr, (j + .5f) / nstr, .015);
                scorevalues[lcp[j * n + i]] = {i, j};
                ++needed;
            }
        }
    }
    for (auto &p : scorevalues) {
        cout << "(" << p.second.first << ", " << p.second.second << "): " << p.first << endl;
    }
    int incorrect = 0, missing = 0;
    auto iter = scorevalues.end();
    iter--;
    float highestgood = iter->first;
    highestgood = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (lcp[j * n + i] <= highestgood && actual[i] != actual[j]) {
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
    
    LogisticRegression(features, classes);
    char dummy;
    cin >> dummy;
    for (int i = 0; i < glyphVec.size(); ++i) {
        Glyph &g = glyphVec[i];
        if (g.points.empty()) cipher.push_back(string());
        else {
            cout << "Glyph number " << i << endl;
            Glyph smoothed = SmoothGlyph(g, 1000, 100, .01);//SmoothGlyph(g, UPSCALE, NRES, STDEV);
            Cluster *best = nullptr;
            float score = INFINITY;
            for (Cluster &c : localClusters) {
                float sc = c.GetDistanceThreaded(smoothed);
                if (sc < score) {
                    best = &c;
                    score = sc;
                }
            }
            if (score > 150) {
                if (best)
                    cout << "Best cluster was " << best->GetLabels() << ", score " << score << endl;
                localClusters.push_back(Cluster());
                best = &localClusters.back();
                best->labels = currlabel++;
                cout << "New cluster. " << endl;
            }
            cout << "Labelled as " << best->labels[0] << " (score " << score << ")" << endl;
            best->Insert(glyphVec[i], best->labels[0]);
            cipher.back() += best->labels[0];
        }
    }
    cout << "Trying to solve: " << endl;
    for (string s : cipher)
        cout << s << " ";
    cout << endl;
    glyphVec.clear();
    
    map<string, int> freqs;
    map<string, int> ngrams;
    int ng = 3;
    
    ReadCorpus("/Users/ben/Documents/projects/corpus", freqs);
    cout << "Found " << freqs.size() << " words" << endl;
    GetNGrams(freqs, ngrams, ng);
    
    RunInference(50, 5000, freqs, ngrams, cipher, ng);

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
        case 'z':
            if (currGlyph.points.size()) {
                cout << "Normalizing..." << endl;
                currGlyph.Normalize();
                DrawClear();
                currGlyph.Draw();
                DrawSwap();
                glyphVec.push_back(currGlyph);
                currGlyph.Empty();
                cout << "Recorded glyph, number " << glyphVec.size() << endl;
            } else {
                glyphVec.push_back(Glyph());
                cout << "Recorded blank glyph, number " << glyphVec.size() << endl;
            }
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
        case 'm':
            SolveCipher();
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
            CreateClusters(true);
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

void DisplayCallback() {
    glClear(GL_COLOR_BUFFER_BIT);
    
    currGlyph.Draw();
    
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

