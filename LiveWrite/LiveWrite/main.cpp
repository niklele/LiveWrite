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


static float mouseX, mouseY;
static Point mouse, mouseOld;
const static int NRES = 100;

void DisplayCallback();

float checker = 1.f;
struct TrainEx {
    Glyph g;
    char label;
};
static char currlabel = 'A';
static map<char, vector<Glyph>> glyphs;
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
void Analyze(bool verbose = false);
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
    outfile << "Number_of_examples: " << nex << endl;
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
}

void GetCurvature(Glyph &g, vector<float> &curv, bool accum = false) {
    if (g.times.size() < 2) return;
    Point prev = g.points[1] - g.points[0];
    float acc = 0;

    for (int i = 2; i < g.times.size(); ++i) {
        Point next = g.points[i] - g.points[i-1];
        acc = Bearing(prev, next) / M_PI / 2.f + (accum ? acc : 0);
        curv.push_back(acc);
        prev = next;
    }
}

void GetCurvReverse(Glyph &g, vector<float> &curv, bool accum = false) {
    if (g.times.size() < 2) return;
    int last = (int)g.points.size() - 1;
    Point prev = g.points[last-1] - g.points[last];
    float acc = 0;
    for (int i = last - 1; i >= 1; --i) {
        Point next = g.points[i-1] - g.points[i];
        acc = Bearing(prev, next) / M_PI / 2.f + (accum ? acc : 0);
        curv.push_back(acc);
        prev = next;
    }
}

/*
void DrawCurvature(Glyph &g, float time) {
    if (g.times.size() < 2) return;
    LineCon(0, .5, 1, .5);
    Color(1,0,0);
    Point prev = g.points[1] - g.points[0];
    Point graphpt(0,.5);
    for (int i = 2; i < g.times.size(); ++i) {
        if (g.times[i] >= time) Color(0,1,0);
        Point next = g.points[i] - g.points[i-1];
        Point ngpt(g.times[i], .5f + Bearing(prev, next) / M_PI * .5f);
        LineCon(graphpt, ngpt);
        graphpt = ngpt;
        prev = next;
    }
    Color(1,1,1);
}
 */


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
    for (int i = 0; i < g.times.size(); ++i) {
        if (g.times[i] < t) Circle(g.points[i], .001);
    }
}


class Cluster {
public:
    float GetDistance(Glyph &g, bool verbose = false) {
        //verbose = true;
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
    void Insert(Glyph &g, char label) {
        glyphs.push_back(g);
        labels += label;
        composite.points.insert(composite.points.end(), g.points.begin(), g.points.end());
        composite.times.insert(composite.times.end(), g.times.begin(), g.times.end());
        average = ReduceGlyph(composite, NRES, .01);
        cf.clear(); cr.clear(); cif.clear(); cir.clear();
        GetCurvature(average, cf, false);
        GetCurvReverse(average, cr, false);
        GetCurvature(average, cif, true);
        GetCurvReverse(average, cir, true);
    }
    string GetLabels() {
        return labels;
    }
    float AverageScore() {
        float ret = 0;
        for (Glyph &g : glyphs) ret += GetDistance(g);
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
    int nexamples;
    infile >> dummy >> nexamples;
    cout << "Reading in " << nexamples << " data examples" << endl;
    for (int i = 0; i < nexamples; ++i) {
        TrainEx curr;
        int npoints;
        infile >> dummy >> curr.label >> dummy >> npoints;
        cout << "Label " << curr.label << " with " << npoints << " points" << endl;
        infile >> dummy >> curr.g.width >> dummy >> curr.g.duration;
        infile >> dummy >> curr.g.centroid[0] >> curr.g.centroid[1];
        for (int j = 0; j < npoints; ++j) {
            float t;
            Point p;
            infile >> t >> p[0] >> p[1];
            curr.g.times.push_back(t);
            curr.g.points.push_back(p);
        }
        data[curr.label].push_back(curr.g);
    }
}

/*
void ReadProcessData(string filename, map<char, vector<Glyph>> &data) {
    string path = wdir + filename;
    ifstream infile(path);
    string dummy;
    int nexamples;
    infile >> dummy >> nexamples;
    cout << "Reading in " << nexamples << " data examples" << endl;
    for (int i = 0; i < nexamples; ++i) {
        TrainEx curr;
        int npoints;
        infile >> dummy >> curr.label >> dummy >> npoints;
        cout << "Label " << curr.label << " with " << npoints << " points" << endl;
        infile >> dummy >> curr.g.width >> dummy >> curr.g.duration;
        infile >> dummy >> curr.g.centroid[0] >> curr.g.centroid[1];
        for (int j = 0; j < npoints; ++j) {
            float t;
            Point p;
            infile >> t >> p[0] >> p[1];
            curr.g.times.push_back(t);
            curr.g.points.push_back(p);
        }
        
        Glyph actual = InterpGlyph(curr.g, 500);
        actual = ReduceGlyph(actual, 100, .01);
        data[curr.label].push_back(actual);
        Cluster *best = nullptr;
        float score = INFINITY;
        for (Cluster &c : clusters) {
            if (c.GetLabels()[0] != curr.label) continue;
            float sc = c.GetDistance(actual);
            if (sc < score) {
                best = &c;
                score = sc;
            }
        }
        if (score > 80) {
            if (best)
                cout << "Best cluster was " << best->GetLabels() << ", score " << score << endl;
            clusters.push_back(Cluster());
            best = &clusters.back();
            cout << "New cluster. ";
        }
        cout << "Adding a \'" << curr.label << "\' to cluster " << best->GetLabels() << endl;
        best->Insert(actual, curr.label);
        //data[curr.label].push_back(SmoothGlyph(curr.g, 1000, 100, .05));
    }
    
    for (int i = 0; i < clusters.size(); ++i) {
        if (clusters[i].glyphs.size() == 1) clusters.erase(clusters.begin() + i--);
    }
    
}
*/

void CreateClusters() {
    clusters.clear();
    for (auto &p : glyphs) {
        char label = p.first;
        vector<Glyph> &vec = p.second;
        for (Glyph &g : vec) {
            Glyph smoothed = InterpGlyph(g, 500);
            smoothed = ReduceGlyph(smoothed, 100, .01);
            Cluster *best = nullptr;
            float score = INFINITY;
            for (Cluster &c : clusters) {
                if (c.GetLabels()[0] != label) continue;
                float sc = c.GetDistance(smoothed);
                if (sc < score) {
                    best = &c;
                    score = sc;
                }
            }
            if (score > 80) {
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
    for (int i = 0; i < clusters.size(); ++i) {
        cout << "#" << i << " (" << clusters[i].AverageScore() << "): " << clusters[i].GetLabels() << endl;
    }
}

/*
void PostPhiVals(vector<float> &phis, vector<int> &counts, int nres) {
    int denom = 0;
    phis.resize(nres*nres);
    for (int &i : counts) denom += i;
    for (int i = 0; i < nres * nres; ++i) {
        phis[i] = log((float)counts[i]) - log((float)denom);
        //cout << "putting in " << phis[i] << endl;
    }
}
 */

static Point dumbpoint;

void Analyze(bool verbose) { // LIEZ
    if (currGlyph.points.size()) {
        currGlyph.Normalize();
        Glyph currInt = SmoothGlyph(currGlyph, 500, 100, .05);
        dumbpoint = currInt.points[10]; // Debugging crap
        //cout << dumbpoint << endl;
        currGlyph.Empty();
        if (clusters.empty()) return;
        map<float, char> scores;
        for (Cluster &c : clusters) {
            float curr = c.GetDistance(currInt, false);
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
        else cout << scores.begin()->second;
    }
}


static char cOne, cTwo;
static int coneIndex, ctwoIndex;
void ShowMe() {
    DrawClear();
    LeastCostFloat(curvint[cOne][coneIndex], curvint[cTwo][ctwoIndex], 1.3f);
    PlotCurvature(glyphs[cOne][coneIndex], 1, true);
    PlotCurvature(glyphs[cTwo][ctwoIndex], 0, true);
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

static bool WaitingForChar = false;
static int gindex = 0;
void KeyCallback(unsigned char key, int x, int y) {

    UpdateMouse(x, y);
    if (WaitingForChar) {
        currlabel = key;
        cout << "Changed label to " << currlabel << endl;
        WaitingForChar = false;
        return;
    }

    switch(key) {
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
            Analyze(true);
            break;
        case '/':
            Analyze(false);
            break;
        case 'f':
            SimulateGlyph(glyphs[cOne][coneIndex], 2000, 1.4);
            SimulateGlyph(glyphs[cTwo][ctwoIndex], 2000, 1.4);
            break;
        case 'd':
            coneIndex = rng.GetInt(glyphs[cOne].size());
            ctwoIndex = rng.GetInt(glyphs[cTwo].size());
            ShowMeSpace();
            break;
        case 'q':
            exit(0);
            break;
        case '4':
            CreateClusters();
            cout << "\n\n\n\n\n\n";
            break;
        case '5': {
            string inputfile;
            cout << "Which file to read? ";
            cin >> inputfile;
            ReadData(inputfile, glyphs);
            cout << "\n\n\n\n\n\n"; }
            break;
        case 'x':
            currGlyph.Empty();
            DisplayCallback();
            break;
            
        case 'r': {
            cout << "Which chars? " << endl;
            cin >> cOne >> cTwo;
            cout << "Ok." << endl;
            
        }
            break;
        case 's': {
            coneIndex = rng.GetInt(glyphs[cOne].size());
            ctwoIndex = rng.GetInt(glyphs[cTwo].size());
            ShowMe();
        }
            break;
        case '1':
            ReadData("data.dat", glyphs);
            break;
            
            /*
        case '2': {
            phivals.clear();
            for (auto &p : glyphs) {
                Glyph composite;
                int nres = 64;
                vector<int> curvs(nres * nres, 1);
                for (Glyph &g : p.second) {
                    composite.points.insert(composite.points.end(), g.points.begin(), g.points.end());
                    composite.times.insert(composite.times.end(), g.times.begin(), g.times.end());
                    DumpCurvatures(g, curvs, nres);
                    curvint[p.first].push_back(vector<float>());
                    GetCurvReverse(g, curvint[p.first].back());
                }
                PostPhiVals(phivals[p.first], curvs, nres);
                cout << "Posting phis for " << p.first << endl;
                Glyph reduced = ReduceGlyph(composite, 100, .05);
                
                //SimulateGlyph(reduced, 2000, 1.5);
                //SimulateBigGlyph(composite, 3000, 1.5, curvs, nres);
            }
            cout << "Done posting phis" << endl;

        }
            break;*/
        case '3':    // DISPLAYS ALL CHARACTERS YOU'vE WRITTEN
            for (auto &p : glyphs) {
                for (int i = 0; i < p.second.size(); ++i) {
                    cout << "Simulating " << p.first << "...";
                    SimulateGlyph(p.second[i], 1000, 1.5);
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
        Analyze(false);
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

