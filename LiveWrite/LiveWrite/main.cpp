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
//static float start = 0;


void DisplayCallback();

static int highlight = 0;
static int bezhigh = 1;

float checker = 1.f;
struct TrainEx {
    Glyph g;
    char label;
};
static char currlabel = 'A';
static map<char, vector<Glyph>> glyphs;
static Glyph currGlyph;

static bool CirclesOn = true;


Point RandDir(Point dir, float scale = 1.f) {
    Point ret(rng.GetNudge(), rng.GetNudge());
    while (Dot(ret, dir) < 0) {
        ret = Point(rng.GetNudge(), rng.GetNudge());
    }
    ret *= scale;
    return ret;
}

void FillRand(float *f, int num) {
    for (int i = 0; i < num; ++i)
        f[i] = rng();
}

void RandCol(float col[], float low = 0.f) {
    FillRand(col, 3);
    float mag = col[0]*col[0] + col[1]*col[1] + col[2]*col[2];
    while (mag < low) {
        FillRand(col, 3);
        mag = col[0]*col[0] + col[1]*col[1] + col[2]*col[2];
    }
}

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

void MouseCallback(int button, int state, int x, int y) {
    //button = GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON
    //state = GLUT_UP or GLUT_DOWN
    UpdateMouse(x,y);
    if (state == GLUT_UP) currGlyph.MouseUp();
    if (state == GLUT_DOWN) currGlyph.MouseDown();
}

vector<Point> mice;
vector<float> thetas(100, 0);

void MouseMotionCallback(int x, int y) {
    /*
    micePix[3 * ((700-y) * 700 + x) + 0] = 255;
    Point delta = MouseDelta(x, y);
    PostElapsed();
    if (Length(delta) < .01) return;
    cout << "(actually used..)\n";
    Point prev = mouse - mouseOld;
    float cosine = Dot(prev, delta) / Length(prev) / Length(delta);
    float thet = acos(cosine) / 3.1415f;
    thetas.push_back(Bearing(prev, delta));
    
    micePix[3 * ((700-y) * 700 + x) + 0] = 0;
    micePix[3 * ((700-y) * 700 + x) + 1] = 255;
     */
    UpdateMouse(x,y);
    
    currGlyph.AddPoint(mouse, ElapsedMillis());
    
    DisplayCallback();
    
}

void PassiveMotionCallback(int x, int y) {
    UpdateMouse(x,y);
    
}

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

void SimulateGlyph(Glyph &g, float totaltime = 1000.f, float extend = 1.5f) {
    SeedStartTime();
    float time;
    while ((time = ElapsedMillis()) <= totaltime * extend) {
        DrawClear();
        for (int i = 0; i < g.times.size(); ++i) {
            if (totaltime * g.times[i] < time) Circle(g.points[i], .001);
        }
        DrawSwap();
        //usleep(10000);
    }
}

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
        case 'd':
            DrawClear();
            
            DrawSwap();
            gindex++;
            break;
        case 'q':
            exit(0);
            break;
        case '1':
            ReadData("data.dat", glyphs);
            break;
        case '2':
            for (auto &p : glyphs) {
                Glyph composite;
                for (Glyph &g : p.second) {
                    composite.points.insert(composite.points.end(), g.points.begin(), g.points.end());
                    composite.times.insert(composite.times.end(), g.times.begin(), g.times.end());

                }
                SimulateGlyph(composite);
            }
            break;
        case '3':
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

void IdleCallback() {
    //DisplayCallback();
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

void DrawThetas(int num) {
    int subdivs = 100;
    int avg = num;
    for (int i = max((int)(thetas.size() - subdivs), 1); i < thetas.size(); ++i) {
        int off = i - max((int)(thetas.size() - subdivs), 1);
        float oldy = 0, newy = 0, weights = 0;
        for (int j = 0; j < avg; ++j) {
            float w = (avg/2 - abs(j - avg/2));
            oldy += thetas[i - 1 - j] * w;
            newy += thetas[i - j] * w;
            weights += w;
        }
        oldy /= weights;
        newy /= weights;
        LineCon((off - 1.f) / subdivs, oldy + .5f, (off - 0.f) / subdivs, newy + .5f);
    }
}

void DisplayCallback() {
    glClear(GL_COLOR_BUFFER_BIT);
    
    //glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, micePix);
    currGlyph.Draw();
    
    //for (int i = 0; i < mice.size(); ++i) Circle(mice[i], .001);
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

