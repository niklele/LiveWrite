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
static char currlabel = 'A' - 1;
static vector<TrainEx> glyphs;

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
    if (state == GLUT_UP) glyphs.back().g.MouseUp();
    if (state == GLUT_DOWN) glyphs.back().g.MouseDown();
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
    
    glyphs.back().g.AddPoint(mouse, ElapsedMillis());
    
    DisplayCallback();
    
}

void PassiveMotionCallback(int x, int y) {
    UpdateMouse(x,y);
    
}

static int gindex = 0;
void KeyCallback(unsigned char key, int x, int y) {
    UpdateMouse(x, y);
    switch(key) {
        case 'r':
            currlabel++;
        case ' ':

            glyphs.push_back(TrainEx({Glyph(), currlabel}));
            cout << "Started new glyph with label " << currlabel << endl;
            break;
        case 'd':
            DrawClear();
            
            DrawSwap();
            gindex++;
            break;
        case 'q':
            exit(0);
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
    if (glyphs.size()) glyphs.back().g.Draw();
    
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

