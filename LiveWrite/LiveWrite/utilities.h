//
//  utilities.h
//  Games
//
//  Created by ben on 7/2/13.
//  Copyright (c) 2013 Ben Mildenhall. All rights reserved.
//

#ifndef Games_utilities_h
#define Games_utilities_h

#include <GLUT/glut.h>
#include <sys/time.h>
#include <math.h>
#include <random>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <thread>
#include <vector>
#include <map>
#include <iomanip>
#include <cstring>
#include <execinfo.h>
#include <functional>

using namespace std;

static int width = 800;
static int height = width;
const static float EPSILON = .000001f;

enum Key {
    UP = 1 << 0,
    DOWN = 1 << 1,
    LEFT = 1 << 2,
    RIGHT = 1 << 3
};

static long StartSeconds;

inline void DrawClear() { glClear(GL_COLOR_BUFFER_BIT); }
inline void DrawSwap()  { glutSwapBuffers(); }

float max2f(float a, float b) {
    return std::max(a, b);
}

float max3f(float a, float b, float c) {
    return std::max(a, std::max(b, c));
}

float min2f(float a, float b) {
    return std::min(a,b);
}

float min3f(float a, float b, float c) {
    return std::min(a, std::min(b, c));
}

float Luminance(unsigned char *c) {
    return (c[0]*.2126f+c[1]*.7152f+c[2]*.0722f) / 255.f;
}

float Luminance(float *c) {
    return c[0]*.2126f+c[1]*.7152f+c[2]*.0722f;
}

void WriteScreen(string filename) {
    glReadBuffer(GL_FRONT_LEFT);
    ofstream file;
    file.open(filename);
    file << "P6" << endl << width << " ";
    file << height << endl << "255" << endl;
    char *pix = new char[3 * width * height];
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pix);
    for (int h = 0; h < height; ++h)
        file.write(pix + 3 * width * (height - 1 - h), 3 * width);
    delete pix;
    file.close();
}

// Writes, assuming that pix stores the 0-255 RGB values of each pixel
int WriteImage(string filename, int width, int height, char *pix) {
    ofstream file;
    file.open(filename);
    if (file.fail()) {
        std::cout << "Could not write image file" << filename << endl;
        return -1;
    }
    file << "P6" << endl << width << " ";
    file << height << endl << "255" << endl;
    file.write(pix, 3 * width * height);
    file.close();
    return 0;
}

int ReadImage(string filename, int &width, int &height, char *pix) {
    std::ifstream file(filename.c_str());
    if (file.fail()) {
        std::cout << "Could not read image file " << filename << endl;
        return -1;
    }
    cout << "Reading image file " << filename << endl;
    char first[20];
    file.getline(first, 19);
    file >> width >> height;
    cout << "Width " << width << ", Height " << height << endl;
    int pixsize;
    file >> pixsize;
    cout << "Color range " << pixsize << endl;
    pix = (char *)malloc(3*width*height);
    file.read(pix, 1);
    file.read(pix, 3*width*height);
    file.close();
    return 0;
}

static mutex coutlock;
void WriteLock() { coutlock.lock(); }
void WriteUnlock() { cout << flush; coutlock.unlock(); }


inline float Dist2d(float a[2], float b[2]) {
    return sqrtf((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]));
}

void WriteStr(const char *mystr, float x = 0, float y = 0) {
    glRasterPos2f(x, y);
    for (int i = 0; mystr[i]; i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, mystr[i]);
}

void WriteStrSm(const char *mystr, float x = 0, float y = 0) {
    glRasterPos2f(x, y);
    for (int i = 0; mystr[i]; i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, mystr[i]);
}

void WriteStr(std::string str, float x, float y) {
    WriteStr(str.c_str(), x, y);
}

void WriteStrSm(std::string str, float x, float y) {
    WriteStrSm(str.c_str(), x, y);
}


bool OnScreen(float x[2], float eps = 0) {
    return x[0] >= -eps && x[1] <= 1+eps && x[0] <= 1+eps && x[1] >= -eps;
}

void PutOnScreen(float x[2]) {
    for (int i = 0; i < 2; i++) {
        float &n = x[i];
        if (n < 0) n =0;//+= 1;
        if (n > 1) n =1;//-= 1;
    }
}

float Clamp(float x) { return std::min(std::max(x, 0.f), 1.f); }

void Color(float r, float g, float b) {
    r = Clamp(r); g = Clamp(g); b = Clamp(b);
    glColor3f(r, g, b);
}

void Color(float *c) {
    glColor3fv(c);
}

void Circle(float x, float y, float rad = .01) {
    glBegin(GL_LINE_LOOP);
    int num = 32;
    for (int i = 0; i < num; ++i) {
        glVertex2f(x+rad * cos(2 * M_PI * i / num),
                   y+rad * sin(2 * M_PI * i / num));
    }
    glEnd();
}


void Circle(float x[], float rad = .01) {
    Circle(x[0],x[1], rad);
}


void SeedStartTime() {
    timeval time;
    gettimeofday(&time, NULL);
    StartSeconds = time.tv_sec;
}

static long PauseTime;
void Pause() {
    timeval time;
    gettimeofday(&time, NULL);
    PauseTime = time.tv_sec;
}

void Unpause() {
    StartSeconds = PauseTime;
}

float ElapsedMillis(float initial = 0) {
    timeval time;
    gettimeofday(&time, NULL);
    return (1000.f * (time.tv_sec - StartSeconds)
            + .001f * time.tv_usec) - initial;
}

void WaitMillis(float wait) {
    float initial = ElapsedMillis();
    while (wait > ElapsedMillis(initial)) ;
    return;
}

void PostElapsed() {
    cout << ElapsedMillis() << " millisec elapsed.\n";
    SeedStartTime();
}

class Timer {
public:
    void SetTimer(float time) {
        set = ElapsedMillis();
        wait = time;
    }
    bool Check() {
        return ElapsedMillis(set) >= wait;
    }
    void Wait(float time) {
        SetTimer(time);
        while (!Check()) ;
    }
private:
    float set, wait;
};

class RNG {
public:
    RNG() {
        timeval time;
        gettimeofday(&time, NULL);
        generator.seed((unsigned int)time.tv_sec);
    }
    float GetUni(float m = 1.f) { return m * udist(generator); }
    int GetInt(int max) { return (int)GetUni((float)max); }
    float GetNudge(float m = 1.f) { return m * 2.f * (GetUni() - .5f); }
    float Trial(float p) { return GetUni() < p; }
    float GetNormal() { return .33f * ndist(generator); }
    int PMO() { return GetInt(3) - 1; }
    float operator()(float p) { return Trial(p); }
    float operator()() { return GetUni(); }
private:
    std::default_random_engine generator;
    std::uniform_real_distribution<float> udist;
    std::normal_distribution<float> ndist;
};

static RNG rng;


class Point {
public:
    Point() { data[0] = data[1] = 0.f; }
    Point(float x, float y) {
        data[0] = x;
        data[1] = y;
    }
    Point(float *v) {
        data[0] = v[0];
        data[1] = v[1];
    }
    Point(const Point& p) {
        data[0] = p[0];
        data[1] = p[1];
    }
    string Print() {
        stringstream ret;
        ret.precision(3);
        ret << "(";
        for (int i = 0; i < 2; ++i) {
            ret << data[i];
            if (i < 2 - 1) ret << ",";
            else ret << ")";
        }
        return ret.str();
    }
    Point &operator=(const Point &p) {
        data[0] = p[0];
        data[1] = p[1];
        return *this;
    }
    Point &operator+=(const Point &p) {
        data[0] += p[0];
        data[1] += p[1];
        return *this;
    }
    Point operator+(const Point &p) const {
        Point ret(p);
        ret += *this;
        return ret;
    }
    Point operator-() const {
        Point ret(*this);
        ret[0] = -ret[0];
        ret[1] = -ret[1];
        return ret;
    }
    Point &operator-=(const Point &p) {
        *this += -p;
        return *this;
    }
    Point operator-(const Point &p) const {
        return *this + -p;
    }
    Point &operator*=(const float f) {
        data[0] *= f;
        data[1] *= f;
        return *this;
    }
    Point operator*(const float f) const {
        Point ret(*this);
        ret *= f;
        return ret;
    }
    Point operator*(const Point &p) {
        Point ret(*this);
        ret[0] *= p[0];
        ret[1] *= p[1];
        return ret;
    }
    Point operator/(const Point &p) {
        Point ret(*this);
        ret[0] /= p[0];
        ret[1] /= p[1];
        return ret;
    }
    Point operator/=(const float f) {
        return operator*=(1.f / f);
    }
    Point operator/(const float f) const {
        Point ret(*this);
        ret /= f;
        return ret;
    }
    bool operator==(const Point &p) {
        return p[0]==data[0] && p[1]==data[1];
    }
    friend ostream &operator<<(ostream &s, const Point &p) {
        s << "(" << p[0] << ", " << p[1] << ")";
        return s;
    }
    inline const float &operator[](int i) const { return data[i]; }
    inline float &operator[](int i) { return data[i]; }

    float data[2];
};

Point DirPt(float theta) {
    return Point(cos(theta), sin(theta));
}

bool OnScreen(Point p, float eps = 0) {
    return OnScreen(p.data, eps);
}


inline float Dot(const Point &p, const Point &q) {
    return p[0]*q[0] + p[1]*q[1];
}

inline float LengthSq(const Point &p) {
    return Dot(p,p);
}

inline float Length(const Point &p) {
    return sqrtf(LengthSq(p));
}

inline float Det(const Point &u, const Point &v) {
    return u[0]*v[1] - u[1]*v[0];
}

inline float pCos(const Point &u, const Point &v) {
    float ret = Dot(u, v) / (Length(u) * Length(v));
    ret = std::max(-1.f, std::min(1.f,ret));
    return ret;
}

inline float pSin(const Point &u, const Point &v) {
    float ret = Det(u, v) / (Length(u) * Length(v));
    ret = std::max(-1.f, std::min(1.f,ret));
    return ret;
}

inline float Bearing(const Point &u, const Point &v) {
    if (pCos(u,v) >= 0) {
        return asin(pSin(u,v));
    } else {
        return M_PI - asin(pSin(u,v));
    }
}

inline float Angle(const Point &p) {
    float l = Length(p);
    if (!l) return 0;
    Point norm = p / l;
    float tc = acosf(norm[0]), ts = asinf(norm[1]);
    return ts >= 0 ? tc : 2.f * M_PI - tc;
}

inline pair<float, float> RTheta(const Point &p) {
    return { Length(p), Angle(p) };
}


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


inline float Dist(const Point &p, const Point &u) {
    return Length(p-u);
}

Point Normalize(const Point &p) {
    float L = Length(p);
    if (L) return p / L;
    else return p;
}


void Circle(const Point &p, float rad = .01f) {
    Circle(p[0], p[1], rad);
}


int PutOnScreen(Point &p) {
    int ret = 0;
    for (int i = 0; i < 2; i++) {
        float &n = p[i];
        if (n < 0) { n += 1; ret++; }
        if (n > 1) { n -= 1; ret++; }
    }
    return ret;
}

int DragOnScreen(Point &p) {
    int ret = 0;
    for (int i = 0; i < 2; i++) {
        float &n = p[i];
        if (n < 0) { n = 0; ret++; }
        if (n > 1) { n = 1; ret++; }
    }
    return ret;
}

int DragOnScreen(Point &p, float buffer) {
    int ret = 0;
    for (int i = 0; i < 2; i++) {
        float &n = p[i];
        if (n < buffer) { n = buffer; ret++; }
        if (n > 1 - buffer) { n = 1 - buffer; ret++; }
    }
    return ret;
}

void BBox(Point p[], int numpts, Point &min, Point &max) {
    min = Point(INFINITY, INFINITY);
    max = Point(-INFINITY, -INFINITY);
    for (int i = 0; i < numpts; ++i) {
        Point &pt = p[i];
        for (int j = 0; j < 2; ++j) {
            min[j] = std::min(min[j], pt[j]);
            max[j] = std::max(max[j], pt[j]);
        }
    }
}

void BBox(vector<Point> p, Point &min, Point &max) {
    min = Point(INFINITY, INFINITY);
    max = Point(-INFINITY, -INFINITY);
    for (Point &pt : p) {
        for (int j = 0; j < 2; ++j) {
            min[j] = std::min(min[j], pt[j]);
            max[j] = std::max(max[j], pt[j]);
        }
    }
}

Point Centroid(vector<Point> p) {
    Point ret;
    for (Point &pt : p)
        ret += pt;
    return ret / (float)p.size();
}

float SqCenterBBox(vector<Point> &p, Point &min, Point &max, Point &centroid) {
    centroid = Centroid(p);
    float maxdist;
    for (Point &pt : p) {
        Point d(pt - centroid);
        maxdist = max3f(maxdist, abs(d[0]), abs(d[1]));
    }
    min = centroid - Point(maxdist, maxdist);
    max = centroid + Point(maxdist, maxdist);
    return maxdist;
}

float NormalizePointSet(vector<Point> &p, Point &centroid) {
    Point min, max;
    float ret = SqCenterBBox(p, min, max, centroid);
    for (Point &pt : p) {
        pt = (pt - min) / (max - min);
    }
    return ret;
}

float NormalizeFloatSet(vector<float> &f) {
    float min = INFINITY, max = -INFINITY;
    for (float n : f) {
        min = min2f(n,min);
        max = max2f(n,max);
    }
    for (float &n : f)
        n = (n - min) / (max - min);
    return .5f * (max - min);
}

float SquareDif(vector<float> &a, vector<float> &b) {
    float res = 0;
    int min = (int)std::min(a.size(), b.size());
    for (int i = 0; i < min; ++i)
        res += (a[i]-b[i])*(a[i]-b[i]);
    return res;
}

float Correlation(vector<float> &a, vector<float> &b) {
    float res = 0;
    int min = (int)std::min(a.size(), b.size());
    for (int i = 0; i < min; ++i)
        res += a[i] * b[i];
    return res;
}

float NormalPDF(float x, float u, float s) {
    return exp(-(x-u) * (x-u) / s / s * .5f) / s / sqrtf(2.f * M_PI);
}

void LineCon(float x1, float y1, float x2, float y2, float frac = 0.f) {
    float a[] = {x1,y1};
    float b[] = {x2,y2};
    Point u(a), v(b);
    Point diff = v - u;
    //float L = Length(diff);
    Point uu = u + diff * frac;
    Point vv = u + diff * (1.-frac);
    glBegin(GL_LINES);
    glVertex2f(uu[0], uu[1]);
    glVertex2f(vv[0], vv[1]);
    glEnd();
}

void DrawBox(float x, float y, float dx, float dy) {
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y);
    glVertex2f(x, y + dy);
    glVertex2f(x + dx, y + dy);
    glVertex2f(x + dx, y);
    glEnd();
}

void LineCon(Point a, Point b, float frac = 0.f) { LineCon(a[0], a[1], b[0], b[1], frac); }

void WriteStr(char *s, const Point &p) {
    WriteStr(s, p[0], p[1]);
}

void WriteStr(string s, const Point &p) {
    WriteStr(s, p[0], p[1]);
}

float ExtractFloat(string s) {
    stringstream ss;
    float ret;
    ss << s;
    ss >> ret;
    return ret;
}

void SetBrightColor(float thresh) {
    float c[] = {rng(), rng(), rng()};
    while (Luminance(c) < thresh)
        for (int i = 0; i < 3; ++i)
            c[i] = rng();
    Color(c);
}

void DrawCosts(vector<float> &costs, int n) {
    float max = costs[n*n-1] * 2.f;
    //for (float f : costs) max = max2f(max, f);
    cout << "Path is " << max / 2 << endl;

    for (int x = 0; x < n * n; ++x) {
        float col = costs[x] / max;
        Color(col,col,col);
        DrawBox((float)(x % n) / n, (float)(x / n) / n, 1.f / n, 1.f / n);
    }
}

void DrawCosts(vector<float> &costs, int m, int n) {
    float max = costs[m*n-1] * 2.f;
    
    cout << "Path is " << max / 2 << endl;
    
    for (int x = 0; x < m * n; ++x) {
        float col = costs[x] / max;
        Color(col,col,col);
        DrawBox((float)(x % m) / m, (float)(x / m) / n, 1.f / m, 1.f / n);
    }
}

float LCP(vector<float> costmat, int n, float dbias = 1.f) {
    vector<float> costs(n * n, 0);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            float cost = (!i && !j ? 0 : INFINITY);
            float multiplier = 1.f;
            if (i > 0 && cost > costs[j * n + (i-1)]) {
                cost = costs[j * n + (i-1)];
                multiplier = dbias;
            }
            if (j > 0 && cost > costs[(j-1) * n + i]) {
                cost = costs[(j-1) * n + i];
                multiplier = dbias;
            }
            if (i && j && cost > costs[(j-1) * n + (i-1)]) {
                cost = costs[(j-1) * n + (i-1)];
                multiplier = 1.f;
            }
            costs[j * n + i] = cost + costmat[j * n + i] * multiplier;
        }
    }
    return costs[n * n - 1];
}

template <class T>
float LeastCost(vector<T> &a, vector<T> &b, function<float(const T&, const T&)> costfunc,
                float dbias = 1.f, bool draw = true) {
    int n = (int)std::min(a.size(), b.size());

    vector<float> costmat(n * n, 0);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            costmat[j * n + i] = costfunc(a[i], b[j]);
    
    vector<float> costs(n * n, 0);

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            float cost = (!i && !j ? 0 : INFINITY);
            float multiplier = 1.f;
            if (i > 0 && cost > costs[j * n + (i-1)]) {
                cost = costs[j * n + (i-1)];
                multiplier = dbias;
            }
            if (j > 0 && cost > costs[(j-1) * n + i]) {
                cost = costs[(j-1) * n + i];
                multiplier = dbias;
            }
            if (i && j && cost > costs[(j-1) * n + (i-1)]) {
                cost = costs[(j-1) * n + (i-1)];
                multiplier = 1.f;
            }
            costs[j * n + i] = cost + costmat[j * n + i] * multiplier;
        }
    }

    if (draw) DrawCosts(costs, n);
    return costs[n * n - 1];
}

float LeastCostFloat(vector<float> &a, vector<float> &b, float dbias, bool draw = true) {
    return LeastCost<float>(a, b, [](const float &x, const float &y) { return (x-y)*(x-y); }, dbias, draw );
}

template <class T>
float LeastCostInside(vector<T> &a, vector<T> &b, function<float(const T&, const T&)> costfunc,
                      vector<pair<int, int>> &indices, float dbias = 1.f, bool draw = true) {
    bool verbose = false;
    dbias = 1.f;
    int m = (int)a.size();
    int n = (int)b.size();
    if (verbose) cout << "Finding " << m << " in " << n << endl;
    //cout << "n is " << n << endl << "Making cost mat" << endl;
    vector<float> costmat(n * (m + 2), 0);
    for (int j = 0; j < n; ++j)
        for (int i = 1; i < m + 1; ++i)
            costmat[j * (m+2) + i] = costfunc(a[i - 1], b[j]);
    if (verbose) cout << "Created costmat" << endl;
    vector<float> costs(n * (m + 2), 0);
    vector<int> previndex(n * (m + 2), 0);
    if (verbose) cout << "Dynprog the costs..." << endl;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m + 2; ++i) {
            float cost = (!i && !j ? 0 : INFINITY);
            previndex[j * (m+2) + i] = -1;
            if (i > 0 && cost > costs[j * (m+2) + (i-1)]) {
                cost = costs[j * (m+2) + (i-1)];
                previndex[j * (m+2) + i] = j * (m+2) + (i-1);
            }   //cost = min2f(cost, costs[j * (m+2) + (i-1)] * dbias);
            if (j > 0 && cost > costs[(j-1) * (m+2) + i]) {
                cost = costs[(j-1) * (m+2) + i];
                previndex[j * (m+2) + i] = (j-1) * (m+2) + i;
            }   //cost = min2f(cost, costs[(j-1) * (m+2) + i] * dbias);
            if (i && j && cost > costs[(j-1) * (m+2) + (i-1)]) {
                cost = costs[(j-1) * (m+2) + (i-1)];
                previndex[j * (m+2) + i] = (j-1) * (m+2) + (i-1);
            }   //cost = min2f(cost, costs[(j-1) * (m+2) + (i-1)]);
            costs[j * (m+2) + i] = cost + costmat[j * (m+2) + i];
        }
    }

    
    if (verbose) cout << "Posting result." << endl;
    if (draw)
        DrawCosts(costs, m + 2, n);
    
    
    vector<int> path;
    int index = n * (m + 2) - 1;
    Point prev(1,1);
    Color(1,0,0);
    m += 2;
    while (previndex[index] >= 0) {
        path.push_back(index);
        index = previndex[index];
        Point next((float)(index % m) / m + .5f / m, (float)(index / m) / n + .5f / n);
        if (draw) {
            LineCon(prev, next);
            prev = next;
        }
    }
    for (int i = (int)path.size() - 1; i >= 0; --i)
        indices.push_back({path[i] % m - 1, path[i] / m});
    return costs[n * m - 1];
}

float LeastCostFloatInside(vector<float> &a, vector<float> &b,
                           vector<pair<int, int>> &indices, float dbias, bool draw = true) {
    return LeastCostInside<float>(a, b, [](const float &x, const float &y) { return (x-y)*(x-y); },
                                  indices, dbias, draw );
}

float vecmax(const vector<float> &v) {
    float r = -INFINITY;
    for (float f : v)
        if (f > r)
            r = f;
    return r;
}

void DrawMatrix(vector<float> &vals, int n, float max) {
    for (int x = 0; x < n * n; ++x) {
        float col = vals[x] / max;
        Color(col,col,col);
        DrawBox((float)(x % n) / n, (float)(x / n) / n, 1.f / n, 1.f / n);
    }
}

void DrawMatrix(vector<float> &vals, int m, int n, float max) {
    for (int x = 0; x < m * n; ++x) {
        float col = vals[x] / max;
        Color(col,col,col);
        DrawBox((float)(x % m) / m, (float)(x / m) / n, 1.f / m, 1.f / n);
    }
}



ostream &operator<<(ostream &os, const vector<int> &v) {
    os << "(";
    for (int i = 0; i < v.size() - 1; ++i)
        os << v[i] << ", ";
    os << v.back() << ")";
    return os;
}
ostream &operator<<(ostream &os, const vector<float> &v) {
    os << "(";
    for (int i = 0; i < v.size() - 1; ++i)
        os << v[i] << ", ";
    os << v.back() << ")";
    return os;
}

void ReadStringIntMap(const string &filename, map<string, int> &simap) {
    ifstream infile(filename);
    cout << "Reading in string-int map " << filename << endl;
    if (infile.fail()) return;
    while (!infile.eof()) {
        string s;
        int i;
        infile >> s >> i;
        simap[s] = i;
    }
    infile.close();
}

#endif

