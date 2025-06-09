#include <windows.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <list>
#include <stack>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

// Global variables (modified)
static int selectedColor = RGB(0, 0, 0);
static int clippedColor = RGB(255, 0, 0);
static int circleColor = RGB(0, 0, 255);
static int fillColor = RGB(255, 255, 0);
static int borderColor = RGB(0, 0, 0); // Added for border color
static int selectedShape = 0;
static POINT startPt = { -1, -1 }, endPt = { -1, -1 };
static int clickCount = 0;
static bool shapeDrawn = false;
static int xc, yc, R, A, B;
static vector<POINT> splinePoints;
static enum CircleAlgorithm { DIRECT_CIRCLE = 4, POLAR, ITERATIVE_POLAR, MIDPOINT_CIRCLE, MODIFIED_MIDPOINT } currentCircleAlgorithm;
static enum EllipseAlgorithm { DIRECT_ELLIPSE = 23, POLAR_ELLIPSE, MIDPOINT_ELLIPSE } currentEllipseAlgorithm;
static enum ClippingMode {
    NONE,
    POINT_CLIP,
    LINE_CLIP,
    POLYGON_CLIP,
    DEFINE_RECT,
    DEFINE_SQUARE,
    DEFINE_CIRCLE
} clippingMode; // Removed FILL_RECURSIVE and FILL_NONRECURSIVE from ClippingMode
static vector<POINT> points;
static vector<POINT> clippingWindow;
static bool clippingWindowDefined = false;
static bool showClippedResult = false;
static vector<POINT> originalPoints;
static vector<POINT> clippedPolygonResult;
struct Square {
    int centerX, centerY;
    int size;
    int left, right, top, bottom;
    Square() : centerX(0), centerY(0), size(0) {}
    Square(int cx, int cy, int s) : centerX(cx), centerY(cy), size(s) { updateBounds(); }
    void updateBounds() {
        int halfSize = size / 2;
        left = centerX - halfSize;
        right = centerX + halfSize;
        top = centerY - halfSize;
        bottom = centerY + halfSize;
    }
};
static Square clippingSquare;
static bool clippingSquareDefined = false;
static int squareDefinitionStep = 0;
struct Circle {
    double centerX, centerY;
    double radius;
    Circle() : centerX(0), centerY(0), radius(0) {}
    Circle(double cx, double cy, double r) : centerX(cx), centerY(cy), radius(r) {}
};
static Circle clippingCircle;
static bool clippingCircleDefined = false;
static int circleDefinitionStep = 0;
static vector<POINT> fillPolygonPoints;
static bool fillingMode = false;
struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};
static vector<Point> Saved;
struct EdgeTable {
    int xleft, xright;
    EdgeTable() : xleft(10000), xright(-10000) {}
};
struct Node {
    double x;
    int ymin, ymax;
    double minv;
    Node* next;
    Node(double x, int ymin, int ymax, double minv) : x(x), ymin(ymin), ymax(ymax), minv(minv), next(nullptr) {}
};
typedef list<Node*> Table[800];
static string FileName;
static vector<int>LoadedPoints;

void LoadPointsFromFile(HWND hwnd, string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "can't open File\n";
        return;
    }

    HDC hdc = GetDC(hwnd);
    string line;
    while (getline(file, line)) {
        int x, y, r = 0, g = 0, b = 0;
        istringstream iss(line);
        if (iss >> x >> y) {
            if (iss >> r >> g >> b) {
                SetPixel(hdc, x, y, RGB(r, g, b));
                LoadedPoints.push_back(x);
                LoadedPoints.push_back(y);
                LoadedPoints.push_back(r);
                LoadedPoints.push_back(g);
                LoadedPoints.push_back(b);

            }
            else {
                LoadedPoints.push_back(x);
                LoadedPoints.push_back(y);
                LoadedPoints.push_back(r);
                LoadedPoints.push_back(g);
                LoadedPoints.push_back(b);
            }
        }
    }
    file.close();
    ReleaseDC(hwnd, hdc);
    InvalidateRect(hwnd, NULL, false); 
    MessageBox(hwnd, L"Points loaded successfully!", L"Success", MB_OK);
}

void SavePointsToFile(HWND hwnd, vector<Point>& points, string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cout << "Failed to open file";
        return;
    }

    HDC hdc = GetDC(hwnd);
    for (const auto& pt : points) {
        COLORREF color = GetPixel(hdc, pt.x, pt.y);
        int r = GetRValue(color);
        int g = GetGValue(color);
        int b = GetBValue(color);
        file << pt.x << " " << pt.y << " " << r << " " << g << " " << b << "\n";
    }
    file.close();
    ReleaseDC(hwnd, hdc);
    cout << "Points Saved Successfully\n";
}



// New function to get color from console
COLORREF GetColorFromConsole(const string& message) {
    cout << message << "\n";
    cout << "1. Black (RGB: 0, 0, 0)\n";
    cout << "2. Red (RGB: 255, 0, 0)\n";
    cout << "3. Green (RGB: 0, 255, 0)\n";
    cout << "4. Blue (RGB: 0, 0, 255)\n";
    cout << "5. Yellow (RGB: 255, 255, 0)\n";
    cout << "Enter choice (1-5): ";
    int choice;
    cin >> choice;
    switch (choice) {
    case 1: return RGB(0, 0, 0);
    case 2: return RGB(255, 0, 0);
    case 3: return RGB(0, 255, 0);
    case 4: return RGB(0, 0, 255);
    case 5: return RGB(255, 255, 0);
    default: cout << "Invalid choice, defaulting to Black\n"; return RGB(0, 0, 0);
    }
}

void FloodFillNonRecursive(HDC hdc, int x, int y, COLORREF fillColor, COLORREF targetColor) {
    if (fillColor == targetColor) return;

    stack<POINT> points;
    POINT p = { x, y };
    points.push(p);

    while (!points.empty()) {
        POINT current = points.top();
        points.pop();

        if (current.x < 0 || current.y < 0) continue;

        COLORREF currentColor = GetPixel(hdc, current.x, current.y);
        if (currentColor != targetColor || currentColor == fillColor) continue;

        SetPixel(hdc, current.x, current.y, fillColor);
        Saved.push_back(Point(x, y));


        points.push({ current.x + 1, current.y });
        points.push({ current.x - 1, current.y });
        points.push({ current.x, current.y + 1 });
        points.push({ current.x, current.y - 1 });
    }
}

void FloodFillRecursive(HDC hdc, int x, int y, COLORREF fillColor, COLORREF targetColor) {
    if (x < 0 || y < 0) return;

    COLORREF currentColor = GetPixel(hdc, x, y);
    if (currentColor != targetColor || currentColor == fillColor) return;

    SetPixel(hdc, x, y, fillColor);
    Saved.push_back(Point(x, y));

    FloodFillRecursive(hdc, x + 1, y, fillColor, targetColor);
    FloodFillRecursive(hdc, x - 1, y, fillColor, targetColor);
    FloodFillRecursive(hdc, x, y + 1, fillColor, targetColor);
    FloodFillRecursive(hdc, x, y - 1, fillColor, targetColor);
}



// Convex fill functions
void init(EdgeTable tbl[]) {
    for (int i = 0; i < 800; i++) {
        tbl[i].xleft = 10000;
        tbl[i].xright = -10000;
    }
}

void Edge2Table(Point v1, Point v2, EdgeTable tbl[]) {
    if (v1.y == v2.y) return;
    if (v1.y > v2.y) swap(v1, v2);

    int y = static_cast<int>(v1.y + 0.5);
    double x = v1.x;
    double minv = (v2.x - v1.x) / (v2.y - v1.y);

    while (y < v2.y) {
        if (x < tbl[y].xleft) tbl[y].xleft = static_cast<int>(ceil(x));
        if (x > tbl[y].xright) tbl[y].xright = static_cast<int>(floor(x));
        y++;
        x += minv;
    }
}

void Polygon2Table(Point p[], int n, EdgeTable tbl[]) {
    Point v1 = p[n - 1];
    for (int i = 0; i < n; i++) {
        Point v2 = p[i];
        Edge2Table(v1, v2, tbl);
        v1 = p[i];
    }
}

void Table2Screen(HDC hdc, EdgeTable tbl[], COLORREF c) {
    for (int i = 0; i < 800; i++) {
        if (tbl[i].xleft < tbl[i].xright) {
            for (int x = tbl[i].xleft; x <= tbl[i].xright; x++) {
                SetPixel(hdc, x, i, c);
                Saved.push_back(Point(x, i));

            }
        }
    }
}

void ConvexFill(HDC hdc, Point p[], int n, COLORREF c) {
    EdgeTable tbl[800];
    init(tbl);
    Polygon2Table(p, n, tbl);
    Table2Screen(hdc, tbl, c);
}

// Non-convex fill functions
void initET(Table& et) {
    for (int i = 0; i < 800; i++) {
        et[i].clear();
    }
}

void EdgeToET(Point v1, Point v2, Table& et) {
    if (abs(v1.y - v2.y) <= 0.5) return;
    if (v1.y > v2.y) swap(v1, v2);

    double x = v1.x;
    int ymin = static_cast<int>(v1.y + 0.5);
    int ymax = static_cast<int>(v2.y + 0.5);
    double minv = (v2.x - v1.x) / (v2.y - v1.y);

    if (abs(minv) > 100.0 || isnan(minv) || isinf(minv)) return;

    Node* node = new Node(x, ymin, ymax, minv);
    et[ymin].push_back(node);
}

void PolygonToET(Point p[], int n, Table& et) {
    Point v1 = p[n - 1];
    for (int i = 0; i < n; i++) {
        Point v2 = p[i];
        EdgeToET(v1, v2, et);
        v1 = p[i];
    }
}

void GeneralPolygonFill(HDC hdc, Point p[], int n, COLORREF c) {
    Table et;
    initET(et);
    PolygonToET(p, n, et);

    int ymin = 800, ymax = 0;
    for (int i = 0; i < n; i++) {
        int y = static_cast<int>(p[i].y + 0.5);
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }

    list<Node*> aet; // Active Edge Table

    for (int y = ymin; y <= ymax; y++) {
        for (Node* node : et[y]) {
            aet.push_back(node);
        }

        aet.sort([](Node* a, Node* b) { return a->x < b->x; });

        auto it = aet.begin();
        while (it != aet.end()) {
            auto it2 = it;
            it2++;
            if (it2 == aet.end()) break;

            Node* edge1 = *it;
            Node* edge2 = *it2;

            int x1 = static_cast<int>(edge1->x + 0.5);
            int x2 = static_cast<int>(edge2->x + 0.5);

            if (x1 != x2) {
                for (int x = x1; x <= x2; x++) {
                    SetPixel(hdc, x, y, c);
                    Saved.push_back(Point(x, y));
                }
            }

            it++;
            it++;
        }

        aet.remove_if([y](Node* node) { return node->ymax == y; });

        for (Node* node : aet) {
            node->x += node->minv;
        }
    }

    // Clean up memory
    for (int i = 0; i < 800; i++) {
        for (Node* node : et[i]) {
            delete node;
        }
    }
}

int Round(double x) { return (int)(x + 0.5); }

void Draw8Points(HDC hdc, int xc, int yc, int x, int y, COLORREF color) {
    SetPixel(hdc, xc + x, yc + y, color);
    SetPixel(hdc, xc - x, yc + y, color);
    SetPixel(hdc, xc + x, yc - y, color);
    SetPixel(hdc, xc - x, yc - y, color);
    SetPixel(hdc, xc + y, yc + x, color);
    SetPixel(hdc, xc - y, yc + x, color);
    SetPixel(hdc, xc + y, yc - x, color);
    SetPixel(hdc, xc - y, yc - x, color);
    Saved.push_back(Point(xc + x, yc + y));
    Saved.push_back(Point(xc - x, yc + y));
    Saved.push_back(Point(xc + x, yc - y));
    Saved.push_back(Point(xc - x, yc - y));
    Saved.push_back(Point(xc + y, yc + x));
    Saved.push_back(Point(xc - y, yc + x));
    Saved.push_back(Point(xc + y, yc - x));
    Saved.push_back(Point(xc - y, yc - x));
}

void DrawCircleDirect(HDC hdc, int xc, int yc, int R, COLORREF color) {
    for (int x = 0; x <= R / sqrt(2); ++x) {
        int y = round(sqrt(R * R - x * x));
        Draw8Points(hdc, xc, yc, x, y, color);
    }
}

void DrawCirclePolar(HDC hdc, int xc, int yc, int R, COLORREF color) {
    double dTheta = 1.0 / R;
    for (double theta = 0; theta <= 3.14159 / 4; theta += dTheta) {
        int x = round(R * cos(theta));
        int y = round(R * sin(theta));
        Draw8Points(hdc, xc, yc, x, y, color);
    }
}

void DrawCircleIterativePolar(HDC hdc, int xc, int yc, int R, COLORREF color) {
    double x = R, y = 0;
    double dTheta = 1.0 / R;
    double cos_d = cos(dTheta), sin_d = sin(dTheta);
    for (int i = 0; i <= R * 3.14159 / 4; ++i) {
        Draw8Points(hdc, xc, yc, round(x), round(y), color);
        double xNew = x * cos_d - y * sin_d;
        y = x * sin_d + y * cos_d;
        x = xNew;
    }
}

void DrawCircleMidpoint(HDC hdc, int xc, int yc, int R, COLORREF color) {
    int x = 0, y = R;
    int d = 1 - R;
    Draw8Points(hdc, xc, yc, x, y, color);
    while (x < y) {
        if (d < 0)
            d += 2 * x + 3;
        else {
            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
        Draw8Points(hdc, xc, yc, x, y, color);
    }
}

void DrawCircleModifiedMidpoint(HDC hdc, int xc, int yc, int R, COLORREF color) {
    int x = 0, y = R;
    int d = 1 - R;
    int dE = 3, dSE = 5 - 2 * R;
    Draw8Points(hdc, xc, yc, x, y, color);
    while (x < y) {
        if (d < 0) {
            d += dE;
            dE += 2;
            dSE += 2;
        }
        else {
            d += dSE;
            dE += 2;
            dSE += 4;
            y--;
        }
        x++;
        Draw8Points(hdc, xc, yc, x, y, color);
    }
}

// DDA Line Drawing Function
void drawDDA(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));
    if (steps == 0) {
        SetPixel(hdc, x1, y1, c);
        Saved.push_back(Point(x1, y1));
        return;
    }
    float xInc = (float)dx / steps;
    float yInc = (float)dy / steps;
    float x = x1, y = y1;
    for (int i = 0; i <= steps; i++) {
        SetPixel(hdc, round(x), round(y), c);
        Saved.push_back(Point(x, y));

        x += xInc;
        y += yInc;
    }
}

// Midpoint Line Drawing Function
void DrawLineMidpoint(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int sx = (dx >= 0) ? 1 : -1;
    int sy = (dy >= 0) ? 1 : -1;
    dx = abs(dx);
    dy = abs(dy);

    int x = x1;
    int y = y1;
    SetPixel(hdc, x, y, c);
    Saved.push_back(Point(x, y));

    if (dx > dy) {
        int d = (2 * dy) - dx;
        int incrE = 2 * dy;
        int incrNE = 2 * (dy - dx);

        while (x != x2) {
            x += sx;
            if (d <= 0) {
                d += incrE;
            }
            else {
                y += sy;
                d += incrNE;
            }
            SetPixel(hdc, x, y, c);
            Saved.push_back(Point(x, y));
        }
    }
    else {
        int d = (2 * dx) - dy;
        int incrN = 2 * dx;
        int incrNE = 2 * (dx - dy);

        while (y != y2) {
            y += sy;
            if (d <= 0) {
                d += incrN;
            }
            else {
                x += sx;
                d += incrNE;
            }
            SetPixel(hdc, x, y, c);
            Saved.push_back(Point(x, y));

        }
    }
}

// Parametric Line Drawing Function
void DrawLineParametric(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    for (float t = 0; t <= 1; t += 0.001f) {
        int x = round(x1 + t * (x2 - x1));
        int y = round(y1 + t * (y2 - y1));
        SetPixel(hdc, x, y, color);
        Saved.push_back(Point(x, y));

    }
}

// Hermite Curve Drawing Function
void DrawHermiteCurve(HDC hdc, int x1, int y1, int u1, int v1, int x2, int y2, int u2, int v2, COLORREF color, int numPoints = 100000) {
    const double H[4][4] = {
            { 2,  1, -2,  1},
            {-3, -2,  3, -1},
            { 0,  1,  0,  0},
            { 1,  0,  0,  0}
    };

    double gx[4] = { x1, u1, x2, u2 };
    double gy[4] = { y1, v1, y2, v2 };

    double cx[4] = { 0 };
    double cy[4] = { 0 };

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cx[i] += H[i][j] * gx[j];
            cy[i] += H[i][j] * gy[j];
        }
    }

    double step = 1.0 / numPoints;
    for (double t = 0; t <= 1.0; t += step) {
        double t2 = t * t;
        double t3 = t2 * t;
        int x = static_cast<int>(cx[0] * t3 + cx[1] * t2 + cx[2] * t + cx[3]);
        int y = static_cast<int>(cy[0] * t3 + cy[1] * t2 + cy[2] * t + cy[3]);
        SetPixel(hdc, x, y, color);
        Saved.push_back(Point(x, y));

    }
}

// Bezier Curve Drawing Function
int binomialCoefficient(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    k = min(k, n - k);
    int result = 1;
    for (int i = 1; i <= k; i++) {
        result *= (n - k + i);
        result /= i;
    }
    return result;
}

Point calculateBezierPoint(vector<Point>& controlPoints, double t) {
    int n = controlPoints.size() - 1;
    Point result(0, 0);
    for (int i = 0; i <= n; i++) {
        double blend = binomialCoefficient(n, i) * pow(t, i) * pow(1 - t, n - i);
        result.x += controlPoints[i].x * blend;
        result.y += controlPoints[i].y * blend;
    }
    return result;
}

void DrawBezierCurve(HDC hdc, vector<Point>& controlPoints, COLORREF color, int numPoints = 1000) {
    if (controlPoints.size() < 2) return;
    for (const auto& point : controlPoints) {
        int x = static_cast<int>(point.x);
        int y = static_cast<int>(point.y);
        SetPixel(hdc, x, y, RGB(255, 0, 0));
        Saved.push_back(Point(x, y));

    }
    double step = 1.0 / numPoints;
    for (double t = 0; t <= 1.0; t += step) {
        Point p = calculateBezierPoint(controlPoints, t);
        SetPixel(hdc, static_cast<int>(p.x), static_cast<int>(p.y), color);
        Saved.push_back(Point(p.x, p.y));

    }
}

// Optimized Hermite Curve Drawing Function
void DrawHermiteCurveOptimized(HDC hdc, double x1, double y1, double u1, double v1, double x2, double y2, double u2, double v2, COLORREF color, int numPoints = 1000) {
    double ax = 2 * x1 + u1 - 2 * x2 + u2;
    double bx = -3 * x1 - 2 * u1 + 3 * x2 - u2;
    double cx = u1;
    double dx = x1;
    double ay = 2 * y1 + v1 - 2 * y2 + v2;
    double by = -3 * y1 - 2 * v1 + 3 * y2 - v2;
    double cy = v1;
    double dy = y1;
    double step = 1.0 / numPoints;
    for (double t = 0; t <= 1.0; t += step) {
        double t2 = t * t;
        double t3 = t2 * t;
        int x = static_cast<int>(ax * t3 + bx * t2 + cx * t + dx + 0.5);
        int y = static_cast<int>(ay * t3 + by * t2 + cy * t + dy + 0.5);
        SetPixel(hdc, x, y, color);
        Saved.push_back(Point(x, y));

    }
}

// Cardinal Spline Drawing Function
void DrawCardinalSpline(HDC hdc, vector<POINT>& points, COLORREF color, double tension = 0.0) {
    if (points.size() < 4) return;
    double s = (1 - tension) / 2.0;
    for (size_t i = 1; i + 2 < points.size(); ++i) {
        POINT p0 = points[i - 1];
        POINT p1 = points[i];
        POINT p2 = points[i + 1];
        POINT p3 = points[i + 2];
        double u1 = s * (p2.x - p0.x);
        double v1 = s * (p2.y - p0.y);
        double u2 = s * (p3.x - p1.x);
        double v2 = s * (p3.y - p1.y);
        DrawHermiteCurveOptimized(hdc, p1.x, p1.y, u1, v1, p2.x, p2.y, u2, v2, color);
    }
}

// Fill Square with Vertical Hermite Curves
void FillSquareWithHermiteVertical(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    int side = min(abs(x2 - x1), abs(y2 - y1));
    int left = x1;
    int top = y1;
    int right = left + side;
    int bottom = top + side;
    if (x2 < x1) {
        left = x1 - side;
        right = x1;
    }
    if (y2 < y1) {
        top = y1 - side;
        bottom = y1;
    }
    for (int x = left; x <= right; x++) {
        SetPixel(hdc, x, top, color);
        SetPixel(hdc, x, bottom, color);
        Saved.push_back(Point(x, top));
        Saved.push_back(Point(left, bottom));
    }
    for (int y = top; y <= bottom; y++) {
        SetPixel(hdc, left, y, color);
        SetPixel(hdc, right, y, color);
        Saved.push_back(Point(left, y));
        Saved.push_back(Point(right, y));

    }
    const int baseSpacing = max(5, side / 20);
    const int pairSpacing = baseSpacing * 2;
    const int tangentLength = min(side / 4, 50);
    int pairCount = 0;
    for (int x = left + baseSpacing; x < right; x += baseSpacing) {
        const int px1 = x;
        const int py1 = top;
        const int px2 = x;
        const int py2 = bottom;
        bool alternateDirection = (pairCount % 2) == 0;
        const int tx1 = alternateDirection ? tangentLength : -tangentLength;
        const int ty1 = 0;
        const int tx2 = alternateDirection ? -tangentLength : tangentLength;
        const int ty2 = 0;
        const double ax = 2 * px1 + tx1 - 2 * px2 + tx2;
        const double bx = -3 * px1 - 2 * tx1 + 3 * px2 - tx2;
        const double cx = tx1;
        const double dx = px1;
        const double ay = 2 * py1 + ty1 - 2 * py2 + ty2;
        const double by = -3 * py1 - 2 * ty1 + 3 * py2 - ty2;
        const double cy = ty1;
        const double dy = py1;
        const double step = 0.001;
        for (double t = 0; t <= 1.0; t += step) {
            double t2 = t * t;
            double t3 = t2 * t;
            int curveX = static_cast<int>(ax * t3 + bx * t2 + cx * t + dx);
            int curveY = static_cast<int>(ay * t3 + by * t2 + cy * t + dy);
            if (curveX >= left && curveX <= right && curveY >= top && curveY <= bottom) {
                SetPixel(hdc, curveX, curveY, color);
                Saved.push_back(Point(curveX, curveY));
            }
        }
        pairCount++;
        if (pairCount % 2 == 0) {
            x += baseSpacing;
        }
    }
}

// Fill Rectangle with Horizontal Bezier Curves
void FillRectangleWithBezierHorizontal(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    const int left = min(x1, x2);
    const int right = max(x1, x2);
    const int top = min(y1, y2);
    const int bottom = max(y1, y2);
    const int width = right - left;
    const int height = bottom - top;
    drawDDA(hdc, left, top, right, top, color);
    drawDDA(hdc, right, top, right, bottom, color);
    drawDDA(hdc, right, bottom, left, bottom, color);
    drawDDA(hdc, left, bottom, left, top, color);
    const int spacing = max(5, height / 15);
    const int waveAmplitude = max(3, min(height / 10, (bottom - top - spacing) / 2));
    HRGN clipRegion = CreateRectRgn(left, top, right, bottom);
    SelectClipRgn(hdc, clipRegion);
    for (int y = top + spacing; y < bottom; y += spacing) {
        int boundedY = min(y, bottom - 1);
        vector<Point> controlPoints;
        controlPoints.push_back(Point(left, boundedY));
        int controlY1 = boundedY + waveAmplitude;
        controlY1 = max(top + 1, min(bottom - 1, controlY1));
        controlPoints.push_back(Point(left + width / 3, controlY1));
        int controlY2 = boundedY - waveAmplitude;
        controlY2 = max(top + 1, min(bottom - 1, controlY2));
        controlPoints.push_back(Point(left + 2 * width / 3, controlY2));
        controlPoints.push_back(Point(right, boundedY));
        DrawBezierCurve(hdc, controlPoints, color, 1000);
    }
    SelectClipRgn(hdc, NULL);
    DeleteObject(clipRegion);
}

// Detect Quarter for Circle Filling
int DetectQuarter(int x, int y) {
    if (x > xc && y < yc) return 1;
    if (x < xc && y < yc) return 2;
    if (x < xc && y > yc) return 3;
    if (x > xc && y > yc) return 4;
    return 0;
}

// Fill Circle Quarter with Spaced Lines
void FillQuarterWithSpacedLines(HDC hdc, int quarter, COLORREF color, int spacing = 5) {
    int startAngle, endAngle;
    switch (quarter) {
    case 1: startAngle = 0; endAngle = 90; break;
    case 2: startAngle = 90; endAngle = 180; break;
    case 3: startAngle = 180; endAngle = 270; break;
    case 4: startAngle = 270; endAngle = 360; break;
    default: return;
    }
    for (int angle = startAngle; angle < endAngle; angle += spacing) {
        double radian = angle * 3.141592653589793 / 180.0;
        int x2 = xc + R * cos(radian);
        int y2 = yc - R * sin(radian);
        drawDDA(hdc, xc, yc, x2, y2, color);
    }
    HPEN hPen = CreatePen(PS_SOLID, 1, color);
    HPEN hOldPen = (HPEN)SelectObject(hdc, hPen);
    switch (quarter) {
    case 1:
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc + R, yc);
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc, yc - R);
        break;
    case 2:
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc - R, yc);
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc, yc - R);
        break;
    case 3:
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc - R, yc);
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc, yc + R);
        break;
    case 4:
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc + R, yc);
        MoveToEx(hdc, xc, yc, NULL);
        LineTo(hdc, xc, yc + R);
        break;
    }
    SelectObject(hdc, hOldPen);
    DeleteObject(hPen);
}

// Fill Circle Quarter with Smaller Circles
void FillQuarterWithCircles(HDC hdc, int quarter, COLORREF color) {
    int smallRadius = R / 20;
    if (smallRadius < 3) smallRadius = 3;
    int spacing = smallRadius * 2 + 2;
    int startX = xc - R, endX = xc + R;
    int startY = yc - R, endY = yc + R;
    for (int y = startY; y <= endY; y += spacing) {
        for (int x = startX; x <= endX; x += spacing) {
            double dx = x - xc;
            double dy = y - yc;
            double distance = sqrt(dx * dx + dy * dy);
            if (distance + smallRadius > R) continue;
            double angle = atan2(-dy, dx) * 180.0 / 3.14159265;
            if (angle < 0) angle += 360;
            bool inQuarter = false;
            switch (quarter) {
            case 1: inQuarter = (angle >= 0 && angle < 90); break;
            case 2: inQuarter = (angle >= 90 && angle < 180); break;
            case 3: inQuarter = (angle >= 180 && angle < 270); break;
            case 4: inQuarter = (angle >= 270 && angle < 360); break;
            }
            if (inQuarter) {
                DrawCircleModifiedMidpoint(hdc, x, y, smallRadius, color);
            }
        }
    }
}

// Ellipse Drawing Functions
void Draw4Points(HDC hdc, int xc, int yc, int x, int y, COLORREF c) {
    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc - x, yc - y, c);
    Saved.push_back(Point(xc + x, yc + y));
    Saved.push_back(Point(xc - x, yc + y));
    Saved.push_back(Point(xc + x, yc - y));
    Saved.push_back(Point(xc - x, yc - y));

}

void DirectEllipse(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {
    for (int x = 0; x <= A; ++x) {
        int y = round(B * sqrt(max(0.0, 1.0 - (double)x * x / (A * A))));
        Draw4Points(hdc, xc, yc, x, y, c);
    }
    for (int y = 0; y <= B; ++y) {
        int x = round(A * sqrt(max(0.0, 1.0 - (double)y * y / (B * B))));
        Draw4Points(hdc, xc, yc, x, y, c);
    }
}

void PolarEllipse(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {
    double deltaTheta = 1.0 / max(A, B);
    for (double theta = 0; theta <= 3.14159 / 2; theta += deltaTheta) {
        int x = round(A * cos(theta));
        int y = round(B * sin(theta));
        Draw4Points(hdc, xc, yc, x, y, c);
    }
}

void MidPointEllipse(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {
    double x = 0, y = B;
    double A2 = A * A, B2 = B * B;
    double d1 = B2 - A2 * B + 0.25 * A2;
    Draw4Points(hdc, xc, yc, x, y, c);
    while (B2 * (x + 1) < A2 * y) {
        x++;
        if (d1 < 0) {
            d1 += B2 * (2 * x + 1);
        }
        else {
            y--;
            d1 += B2 * (2 * x + 1) - 2 * A2 * y;
        }
        Draw4Points(hdc, xc, yc, round(x), round(y), c);
    }
    double d2 = B2 * (x + 0.5) * (x + 0.5) + A2 * (y - 1) * (y - 1) - A2 * B2;
    while (y > 0) {
        y--;
        if (d2 > 0) {
            d2 -= A2 * (2 * y + 1);
        }
        else {
            x++;
            d2 += B2 * (2 * x) - A2 * (2 * y + 1);
        }
        Draw4Points(hdc, xc, yc, round(x), round(y), c);
    }
}

// Clipping Functions
void drawPoint(HDC hdc, int x, int y, COLORREF color) {
    for (int dx = -3; dx <= 3; dx++) {
        for (int dy = -3; dy <= 3; dy++) {
            if (dx * dx + dy * dy <= 9) {
                SetPixel(hdc, x + dx, y + dy, color);
                Saved.push_back(Point(x, y));
            }
        }
    }
}

int computeCode(int x, int y, int xmin, int ymin, int xmax, int ymax) {
    int code = 0;
    if (x < xmin) code |= 1;
    else if (x > xmax) code |= 2;
    if (y < ymin) code |= 4;
    else if (y > ymax) code |= 8;
    return code;
}

bool isPointInside(int x, int y, int xmin, int ymin, int xmax, int ymax) {
    return (x >= xmin && x <= xmax && y >= ymin && y <= ymax);
}

bool clipLine(int x1, int y1, int x2, int y2, int& cx1, int& cy1, int& cx2, int& cy2, int xmin, int ymin, int xmax, int ymax) {
    cx1 = x1; cy1 = y1; cx2 = x2; cy2 = y2;
    int code1 = computeCode(cx1, cy1, xmin, ymin, xmax, ymax);
    int code2 = computeCode(cx2, cy2, xmin, ymin, xmax, ymax);
    bool accept = false;
    while (true) {
        if (!(code1 | code2)) {
            accept = true;
            break;
        }
        else if (code1 & code2) {
            break;
        }
        else {
            int codeOut = code1 ? code1 : code2;
            double x, y;
            if (codeOut & 4) {
                x = cx1 + (double)(cx2 - cx1) * (ymin - cy1) / (cy2 - cy1);
                y = ymin;
            }
            else if (codeOut & 8) {
                x = cx1 + (double)(cx2 - cx1) * (ymax - cy1) / (cy2 - cy1);
                y = ymax;
            }
            else if (codeOut & 2) {
                y = cy1 + (double)(cy2 - cy1) * (xmax - cx1) / (cx2 - cx1);
                x = xmax;
            }
            else if (codeOut & 1) {
                y = cy1 + (double)(cy2 - cy1) * (xmin - cx1) / (cx2 - cx1);
                x = xmin;
            }
            if (codeOut == code1) {
                cx1 = Round(x);
                cy1 = Round(y);
                code1 = computeCode(cx1, cy1, xmin, ymin, xmax, ymax);
            }
            else {
                cx2 = Round(x);
                cy2 = Round(y);
                code2 = computeCode(cx2, cy2, xmin, ymin, xmax, ymax);
            }
        }
    }
    return accept;
}

vector<POINT> clipPolygon(const vector<POINT>& polygon, int xmin, int ymin, int xmax, int ymax) {
    vector<POINT> output;
    if (polygon.size() < 3) return output;
    vector<POINT> current = polygon;
    output.clear();
    if (!current.empty()) {
        for (size_t i = 0; i < current.size(); i++) {
            POINT p1 = current[i];
            POINT p2 = current[(i + 1) % current.size()];
            if (p1.x >= xmin) {
                if (p2.x >= xmin) {
                    output.push_back(p2);
                }
                else {
                    if (p2.x != p1.x) {
                        POINT intersect;
                        intersect.x = xmin;
                        intersect.y = p1.y + (xmin - p1.x) * (p2.y - p1.y) / (p2.x - p1.x);
                        output.push_back(intersect);
                    }
                }
            }
            else {
                if (p2.x >= xmin) {
                    if (p2.x != p1.x) {
                        POINT intersect;
                        intersect.x = xmin;
                        intersect.y = p1.y + (xmin - p1.x) * (p2.y - p1.y) / (p2.x - p1.x);
                        output.push_back(intersect);
                    }
                    output.push_back(p2);
                }
            }
        }
        current = output;
    }
    output.clear();
    if (!current.empty()) {
        for (size_t i = 0; i < current.size(); i++) {
            POINT p1 = current[i];
            POINT p2 = current[(i + 1) % current.size()];
            if (p1.x <= xmax) {
                if (p2.x <= xmax) {
                    output.push_back(p2);
                }
                else {
                    if (p2.x != p1.x) {
                        POINT intersect;
                        intersect.x = xmax;
                        intersect.y = p1.y + (xmax - p1.x) * (p2.y - p1.y) / (p2.x - p1.x);
                        output.push_back(intersect);
                    }
                }
            }
            else {
                if (p2.x <= xmax) {
                    if (p2.x != p1.x) {
                        POINT intersect;
                        intersect.x = xmax;
                        intersect.y = p1.y + (xmax - p1.x) * (p2.y - p1.y) / (p2.x - p1.x);
                        output.push_back(intersect);
                    }
                    output.push_back(p2);
                }
            }
        }
        current = output;
    }
    output.clear();
    if (!current.empty()) {
        for (size_t i = 0; i < current.size(); i++) {
            POINT p1 = current[i];
            POINT p2 = current[(i + 1) % current.size()];
            if (p1.y >= ymin) {
                if (p2.y >= ymin) {
                    output.push_back(p2);
                }
                else {
                    if (p2.y != p1.y) {
                        POINT intersect;
                        intersect.y = ymin;
                        intersect.x = p1.x + (ymin - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
                        output.push_back(intersect);
                    }
                }
            }
            else {
                if (p2.y >= ymin) {
                    if (p2.y != p1.y) {
                        POINT intersect;
                        intersect.y = ymin;
                        intersect.x = p1.x + (ymin - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
                        output.push_back(intersect);
                    }
                    output.push_back(p2);
                }
            }
        }
        current = output;
    }
    output.clear();
    if (!current.empty()) {
        for (size_t i = 0; i < current.size(); i++) {
            POINT p1 = current[i];
            POINT p2 = current[(i + 1) % current.size()];
            if (p1.y <= ymax) {
                if (p2.y <= ymax) {
                    output.push_back(p2);
                }
                else {
                    if (p2.y != p1.y) {
                        POINT intersect;
                        intersect.y = ymax;
                        intersect.x = p1.x + (ymax - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
                        output.push_back(intersect);
                    }
                }
            }
            else {
                if (p2.y <= ymax) {
                    if (p2.y != p1.y) {
                        POINT intersect;
                        intersect.y = ymax;
                        intersect.x = p1.x + (ymax - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
                        output.push_back(intersect);
                    }
                    output.push_back(p2);
                }
            }
        }
    }
    return output;
}

bool clipLineToSquare(int x1, int y1, int x2, int y2, int& cx1, int& cy1, int& cx2, int& cy2, const Square& square) {
    cx1 = x1; cy1 = y1; cx2 = x2; cy2 = y2;
    int code1 = computeCode(cx1, cy1, square.left, square.top, square.right, square.bottom);
    int code2 = computeCode(cx2, cy2, square.left, square.top, square.right, square.bottom);
    bool accept = false;
    while (true) {
        if (!(code1 | code2)) {
            accept = true;
            break;
        }
        else if (code1 & code2) {
            break;
        }
        else {
            int codeOut = code1 ? code1 : code2;
            double x, y;
            if (codeOut & 4) {
                x = cx1 + (double)(cx2 - cx1) * (square.top - cy1) / (cy2 - cy1);
                y = square.top;
            }
            else if (codeOut & 8) {
                x = cx1 + (double)(cx2 - cx1) * (square.bottom - cy1) / (cy2 - cy1);
                y = square.bottom;
            }
            else if (codeOut & 2) {
                y = cy1 + (double)(cy2 - cy1) * (square.right - cx1) / (cx2 - cx1);
                x = square.right;
            }
            else if (codeOut & 1) {
                y = cy1 + (double)(cy2 - cy1) * (square.left - cx1) / (cx2 - cx1);
                x = square.left;
            }
            if (codeOut == code1) {
                cx1 = Round(x);
                cy1 = Round(y);
                code1 = computeCode(cx1, cy1, square.left, square.top, square.right, square.bottom);
            }
            else {
                cx2 = Round(x);
                cy2 = Round(y);
                code2 = computeCode(cx2, cy2, square.left, square.top, square.right, square.bottom);
            }
        }
    }
    return accept;
}

// Circle Clipping Functions
bool isPointInsideCircle(int x, int y, const Circle& circle) {
    double dist = sqrt(pow(x - circle.centerX, 2) + pow(y - circle.centerY, 2));
    return dist <= circle.radius;
}

bool clipLineToCircle(double x1, double y1, double x2, double y2,
    double& cx1, double& cy1, double& cx2, double& cy2,
    const Circle& circle) {

    // Translate line to circle-centered coordinates
    double dx = x2 - x1;
    double dy = y2 - y1;
    double fx = x1 - circle.centerX;
    double fy = y1 - circle.centerY;

    // Quadratic equation coefficients: at² + bt + c = 0
    double a = dx * dx + dy * dy;
    double b = 2 * (fx * dx + fy * dy);
    double c = (fx * fx + fy * fy) - circle.radius * circle.radius;

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        // No intersection
        return false;
    }

    if (discriminant == 0) {
        // One intersection point (tangent)
        double t = -b / (2 * a);
        if (t >= 0 && t <= 1) {
            cx1 = cx2 = x1 + t * dx;
            cy1 = cy2 = y1 + t * dy;
            return true;
        }
        return false;
    }

    // Two intersection points
    double sqrtD = sqrt(discriminant);
    double t1 = (-b - sqrtD) / (2 * a);
    double t2 = (-b + sqrtD) / (2 * a);

    // Ensure t1 <= t2
    if (t1 > t2) {
        double temp = t1;
        t1 = t2;
        t2 = temp;
    }

    // Check if line segment intersects circle
    if (t2 < 0 || t1 > 1) {
        // Line segment is completely outside circle
        return false;
    }

    // Clamp intersection parameters to [0, 1]
    double tStart = max(0.0, t1);
    double tEnd = min(1.0, t2);

    // Calculate intersection points
    cx1 = x1 + tStart * dx;
    cy1 = y1 + tStart * dy;
    cx2 = x1 + tEnd * dx;
    cy2 = y1 + tEnd * dy;

    // Check if we have a valid line segment
    if (tStart >= tEnd) {
        return false;
    }

    return true;
}

void GetFileName(string Message)
{
    cout << Message << endl;
    cin >> FileName;
}

// Clear Window Function
void ClearWindow(HWND hwnd) {
    InvalidateRect(hwnd, NULL, TRUE);
    UpdateWindow(hwnd);
    shapeDrawn = false;
    clickCount = 0;
    splinePoints.clear();
    points.clear();
    originalPoints.clear();
    clippingWindow.clear();
    clippedPolygonResult.clear();
    clippingWindowDefined = false;
    clippingSquareDefined = false;
    clippingSquare = Square();
    squareDefinitionStep = 0;
    clippingCircleDefined = false;
    clippingCircle = Circle();
    circleDefinitionStep = 0;
    clippingMode = NONE;
    showClippedResult = false;
    fillPolygonPoints.clear(); // Clear filling points
    fillingMode = false;
    Saved.clear();
}


// Modified WndProc
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    HDC hdc;
    HMENU hMenu;
    POINT pt;

    switch (msg) {
    case WM_CREATE: {
        hMenu = CreateMenu();
        HMENU hDrawMenu = CreatePopupMenu();
        HMENU hColorMenu = CreatePopupMenu();
        HMENU hOptionsMenu = CreatePopupMenu();
        HMENU hRectClipMenu = CreatePopupMenu();
        HMENU hSquareClipMenu = CreatePopupMenu();
        HMENU hCircleClipMenu = CreatePopupMenu();

        // Draw Menu (added flood fill options)
        AppendMenu(hDrawMenu, MF_STRING, 1, L"Line (DDA)");
        AppendMenu(hDrawMenu, MF_STRING, 2, L"Line (Midpoint)");
        AppendMenu(hDrawMenu, MF_STRING, 3, L"Line (Parametric)");
        AppendMenu(hDrawMenu, MF_STRING, 4, L"Circle (Direct)");
        AppendMenu(hDrawMenu, MF_STRING, 5, L"Circle (Polar)");
        AppendMenu(hDrawMenu, MF_STRING, 6, L"Circle (Iterative Polar)");
        AppendMenu(hDrawMenu, MF_STRING, 7, L"Circle (Midpoint)");
        AppendMenu(hDrawMenu, MF_STRING, 8, L"Circle (Modified Midpoint)");
        AppendMenu(hDrawMenu, MF_STRING, 9, L"Filling Circle with Lines");
        AppendMenu(hDrawMenu, MF_STRING, 10, L"Filling Circle with Circles");
        AppendMenu(hDrawMenu, MF_STRING, 11, L"Filling Square with Hermite Curve (Vertical)");
        AppendMenu(hDrawMenu, MF_STRING, 12, L"Filling Rectangle with Bezier Curve (Horizontal)");
        AppendMenu(hDrawMenu, MF_STRING, 13, L"Cardinal Spline Curve");
        AppendMenu(hDrawMenu, MF_STRING, 23, L"Ellipse (Direct)");
        AppendMenu(hDrawMenu, MF_STRING, 24, L"Ellipse (Polar)");
        AppendMenu(hDrawMenu, MF_STRING, 25, L"Ellipse (Midpoint)");
        AppendMenu(hDrawMenu, MF_STRING, 30, L"Convex Fill");
        AppendMenu(hDrawMenu, MF_STRING, 34, L"Non-Convex Fill");
        AppendMenu(hDrawMenu, MF_STRING, 38, L"Recursive Flood Fill");
        AppendMenu(hDrawMenu, MF_STRING, 39, L"Non-Recursive Flood Fill");

        // Color Menu
        AppendMenu(hColorMenu, MF_STRING, 14, L"Black");
        AppendMenu(hColorMenu, MF_STRING, 15, L"Red");
        AppendMenu(hColorMenu, MF_STRING, 16, L"Green");
        AppendMenu(hColorMenu, MF_STRING, 17, L"Blue");
        AppendMenu(hColorMenu, MF_STRING, 18, L"Yellow");

        // Rectangle Clipping Menu
        AppendMenu(hRectClipMenu, MF_STRING, 26, L"Point Clipping");
        AppendMenu(hRectClipMenu, MF_STRING, 27, L"Line Clipping");
        AppendMenu(hRectClipMenu, MF_STRING, 28, L"Polygon Clipping");
        AppendMenu(hRectClipMenu, MF_STRING, 29, L"Define Rectangle Clipping Window");

        // Square Clipping Menu
        AppendMenu(hSquareClipMenu, MF_STRING, 31, L"Point Clipping");
        AppendMenu(hSquareClipMenu, MF_STRING, 32, L"Line Clipping");
        AppendMenu(hSquareClipMenu, MF_STRING, 33, L"Define Square Clipping Window");

        // Circle Clipping Menu
        AppendMenu(hCircleClipMenu, MF_STRING, 35, L"Define Circle Clipping Window");
        AppendMenu(hCircleClipMenu, MF_SEPARATOR, 0, NULL);
        AppendMenu(hCircleClipMenu, MF_STRING, 36, L"Point Clipping");
        AppendMenu(hCircleClipMenu, MF_STRING, 37, L"Line Clipping");

        // Options Menu
        AppendMenu(hOptionsMenu, MF_STRING, 19, L"Change Shape of Window");
        AppendMenu(hOptionsMenu, MF_STRING, 20, L"Clear Screen");
        AppendMenu(hOptionsMenu, MF_STRING, 21, L"Load Data from File");
        AppendMenu(hOptionsMenu, MF_STRING, 40, L"Save Data To File");
        AppendMenu(hOptionsMenu, MF_POPUP, (UINT_PTR)hRectClipMenu, L"Clipping (Rectangle)");
        AppendMenu(hOptionsMenu, MF_POPUP, (UINT_PTR)hSquareClipMenu, L"Clipping (Square)");
        AppendMenu(hOptionsMenu, MF_POPUP, (UINT_PTR)hCircleClipMenu, L"Clipping (Circle)");

        // Main Menu
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hDrawMenu, L"Draw");
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hColorMenu, L"Color");
        AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hOptionsMenu, L"Options");

        SetMenu(hwnd, hMenu);
    }
                  break;

    case WM_COMMAND: {
        switch (LOWORD(wParam)) {
        case 1: selectedShape = 1; clickCount = 0; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 2: selectedShape = 2; clickCount = 0; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 3: selectedShape = 3; clickCount = 0; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 4: selectedShape = 4; clickCount = 0; currentCircleAlgorithm = DIRECT_CIRCLE; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 5: selectedShape = 5; clickCount = 0; currentCircleAlgorithm = POLAR; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 6: selectedShape = 6; clickCount = 0; currentCircleAlgorithm = ITERATIVE_POLAR; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 7: selectedShape = 7; clickCount = 0; currentCircleAlgorithm = MIDPOINT_CIRCLE; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 8: selectedShape = 8; clickCount = 0; currentCircleAlgorithm = MODIFIED_MIDPOINT; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 9: selectedShape = 9; shapeDrawn = false; clickCount = 0; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 10: selectedShape = 10; shapeDrawn = false; clickCount = 0; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 11: selectedShape = 11; clickCount = 0; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 12: selectedShape = 12; clickCount = 0; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 13: selectedShape = 13; clickCount = 0; shapeDrawn = false; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 14: selectedColor = RGB(0, 0, 0); break;
        case 15: selectedColor = RGB(255, 0, 0); break;
        case 16: selectedColor = RGB(0, 255, 0); break;
        case 17: selectedColor = RGB(0, 0, 255); break;
        case 18: selectedColor = RGB(255, 255, 0); break;
        case 19: selectedShape = 19; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 20: selectedShape = 20; ClearWindow(hwnd); break;
        case 21: selectedShape = 21; GetFileName("Enter the File path"); LoadPointsFromFile(hwnd, FileName); splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 23: selectedShape = 23; clickCount = 0; currentEllipseAlgorithm = DIRECT_ELLIPSE; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 24: selectedShape = 24; clickCount = 0; currentEllipseAlgorithm = POLAR_ELLIPSE; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 25: selectedShape = 25; clickCount = 0; currentEllipseAlgorithm = MIDPOINT_ELLIPSE; splinePoints.clear(); clippingMode = NONE; fillingMode = false; break;
        case 26: selectedShape = 26; clippingMode = POINT_CLIP; points.clear(); originalPoints.clear(); clippedPolygonResult.clear(); showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click a point, then right-click to clip it", L"Point Clipping", MB_OK); break;
        case 27: selectedShape = 27; clippingMode = LINE_CLIP; points.clear(); originalPoints.clear(); clippedPolygonResult.clear(); showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click two points for line, then right-click to clip it", L"Line Clipping", MB_OK); break;
        case 28: selectedShape = 28; clippingMode = POLYGON_CLIP; points.clear(); originalPoints.clear(); clippedPolygonResult.clear(); showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click at least 3 points for polygon, then right-click to clip it", L"Polygon Clipping", MB_OK); break;
        case 29: selectedShape = 29; clippingMode = DEFINE_RECT; points.clear(); originalPoints.clear(); clippingWindow.clear(); clippedPolygonResult.clear(); clippingWindowDefined = false; showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click two points to define rectangle clipping window", L"Define Rectangle Clipping Window", MB_OK); break;
        case 31: selectedShape = 31; clippingMode = POINT_CLIP; points.clear(); originalPoints.clear(); showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click a point, then right-click to clip it", L"Point Clipping", MB_OK); break;
        case 32: selectedShape = 32; clippingMode = LINE_CLIP; points.clear(); originalPoints.clear(); showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click two points for line, then right-click to clip it", L"Line Clipping", MB_OK); break;
        case 33: selectedShape = 33; clippingMode = DEFINE_SQUARE; points.clear(); originalPoints.clear(); clippingSquareDefined = false; clippingSquare = Square(); squareDefinitionStep = 0; showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click to define the center of the square", L"Define Square Clipping Window", MB_OK); break;
        case 35: selectedShape = 35; clippingMode = DEFINE_CIRCLE; points.clear(); originalPoints.clear(); clippingCircleDefined = false; clippingCircle = Circle(); circleDefinitionStep = 0; showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click to define the center of the circle", L"Define Circle Clipping Window", MB_OK); break;
        case 36: selectedShape = 36; clippingMode = POINT_CLIP; points.clear(); originalPoints.clear(); showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click a point, then right-click to clip it", L"Point Clipping", MB_OK); break;
        case 37: selectedShape = 37; clippingMode = LINE_CLIP; points.clear(); originalPoints.clear(); showClippedResult = false; fillingMode = false; InvalidateRect(hwnd, NULL, TRUE); MessageBox(hwnd, L"Click two points for line, then right-click to clip it", L"Line Clipping", MB_OK); break;
        case 38:
            selectedShape = 38;
            fillingMode = true;
            fillColor = GetColorFromConsole("Select fill color for Recursive Flood Fill:");

            break;
        case 39:
            selectedShape = 39;
            fillingMode = true;
            fillColor = GetColorFromConsole("Select fill color for Non-Recursive Flood Fill:");

            break;
        case 30:
            selectedShape = 30;
            fillPolygonPoints.clear();
            fillingMode = true;
            fillColor = GetColorFromConsole("Select fill color for Convex Fill:");
            borderColor = GetColorFromConsole("Select border color for Convex Fill:");
            break;
        case 34:
            selectedShape = 34;
            fillPolygonPoints.clear();
            fillingMode = true;
            fillColor = GetColorFromConsole("Select fill color for Non-Convex Fill:");
            borderColor = GetColorFromConsole("Select border color for Non-Convex Fill:");
            break;
        case 40:
            GetFileName("Enter the file path:");
            SavePointsToFile(hwnd, Saved, FileName);
        }
        break;
    }

    case WM_LBUTTONDOWN: {
        hdc = GetDC(hwnd);
        pt.x = LOWORD(lParam);
        pt.y = HIWORD(lParam);

        // Handle flood fill modes
        if (selectedShape == 38 || selectedShape == 39) {
            COLORREF targetColor = GetPixel(hdc, pt.x, pt.y);
            if (selectedShape == 38) {
                FloodFillRecursive(hdc, pt.x, pt.y, fillColor, targetColor);
            }
            else {
                FloodFillNonRecursive(hdc, pt.x, pt.y, fillColor, targetColor);
            }
            ReleaseDC(hwnd, hdc);
            InvalidateRect(hwnd, NULL, FALSE); // Changed to FALSE to prevent screen clearing
            break;
        }

        // Handle filling mode
        if (fillingMode && (selectedShape == 30 || selectedShape == 34)) {
            fillPolygonPoints.push_back(pt);
            if (fillPolygonPoints.size() > 1) {
                POINT prev = fillPolygonPoints[fillPolygonPoints.size() - 2];
                drawDDA(hdc, prev.x, prev.y, pt.x, pt.y, borderColor);
            }
            else {
                drawPoint(hdc, pt.x, pt.y, borderColor);
            }
            ReleaseDC(hwnd, hdc);
            break;
        }

        if (clippingMode == DEFINE_RECT) {
            if (points.size() < 2) {
                points.push_back(pt);
                drawPoint(hdc, pt.x, pt.y, RGB(0, 0, 255));
                if (points.size() == 2) {
                    clippingWindow = points;
                    clippingWindowDefined = true;
                    points.clear();
                    clippingMode = NONE;
                    InvalidateRect(hwnd, NULL, TRUE);
                    MessageBox(hwnd, L"Rectangle clipping window defined successfully!", L"Rectangle Definition Complete", MB_OK);
                }
            }
            ReleaseDC(hwnd, hdc);
            break;
        }
        else if (clippingMode == DEFINE_SQUARE) {
            if (squareDefinitionStep == 0) {
                clippingSquare.centerX = pt.x;
                clippingSquare.centerY = pt.y;
                clippingSquare.size = 0;
                squareDefinitionStep = 1;
                drawPoint(hdc, pt.x, pt.y, RGB(0, 0, 255));
                InvalidateRect(hwnd, NULL, false);
                MessageBox(hwnd, L"Center defined! Now click to set the size of the square.", L"Square Definition", MB_OK);
            }
            else if (squareDefinitionStep == 1) {
                double dist = sqrt(pow(pt.x - clippingSquare.centerX, 2) + pow(pt.y - clippingSquare.centerY, 2));
                clippingSquare.size = Round(dist * 2.0);
                clippingSquare.updateBounds();
                clippingSquareDefined = true;
                squareDefinitionStep = 0;
                clippingMode = NONE;
                InvalidateRect(hwnd, NULL, false);
                MessageBox(hwnd, L"Square clipping window defined successfully!", L"Square Definition Complete", MB_OK);
            }
            ReleaseDC(hwnd, hdc);
            break;
        }
        else if (clippingMode == DEFINE_CIRCLE) {
            if (circleDefinitionStep == 0) {
                clippingCircle.centerX = pt.x;
                clippingCircle.centerY = pt.y;
                clippingCircle.radius = 0;
                circleDefinitionStep = 1;
                drawPoint(hdc, pt.x, pt.y, RGB(0, 0, 255));
                InvalidateRect(hwnd, NULL, false);
                MessageBox(hwnd, L"Center defined! Now click to set the radius of the circle.", L"Circle Definition", MB_OK);
            }
            else if (circleDefinitionStep == 1) {
                double dist = sqrt(pow(pt.x - clippingCircle.centerX, 2) + pow(pt.y - clippingCircle.centerY, 2));
                clippingCircle.radius = dist;
                clippingCircleDefined = true;
                circleDefinitionStep = 0;
                clippingMode = NONE;
                InvalidateRect(hwnd, NULL, TRUE);
                MessageBox(hwnd, L"Circle clipping window defined successfully!", L"Circle Definition Complete", MB_OK);
            }
            ReleaseDC(hwnd, hdc);
            break;
        }
        else if (clippingMode == POINT_CLIP || clippingMode == LINE_CLIP || clippingMode == POLYGON_CLIP) {
            if ((clippingMode == POINT_CLIP && points.size() >= 1) ||
                (clippingMode == LINE_CLIP && points.size() >= 2) ||
                (clippingMode == POLYGON_CLIP && points.size() >= 3)) {
                MessageBox(hwnd, L"Shape complete! Right-click to clip or clear screen to start over.", L"Info", MB_OK);
            }
            else {
                points.push_back(pt);
                showClippedResult = false;
                drawPoint(hdc, pt.x, pt.y, selectedColor);
                if (clippingMode == LINE_CLIP && points.size() == 2) {
                    drawDDA(hdc, points[0].x, points[0].y, points[1].x, points[1].y, selectedColor);
                    cout << "Line drawn - Start: (" << points[0].x << ", " << points[0].y << "), End: (" << points[1].x << ", " << points[1].y << ")\n";
                }
                else if (clippingMode == POLYGON_CLIP && points.size() >= 2) {
                    drawDDA(hdc, points[points.size() - 2].x, points[points.size() - 2].y,
                        points[points.size() - 1].x, points[points.size() - 1].y, selectedColor);
                }
                InvalidateRect(hwnd, NULL, false);
            }
            ReleaseDC(hwnd, hdc);
            break;
        }

        switch (selectedShape) {
        case 1: case 2: case 3:
            if (clickCount == 0) {
                startPt = pt;
                clickCount = 1;
            }
            else {
                endPt = pt;
                if (selectedShape == 1) drawDDA(hdc, startPt.x, startPt.y, endPt.x, endPt.y, selectedColor);
                else if (selectedShape == 2) DrawLineMidpoint(hdc, startPt.x, startPt.y, endPt.x, endPt.y, selectedColor);
                else if (selectedShape == 3) DrawLineParametric(hdc, startPt.x, startPt.y, endPt.x, endPt.y, selectedColor);
                shapeDrawn = true;
                clickCount = 0;
                cout << "Line drawn - Start: (" << startPt.x << ", " << startPt.y << "), End: (" << endPt.x << ", " << endPt.y << ")\n";
            }
            break;
        case 4: case 5: case 6: case 7: case 8:
            if (clickCount == 0) {
                xc = pt.x;
                yc = pt.y;
                clickCount = 1;
                shapeDrawn = false;
            }
            else {
                R = Round(sqrt(pow(pt.x - xc, 2) + pow(pt.y - yc, 2)));
                switch (currentCircleAlgorithm) {
                case DIRECT_CIRCLE: DrawCircleDirect(hdc, xc, yc, R, selectedColor); break;
                case POLAR: DrawCirclePolar(hdc, xc, yc, R, selectedColor); break;
                case ITERATIVE_POLAR: DrawCircleIterativePolar(hdc, xc, yc, R, selectedColor); break;
                case MIDPOINT_CIRCLE: DrawCircleMidpoint(hdc, xc, yc, R, selectedColor); break;
                case MODIFIED_MIDPOINT: DrawCircleModifiedMidpoint(hdc, xc, yc, R, selectedColor); break;
                }
                shapeDrawn = true;
                clickCount = 0;
                cout << "Circle drawn - Center: (" << xc << ", " << yc << "), Radius: " << R << "\n";
            }
            break;
        case 9: case 10:
            if (!shapeDrawn) {
                if (clickCount == 0) {
                    xc = pt.x;
                    yc = pt.y;
                    clickCount = 1;
                }
                else {
                    R = Round(sqrt(pow(pt.x - xc, 2) + pow(pt.y - yc, 2)));
                    DrawCircleModifiedMidpoint(hdc, xc, yc, R, selectedColor);
                    shapeDrawn = true;
                    clickCount = 0;
                    cout << "Circle drawn - Center: (" << xc << ", " << yc << "), Radius: " << R << "\n";
                }
            }
            else {
                int quarter = DetectQuarter(pt.x, pt.y);
                if (quarter != 0) {
                    if (selectedShape == 9) FillQuarterWithSpacedLines(hdc, quarter, selectedColor, 5);
                    else if (selectedShape == 10) FillQuarterWithCircles(hdc, quarter, selectedColor);
                }
            }
            break;
        case 11: case 12:
            if (clickCount == 0) {
                startPt = pt;
                clickCount = 1;
            }
            else {
                endPt = pt;
                if (selectedShape == 11) FillSquareWithHermiteVertical(hdc, startPt.x, startPt.y, endPt.x, endPt.y, selectedColor);
                else if (selectedShape == 12) FillRectangleWithBezierHorizontal(hdc, startPt.x, startPt.y, endPt.x, endPt.y, selectedColor);
                shapeDrawn = true;
                clickCount = 0;
            }
            break;
        case 13:
            splinePoints.push_back(pt);
            clickCount++;
            DrawCircleModifiedMidpoint(hdc, pt.x, pt.y, 3, selectedColor);
            if (clickCount >= 4) {
                DrawCardinalSpline(hdc, splinePoints, selectedColor, 0.0);
                shapeDrawn = true;
                clickCount = 0;
            }
            break;
        case 23: case 24: case 25:
            if (clickCount == 0) {
                xc = pt.x;
                yc = pt.y;
                clickCount = 1;
                shapeDrawn = false;
            }
            else {
                A = abs(pt.x - xc);
                B = abs(pt.y - yc);
                switch (currentEllipseAlgorithm) {
                case DIRECT_ELLIPSE: DirectEllipse(hdc, xc, yc, A, B, selectedColor); break;
                case POLAR_ELLIPSE: PolarEllipse(hdc, xc, yc, A, B, selectedColor); break;
                case MIDPOINT_ELLIPSE: MidPointEllipse(hdc, xc, yc, A, B, selectedColor); break;
                }
                shapeDrawn = true;
                clickCount = 0;
            }
            break;
        }
        ReleaseDC(hwnd, hdc);
        break;
    }

    case WM_RBUTTONDOWN: {
        // Handle filling mode
        if (fillingMode && (selectedShape == 30 || selectedShape == 34) && fillPolygonPoints.size() >= 3) {
            hdc = GetDC(hwnd);
            // Close the polygon
            POINT first = fillPolygonPoints[0];
            POINT last = fillPolygonPoints[fillPolygonPoints.size() - 1];
            drawDDA(hdc, last.x, last.y, first.x, first.y, borderColor);

            // Print polygon vertices
            cout << "Polygon drawn - Vertices:\n";
            for (size_t i = 0; i < fillPolygonPoints.size(); ++i) {
                cout << "Vertex " << i + 1 << ": (" << fillPolygonPoints[i].x << ", " << fillPolygonPoints[i].y << ")\n";
            }

            vector<Point> polygonD;
            for (const auto& p : fillPolygonPoints) {
                polygonD.push_back(Point(static_cast<double>(p.x), static_cast<double>(p.y)));
            }

            // Perform the filling
            if (selectedShape == 30) {
                ConvexFill(hdc, polygonD.data(), polygonD.size(), fillColor);
            }
            else if (selectedShape == 34) {
                GeneralPolygonFill(hdc, polygonD.data(), polygonD.size(), fillColor);
            }

            fillPolygonPoints.clear();
            fillingMode = false;
            ReleaseDC(hwnd, hdc);
            InvalidateRect(hwnd, NULL, FALSE); // Changed to FALSE to prevent screen clearing
            break;
        }

        if ((selectedShape >= 26 && selectedShape <= 28) ||
            (selectedShape >= 31 && selectedShape <= 32) ||
            (selectedShape >= 36 && selectedShape <= 37)) {
            bool validWindow = false;
            if (selectedShape >= 26 && selectedShape <= 28) {
                validWindow = clippingWindowDefined;
            }
            else if (selectedShape >= 31 && selectedShape <= 32) {
                validWindow = clippingSquareDefined;
            }
            else if (selectedShape >= 36 && selectedShape <= 37) {
                validWindow = clippingCircleDefined;
            }

            if (!validWindow) {
                MessageBox(hwnd, L"Please define a clipping window first!", L"Error", MB_OK);
                break;
            }

            if ((clippingMode == POINT_CLIP && points.empty()) ||
                (clippingMode == LINE_CLIP && points.size() < 2) ||
                (clippingMode == POLYGON_CLIP && points.size() < 3)) {
                MessageBox(hwnd, L"Insufficient points for clipping!", L"Error", MB_OK);
                break;
            }

            hdc = GetDC(hwnd);
            originalPoints = points;
            showClippedResult = true;

            if (selectedShape >= 26 && selectedShape <= 28 && clippingMode == POLYGON_CLIP) {
                int xmin = min(clippingWindow[0].x, clippingWindow[1].x);
                int ymin = min(clippingWindow[0].y, clippingWindow[1].y);
                int xmax = max(clippingWindow[0].x, clippingWindow[1].x);
                int ymax = max(clippingWindow[0].y, clippingWindow[1].y);
                clippedPolygonResult = clipPolygon(points, xmin, ymin, xmax, ymax);
                // Print polygon vertices
                cout << "Polygon clipped - Original vertices:\n";
                for (size_t i = 0; i < points.size(); ++i) {
                    cout << "Vertex " << i + 1 << ": (" << points[i].x << ", " << points[i].y << ")\n";
                }
            }

            InvalidateRect(hwnd, NULL, TRUE);
            ReleaseDC(hwnd, hdc);
        }
        break;
    }

    case WM_LBUTTONUP: {
        if ((selectedShape == 9 || selectedShape == 10) && !shapeDrawn && clickCount == 1) {
            hdc = GetDC(hwnd);
            pt.x = LOWORD(lParam);
            pt.y = HIWORD(lParam);
            R = Round(sqrt(pow(pt.x - xc, 2) + pow(pt.y - yc, 2)));
            DrawCircleModifiedMidpoint(hdc, xc, yc, R, selectedColor);
            shapeDrawn = true;
            clickCount = 0;
            cout << "Circle drawn - Center: (" << xc << ", " << yc << "), Radius: " << R << "\n";
            ReleaseDC(hwnd, hdc);
        }
        break;
    }

    case WM_PAINT: {
        PAINTSTRUCT ps;
        hdc = BeginPaint(hwnd, &ps);

        if (LoadedPoints.size() > 0)
        {
            for (int i = 0;i < LoadedPoints.size();i += 5)
            {
                int x = LoadedPoints[i];
                int y = LoadedPoints[i+1];
                int r = LoadedPoints[i+2];
                int g = LoadedPoints[i+3];
                int b = LoadedPoints[i+4];
                SetPixel(hdc, x, y, RGB(r, g, b));
            }
            LoadedPoints.clear();
        }

        // Draw clipping windows
        if (clippingWindowDefined && clippingWindow.size() == 2) {
            HPEN hPen = CreatePen(PS_DASH, 2, RGB(0, 0, 255));
            HPEN oldPen = (HPEN)SelectObject(hdc, hPen);
            int x1 = clippingWindow[0].x, y1 = clippingWindow[0].y;
            int x2 = clippingWindow[1].x, y2 = clippingWindow[1].y;
            Rectangle(hdc, min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2));
            SelectObject(hdc, oldPen);
            DeleteObject(hPen);
        }
        if (clippingSquareDefined) {
            HPEN hPen = CreatePen(PS_SOLID, 2, RGB(0, 0, 255));
            HPEN oldPen = (HPEN)SelectObject(hdc, hPen);
            drawDDA(hdc, clippingSquare.left, clippingSquare.top, clippingSquare.right, clippingSquare.top, RGB(0, 0, 255));
            drawDDA(hdc, clippingSquare.right, clippingSquare.top, clippingSquare.right, clippingSquare.bottom, RGB(0, 0, 255));
            drawDDA(hdc, clippingSquare.right, clippingSquare.bottom, clippingSquare.left, clippingSquare.bottom, RGB(0, 0, 255));
            drawDDA(hdc, clippingSquare.left, clippingSquare.bottom, clippingSquare.left, clippingSquare.top, RGB(0, 0, 255));
            drawPoint(hdc, clippingSquare.centerX, clippingSquare.centerY, RGB(0, 0, 255));
            SelectObject(hdc, oldPen);
            DeleteObject(hPen);
        }
        else if (clippingMode == DEFINE_SQUARE && squareDefinitionStep == 1) {
            drawPoint(hdc, clippingSquare.centerX, clippingSquare.centerY, RGB(0, 0, 255));
        }

        // Draw circle clipping window
        if (clippingCircleDefined) {
            HPEN hPen = CreatePen(PS_SOLID, 2, RGB(0, 0, 255));
            HPEN oldPen = (HPEN)SelectObject(hdc, hPen);
            DrawCircleMidpoint(hdc, Round(clippingCircle.centerX), Round(clippingCircle.centerY), Round(clippingCircle.radius), RGB(0, 0, 255));
            drawPoint(hdc, Round(clippingCircle.centerX), Round(clippingCircle.centerY), RGB(0, 0, 255));
            SelectObject(hdc, oldPen);
            DeleteObject(hPen);
        }
        else if (clippingMode == DEFINE_CIRCLE && circleDefinitionStep == 1) {
            drawPoint(hdc, Round(clippingCircle.centerX), Round(clippingCircle.centerY), RGB(0, 0, 255));
        }

        // Draw clipping shapes
        if (showClippedResult) {
            if (clippingMode == POINT_CLIP && !originalPoints.empty()) {
                int x = originalPoints[0].x, y = originalPoints[0].y;
                bool inside = false;

                if (selectedShape >= 26 && selectedShape <= 28) {
                    int xmin = min(clippingWindow[0].x, clippingWindow[1].x);
                    int ymin = min(clippingWindow[0].y, clippingWindow[1].y);
                    int xmax = max(clippingWindow[0].x, clippingWindow[1].x);
                    int ymax = max(clippingWindow[0].y, clippingWindow[1].y);
                    inside = isPointInside(x, y, xmin, ymin, xmax, ymax);
                }
                else if (selectedShape >= 31 && selectedShape <= 32) {
                    inside = isPointInside(x, y, clippingSquare.left, clippingSquare.top, clippingSquare.right, clippingSquare.bottom);
                }
                else if (selectedShape >= 36 && selectedShape <= 37) {
                    inside = isPointInsideCircle(x, y, clippingCircle);
                }

                if (inside) drawPoint(hdc, x, y, clippedColor);
            }
            else if (clippingMode == LINE_CLIP && originalPoints.size() >= 2) {
                int cx1, cy1, cx2, cy2;
                bool clipped = false;

                if (selectedShape >= 26 && selectedShape <= 28) {
                    int xmin = min(clippingWindow[0].x, clippingWindow[1].x);
                    int ymin = min(clippingWindow[0].y, clippingWindow[1].y);
                    int xmax = max(clippingWindow[0].x, clippingWindow[1].x);
                    int ymax = max(clippingWindow[0].y, clippingWindow[1].y);
                    clipped = clipLine(originalPoints[0].x, originalPoints[0].y, originalPoints[1].x, originalPoints[1].y,
                        cx1, cy1, cx2, cy2, xmin, ymin, xmax, ymax);
                }
                else if (selectedShape >= 31 && selectedShape <= 32) {
                    clipped = clipLineToSquare(originalPoints[0].x, originalPoints[0].y, originalPoints[1].x, originalPoints[1].y,
                        cx1, cy1, cx2, cy2, clippingSquare);
                }
                else if (selectedShape >= 36 && selectedShape <= 37) {
                    double dcx1, dcy1, dcx2, dcy2;
                    clipped = clipLineToCircle(originalPoints[0].x, originalPoints[0].y, originalPoints[1].x, originalPoints[1].y,
                        dcx1, dcy1, dcx2, dcy2, clippingCircle);
                    cx1 = Round(dcx1);
                    cy1 = Round(dcy1);
                    cx2 = Round(dcx2);
                    cy2 = Round(dcy2);
                }

                if (clipped) drawDDA(hdc, cx1, cy1, cx2, cy2, clippedColor);
            }
            else if (clippingMode == POLYGON_CLIP && selectedShape >= 26 && selectedShape <= 28 && !clippedPolygonResult.empty()) {
                for (size_t i = 0; i < clippedPolygonResult.size(); i++) {
                    POINT p1 = clippedPolygonResult[i];
                    POINT p2 = clippedPolygonResult[(i + 1) % clippedPolygonResult.size()];
                    drawDDA(hdc, p1.x, p1.y, p2.x, p2.y, clippedColor);
                }
            }
        }
        else {
            for (size_t i = 0; i < points.size(); i++) {
                drawPoint(hdc, points[i].x, points[i].y, selectedColor);
                if (clippingMode == LINE_CLIP && i == 1) {
                    drawDDA(hdc, points[0].x, points[0].y, points[1].x, points[1].y, selectedColor);
                }
                else if (clippingMode == POLYGON_CLIP && i > 0) {
                    drawDDA(hdc, points[i - 1].x, points[i - 1].y, points[i].x, points[i].y, selectedColor);
                }
            }
            if (clippingMode == POLYGON_CLIP && points.size() >= 3) {
                drawDDA(hdc, points.back().x, points.back().y, points[0].x, points[0].y, selectedColor);
            }
        }

        // Draw filling polygons
        if (!fillPolygonPoints.empty()) {
            for (size_t i = 0; i < fillPolygonPoints.size(); i++) {
                drawPoint(hdc, fillPolygonPoints[i].x, fillPolygonPoints[i].y, borderColor);
                if (i > 0) {
                    drawDDA(hdc, fillPolygonPoints[i - 1].x, fillPolygonPoints[i - 1].y,
                        fillPolygonPoints[i].x, fillPolygonPoints[i].y, borderColor);
                }
            }
        }

        // Draw other shapes
        if (shapeDrawn) {
            switch (selectedShape) {
            case 4: case 5: case 6: case 7: case 8:
                switch (currentCircleAlgorithm) {
                case DIRECT_CIRCLE: DrawCircleDirect(hdc, xc, yc, R, selectedColor); break;
                case POLAR: DrawCirclePolar(hdc, xc, yc, R, selectedColor); break;
                case ITERATIVE_POLAR: DrawCircleIterativePolar(hdc, xc, yc, R, selectedColor); break;
                case MIDPOINT_CIRCLE: DrawCircleMidpoint(hdc, xc, yc, R, selectedColor); break;
                case MODIFIED_MIDPOINT: DrawCircleModifiedMidpoint(hdc, xc, yc, R, selectedColor); break;
                }
                break;
            case 9: case 10:
                DrawCircleModifiedMidpoint(hdc, xc, yc, R, selectedColor);
                break;
            case 11:
                FillSquareWithHermiteVertical(hdc, startPt.x, startPt.y, endPt.x, endPt.y, selectedColor);
                break;
            case 12:
                FillRectangleWithBezierHorizontal(hdc, startPt.x, startPt.y, endPt.x, endPt.y, selectedColor);
                break;
            case 13:
                if (splinePoints.size() >= 4) {
                    DrawCardinalSpline(hdc, splinePoints, selectedColor, 0.0);
                    for (const auto& point : splinePoints) {
                        DrawCircleModifiedMidpoint(hdc, point.x, point.y, 3, selectedColor);
                    }
                }
                break;
            case 23: case 24: case 25:
                switch (currentEllipseAlgorithm) {
                case DIRECT_ELLIPSE: DirectEllipse(hdc, xc, yc, A, B, selectedColor); break;
                case POLAR_ELLIPSE: PolarEllipse(hdc, xc, yc, A, B, selectedColor); break;
                case MIDPOINT_ELLIPSE: MidPointEllipse(hdc, xc, yc, A, B, selectedColor); break;
                }
                break;
            }
        }

        // Draw instructions on screen
        SetTextColor(hdc, RGB(0, 0, 0));
        SetBkMode(hdc, TRANSPARENT);

        if (clippingMode == DEFINE_CIRCLE) {
            if (circleDefinitionStep == 0) {
                TextOut(hdc, 10, 10, L"Step 1: Click to define circle center", 37);
            }
            else {
                TextOut(hdc, 10, 10, L"Step 2: Click to define circle radius", 37);
            }
        }
        
        else if (fillingMode && (selectedShape == 30 || selectedShape == 34)) {
            TextOut(hdc, 10, 10, L"Mode: Polygon Fill - Click to add vertices, right-click to fill", 63);
        }

        EndPaint(hwnd, &ps);
        break;
    }

    case WM_CLOSE:
        DestroyWindow(hwnd);
        break;

    case WM_DESTROY:
        PostQuitMessage(0);
        break;

    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
    return 0;
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    // Ensure console window is available
    AllocConsole();
    FILE* dummy;
    freopen_s(&dummy, "CONOUT$", "w", stdout);
    freopen_s(&dummy, "CONIN$", "r", stdin);

    WNDCLASS wc = { 0 };
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hIcon = LoadIcon(NULL, IDI_WINLOGO);
    wc.lpszClassName = L"2DDrawingClass";
    if (!RegisterClass(&wc)) {
        MessageBox(NULL, L"Window Registration Failed!", L"Error", MB_ICONERROR | MB_OK);
        return 1;
    }
    HWND hwnd = CreateWindow(
        L"2DDrawingClass",
        L"2D Drawing and Clipping Program",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 800, 600,
        NULL, NULL, hInstance, NULL
    );
    if (!hwnd) {
        MessageBox(NULL, L"Window Creation Failed!", L"Error", MB_ICONERROR | MB_OK);
        return 1;
    }
    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    FreeConsole();
    return (int)msg.wParam;
}
