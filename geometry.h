#include <cmath>
#include <iostream>
#include <vector>

const double eps = 1e-6;

struct Line;

struct Point {
    double x, y;
    Point() = default;
    Point(double x, double y);
    void operator+=(const Point& p);
    void operator-=(const Point& p);
    void operator*=(double k);
    void operator/=(double k);
    void rotate(const Point& center, double angle);
    void reflect(const Point& center);
    void reflect(const Line& axis);
    void scale(const Point& center, double coefficient);
};

struct Line {
    double a, b, c;
    Line(double a, double b, double c);
    Line();
    Line(const Point& p1, const Point& p2);
    Line(double k, double b);
    Line(const Point& p, double k);
};

class Shape {
  public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool containsPoint(const Point& point) const = 0;
    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
    bool isCongruentTo(const Shape& b) const;
    bool isSimilarTo(const Shape& b) const;
    virtual ~Shape() = default;
};

class Polygon : public Shape {
  protected:
    std::vector<Point> v;

  public:
    Polygon();
    Polygon(std::vector<Point>& a);
    template <class... Args>
    Polygon(const Point& p, Args... args) : Polygon(args...) {
        v.push_back(p);
    }
    size_t verticesCount();
    std::vector<Point> getVertices();
    bool isConvex() const;
    bool containsPoint(const Point& point) const override;
    double perimeter() const override;
    double area() const override;
    void rotate(const Point& center, double angle) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;
    Point operator[](int i) const;
    int sz() const;
};

class Ellipse : public Shape {
  protected:
    Point p1, p2;
    double a;

  public:
    Ellipse();
    Ellipse(const Point& p1, const Point& p2, double r);
    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;
    double get_a() const;
    Point get_p1() const;
    Point get_p2() const;
    double c() const;
    double b() const;
    bool containsPoint(const Point& point) const override;
    double perimeter() const override;
    double area() const override;
    void rotate(const Point& center, double angle) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;
};

class Circle : public Ellipse {
  public:
    Circle(const Point& p, double r);
    double radius();
};

class Rectangle : public Polygon {
  public:
    Rectangle();
    Rectangle(const Point& p1, const Point& p2, double q);
    Point center();
    std::pair<Line, Line> diagonals();
};

class Square : public Rectangle {
  public:
    Square(const Point& p1, const Point& p2);
    Circle circumscribedCircle();
    Circle inscribedCircle();
};

class Triangle : public Polygon {
  public:
    Triangle(const Point& a, const Point& b, const Point& c);
    Circle circumscribedCircle();
    Circle inscribedCircle();
    Point centroid();
    Point orthocenter();
    Line EulerLine();
    Circle ninePointsCircle();
};

bool operator==(const Line& l1, const Line& l2);

bool operator==(const Point& p1, const Point& p2);

double distance(const Point& p);

Point operator+(const Point& p1, const Point& p2);

Point operator-(const Point& p);

Point operator-(const Point& p1, const Point& p2);

Point operator*(const Point& p, double k);

Point operator*(double k, const Point& p);

Point operator/(const Point& p, double k);

double distance(const Point& p, const Line& l);

Line bisection(const Point& p1, const Point& p2);

Line perpendicular(const Point& p, const Line& l);

Point cross(const Line& l1, const Line& l2);

bool operator==(const Ellipse& a, const Ellipse& b);

bool operator==(const Polygon& a, const Polygon& b);

bool operator==(const Shape& a, const Shape& b);

bool isCongruentTo(const Ellipse& a, const Ellipse& b);

bool isCongruentTo(const Polygon& a, const Polygon& b);

bool isCongruentTo(const Shape& a, const Shape& b);

bool isSimilarTo(const Ellipse& a, const Ellipse& b);

bool isSimilarTo(const Polygon& a, const Polygon& b);

bool isSimilarTo(const Shape& a, const Shape& b);

double cross_product(const Point& a, const Point& b);

double dot_product(const Point& a, const Point& b);

int sign(double d);

bool between(const Point& p, const Point& p1, const Point& p2);

double det(double a11, double a12, double a21, double a22);

bool equalAngles(const Polygon& a, const Polygon& b, int k, int sign);

bool similarSides(const Polygon& a, const Polygon& b, int k, bool equal,
                  int sign);

bool isSimilarTo(const Polygon& a, const Polygon& b, int k, bool flag,
                 int sign);

bool isSimilarTo(const Polygon& a, const Polygon& b, bool flag);
