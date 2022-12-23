#include "geometry.h"

bool operator==(const Point& p1, const Point& p2) {
    return (std::abs(p1.x - p2.x) < eps) && (std::abs(p1.y - p2.y) < eps);
}

double cr_prod(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

double sc_prod(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

int sign(double d) {
    return d < 0 ? -1 : 1;
}

bool between(const Point& p, const Point& p1, const Point& p2) {
    return sign(sc_prod(p - p2, p - p1)) == -1;
}

double det(double a11, double a12, double a21, double a22) {
    return a11 * a22 - a12 * a21;
}

void Point::operator+=(const Point& p) {
    x += p.x;
    y += p.y;
}

void Point::operator-=(const Point& p) {
    x -= p.x;
    y -= p.y;
}

Point operator+(const Point& p1, const Point& p2) {
    return {p1.x + p2.x, p1.y + p2.y};
}

Point operator-(const Point& p) {
    return {-p.x, -p.y};
}

Point operator-(const Point& p1, const Point& p2) {
    return {p1.x - p2.x, p1.y - p2.y};
}

void Point::operator*=(double k) {
    x *= k;
    y *= k;
}

void Point::operator/=(double k) {
    x /= k;
    y /= k;
}

Point operator*(const Point& p, double k) {
    return {p.x * k, p.y * k};
}

Point operator*(double k, const Point& p) {
    return {p.x * k, p.y * k};
}

Point operator/(const Point& p, double k) {
    return {p.x / k, p.y / k};
}

double d(const Point& p) {
    return std::sqrt(p.x * p.x + p.y * p.y);
}

double d(const Point& p, const Line& l) {
    return std::abs(l.a * p.x + l.b * p.y + l.c) /
           std::sqrt(l.a * l.a + l.b * l.b);
}

Line::Line(const Point& p1, const Point& p2) {
    a = p2.y - p1.y;
    b = p1.x - p2.x;
    c = p1.x * (p1.y - p2.y) + p1.y * (p2.x - p1.x);
}

Line::Line(double k, double b) : a(k), b(-1), c(b) {}

Line::Line(const Point& p, double k) : a(k), b(-1), c(p.y - k * p.x) {}

Line::Line(double a, double b, double c) : a(a), b(b), c(c) {}

bool operator==(const Line& l1, const Line& l2) {
    return (std::abs(l1.a * l2.b - l1.b * l2.a) < eps) &&
           (std::abs(l1.b * l2.c - l1.c * l2.b) < eps);
}

Line::Line() {}

Line perp(const Point& p, const Line& l) {
    return Line(-l.b, l.a, p.x * l.b - p.y * l.a);
}

Line perp(const Point& p1, const Point& p2) {
    double b = (p1.y - p2.y) / (p1.x - p2.x);
    double c = -0.5 * (p1.x + p2.x + b * (p1.y + p2.y));
    return Line(1, b, c);
}

Point cross(const Line& l1, const Line& l2) {
    double d = det(l1.a, l1.b, l2.a, l2.b);
    double dx = det(-l1.c, l1.b, -l2.c, l2.b);
    double dy = det(l1.a, -l1.c, l2.a, -l2.c);
    return {dx / d, dy / d};
}

Circle Triangle::inscribedCircle() {
    Point center =
        v[0] +
        (d(v[2] - v[0]) * (v[1] - v[0]) + d(v[1] - v[0]) * (v[2] - v[0])) /
            (d(v[1] - v[0]) + d(v[2] - v[1]) + d(v[0] - v[2]));
    double radius = d(center, Line(v[0], v[1]));
    return Circle(center, radius);
}

Circle Triangle::circumscribedCircle() {
    Point center = cross(perp(v[0], v[1]), perp(v[1], v[2]));
    double radius = d(v[0] - center);
    return Circle(center, radius);
}

Point Triangle::centroid() {
    return (v[0] + v[1] + v[2]) / 3;
}

Point Triangle::orthocenter() {
    return cross(perp(v[0], Line(v[1], v[2])), perp(v[1], Line(v[0], v[2])));
}

Line Triangle::EulerLine() {
    return Line(orthocenter(), circumscribedCircle().center());
}

Circle Triangle::ninePointsCircle() {
    Circle ans = circumscribedCircle();
    ans.scale(orthocenter(), 0.5);
    return ans;
}

Triangle::Triangle(const Point& a, const Point& b, const Point& c) {
    v = {a, b, c};
}

Ellipse::Ellipse(const Point& p1, const Point& p2, double r)
    : p1(p1), p2(p2), a(r / 2) {}

std::pair<Point, Point> Ellipse::focuses() const {
    return {p1, p2};
}

double Ellipse::c() const {
    return d(p1 - p2) / 2;
}

double Ellipse::eccentricity() const {
    return c() / a;
}

std::pair<Line, Line> Ellipse::directrices() const {
    return {perp(p2 - p1, center() + 0.5 * (p2 - p1) * c() / a),
            perp(p2 - p1, center() - 0.5 * (p2 - p1) * c() / a)};
}

Point Ellipse::center() const {
    return (p1 + p2) / 2;
}

double Ellipse::perimeter() const {
    return std::comp_ellint_2(eccentricity()) * 4 * a;
}

double Ellipse::area() const {
    return M_PI * a * b();
}

void Ellipse::rotate(const Point& center, double angle) {
    p1.rotate(center, angle);
    p2.rotate(center, angle);
}

void Ellipse::reflect(const Point& center) {
    p1.reflect(center);
    p2.reflect(center);
}

void Ellipse::reflect(const Line& axis) {
    p1.reflect(axis);
    p2.reflect(axis);
}

void Ellipse::scale(const Point& center, double coefficient) {
    p1.scale(center, coefficient);
    p2.scale(center, coefficient);
    a *= std::abs(coefficient);
}

Ellipse::Ellipse() {}

bool Ellipse::containsPoint(const Point& point) const {
    return d(point - p1) + d(point - p2) - 2 * a < eps;
}

double Ellipse::b() const {
    return std::sqrt(a * a - c() * c());
}

Point Ellipse::get_p1() const {
    return p1;
}

Point Ellipse::get_p2() const {
    return p2;
}

double Ellipse::get_a() const {
    return a;
}

bool operator==(const Shape& a, const Shape& b) {
    const auto* Polygon_a = dynamic_cast<const Polygon*>(&a);
    const auto* Ellipse_a = dynamic_cast<const Ellipse*>(&a);
    const auto* Polygon_b = dynamic_cast<const Polygon*>(&b);
    const auto* Ellipse_b = dynamic_cast<const Ellipse*>(&b);
    if ((Polygon_a != nullptr) && (Polygon_b != nullptr)) {
        return (*Polygon_a) == (*Polygon_b);
    }
    if ((Ellipse_a != nullptr) && (Ellipse_b != nullptr)) {
        return (*Ellipse_a) == (*Ellipse_b);
    }
    return false;
}

bool operator==(const Ellipse& a, const Ellipse& b) {
    return ((a.get_p1() == b.get_p1() && a.get_p2() == b.get_p2()) ||
            (a.get_p1() == b.get_p2() && a.get_p2() == b.get_p1())) &&
           (a.get_a() == b.get_a());
}

Polygon::Polygon(std::vector<Point>& a) : v(a) {}

size_t Polygon::verticesCount() {
    return v.size();
}

std::vector<Point> Polygon::getVertices() {
    return v;
}

Polygon::Polygon(){};

double Polygon::perimeter() const {
    double ans = d(v.back() - v[0]);
    for (size_t i = 1; i < v.size(); ++i) {
        ans += d(v[i] - v[i - 1]);
    }
    return ans;
}

double Polygon::area() const {
    double ans = v.back().x * v[0].y - v.back().y * v[0].x;
    for (size_t i = 0; i < v.size() - 1; ++i) {
        ans += v[i].x * v[i + 1].y - v[i].y * v[i + 1].x;
    }
    return std::abs(ans) / 2;
}

void Polygon::rotate(const Point& center, double angle) {
    for (Point& p : v) {
        p.rotate(center, angle);
    }
}

void Polygon::reflect(const Point& center) {
    for (Point& p : v) {
        p.reflect(center);
    }
}

void Polygon::reflect(const Line& axis) {
    for (Point& p : v) {
        p.reflect(axis);
    }
}

void Polygon::scale(const Point& center, double coefficient) {
    for (Point& p : v) {
        p.scale(center, coefficient);
    }
}

bool Polygon::containsPoint(const Point& point) const {
    Point tmp(0.354354, 0.734634);
    Line l(point, tmp);
    int sz = this->sz();
    int count = 0;
    for (int i = 0; i < sz; ++i) {
        Point p = cross(l, Line(v[i], v[mod(i + 1, sz)]));
        if (!between(point, p, tmp) && between(p, v[mod(i + 1, sz)], v[i])) {
            ++count;
        }
    }
    return (count % 2) == 1;
}

bool Polygon::isConvex() const {
    int sz = this->sz();
    int x = sign(cr_prod(v[1] - v[0], v[2] - v[1]));
    for (int i = 0; i < sz; ++i) {
        int y = sign(cr_prod(v[mod(i + 1, sz)] - v[i],
                             v[mod(i + 2, sz)] - v[mod(i + 1, sz)]));
        if (y != x) {
            return false;
        }
    }
    return true;
}

Point Polygon::operator[](int i) const {
    return v[mod(i, sz())];
}

int Polygon::sz() const {
    return static_cast<int>(v.size());
}

int mod(int a, int b) {
    return (a % b + b) % b;
}

bool operator==(const Polygon& a, const Polygon& b) {
    if (a.sz() != b.sz()) {
        return false;
    }
    int sz = b.sz();
    for (int i = 0; i < sz; ++i) {
        if (b[i] == a[0]) {
            if (b[i + 1] == a[1]) {
                for (int j = 0; j < sz; ++j) {
                    if (a[j] != b[j + i]) {
                        return false;
                    }
                }
                return true;
            } else if (b[i - 1] == a[1]) {
                for (int j = 0; j < sz; ++j) {
                    if (a[j] != b[j - i]) {
                        return false;
                    }
                }
                return true;
            }
            return false;
        }
    }
    return false;
}

void Point::rotate(const Point& center, double angle) {
    angle *= M_PI / 180;
    Point p1 = *this - center;
    Point p2{p1.x * cos(angle) - p1.y * std::sin(angle),
             p1.x * sin(angle) + p1.y * std::cos(angle)};
    *this = center + p2;
}

void Point::reflect(const Point& center) {
    *this += 2 * (center - *this);
}

void Point::reflect(const Line& axis) {
    this->reflect(cross(axis, perp(*this, axis)));
}

void Point::scale(const Point& center, double coefficient) {
    *this = center + (*this - center) * coefficient;
}

Point::Point(double x, double y) : x(x), y(y) {}

Point::Point() {}

Circle::Circle(const Point& p, double r) : Ellipse(p, p, 2 * r) {}

double Circle::radius() {
    return a;
}

Rectangle::Rectangle(const Point& p1, const Point& p2, double q) {
    double y = d(p1 - p2) / std::sqrt(q * q + 1);
    Point p = p2;
    p.scale(p1, y / d(p1 - p2));
    p.rotate(p1, atan(q) * 180 / M_PI);
    Point pp = p;
    pp.reflect((p1 + p2) / 2);
    v = {p1, pp, p2, p};
}

Point Rectangle::center() {
    return (v[0] + v[2]) / 2;
}

std::pair<Line, Line> Rectangle::diagonals() {
    return {Line(v[0], v[2]), Line(v[1], v[3])};
}

Rectangle::Rectangle() {}

Circle Square::circumscribedCircle() {
    return Circle(center(), d(center() - v[0]));
}

Circle Square::inscribedCircle() {
    return Circle(center(), d(v[0] - v[1]) / 2);
}

Square::Square(const Point& p1, const Point& p2) : Rectangle(p1, p2, 1) {}

bool isCongruentTo(const Shape& a, const Shape& b) {
    const auto* Polygon_a = dynamic_cast<const Polygon*>(&a);
    const auto* Ellipse_a = dynamic_cast<const Ellipse*>(&a);
    const auto* Polygon_b = dynamic_cast<const Polygon*>(&b);
    const auto* Ellipse_b = dynamic_cast<const Ellipse*>(&b);
    if ((Polygon_a != nullptr) && (Polygon_b != nullptr)) {
        return isCongruentTo(*Polygon_a, *Polygon_b);
    }
    if ((Ellipse_a != nullptr) && (Ellipse_b != nullptr)) {
        return isCongruentTo(*Ellipse_a, *Ellipse_b);
    }
    return false;
}

bool isSimilarTo(const Shape& a, const Shape& b) {
    const auto* Polygon_a = dynamic_cast<const Polygon*>(&a);
    const auto* Ellipse_a = dynamic_cast<const Ellipse*>(&a);
    const auto* Polygon_b = dynamic_cast<const Polygon*>(&b);
    const auto* Ellipse_b = dynamic_cast<const Ellipse*>(&b);
    if ((Polygon_a != nullptr) && (Polygon_b != nullptr)) {
        return isSimilarTo(*Polygon_a, *Polygon_b);
    }
    if ((Ellipse_a != nullptr) && (Ellipse_b != nullptr)) {
        return isSimilarTo(*Ellipse_a, *Ellipse_b);
    }
    return false;
}

bool isCongruentTo(const Ellipse& a, const Ellipse& b) {
    return (a.get_a() == b.get_a()) && (a.b() == b.b());
}

bool isSimilarTo(const Ellipse& a, const Ellipse& b) {
    return a.eccentricity() == b.eccentricity();
}

bool equalAngles(const Polygon& a, const Polygon& b, int k, int sign) {
    int sz = a.sz();
    for (int i = 0; i < sz; ++i) {
        Point pa1 = a[i + 1] - a[i];
        Point pa2 = a[i] - a[i - 1];
        Point pb1 = b[(i + 1 + k) * sign] - b[(i + k) * sign];
        Point pb2 = b[(i + k) * sign] - b[(i - 1 + k) * sign];
        double a1 = sc_prod(pa1, pa2) * d(pb1) * d(pb2);
        double b1 = sc_prod(pb1, pb2) * d(pa1) * d(pa2);
        double a2 = cr_prod(pa1, pa2) * d(pb1) * d(pb2);
        double b2 = cr_prod(pb1, pb2) * d(pa1) * d(pa2);
        if (!((std::abs(a1 - b1) < eps) &&
              ((std::abs(a2 - b2) < eps) || (std::abs(a2 + b2) < eps)))) {
            return false;
        }
    }
    return true;
}

bool similarSides(const Polygon& a, const Polygon& b, int k, bool equal,
                  int sign) {
    double q = d(a[1] - a[0]) / d(b[(k + 1) * sign] - b[k * sign]);
    if (equal && (std::abs(q - 1) > eps)) {
        return false;
    }
    int sz = a.sz();
    for (int i = 0; i < sz; ++i) {
        double a1 = d(a[i + 1] - a[i]);
        double b1 = d(b[(i + 1 + k) * sign] - b[(i + k) * sign]);
        if (std::abs(a1 / b1 - q) > eps) {
            return false;
        }
    }
    return true;
}

bool isSimilarTo(const Polygon& a, const Polygon& b, int k, bool flag,
                 int sign) {
    return equalAngles(a, b, k, sign) && similarSides(a, b, k, flag, sign);
}

bool isSimilarTo(const Polygon& a, const Polygon& b, bool flag) {
    if (a.sz() != b.sz()) {
        return false;
    }
    int sz = a.sz();
    for (int k = 0; k < sz; ++k) {
        if (isSimilarTo(a, b, k, flag, 1) || isSimilarTo(a, b, k, flag, -1)) {
            return true;
        }
    }
    return false;
}

bool isSimilarTo(const Polygon& a, const Polygon& b) {
    return isSimilarTo(a, b, false);
}

bool isCongruentTo(const Polygon& a, const Polygon& b) {
    return isSimilarTo(a, b, true);
}

bool Shape::isCongruentTo(const Shape& b) const {
    return ::isCongruentTo(*this, b);
}

bool Shape::isSimilarTo(const Shape& b) const {
    return ::isSimilarTo(*this, b);
}
