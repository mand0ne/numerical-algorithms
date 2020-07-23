#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <map>
#include <complex>
#include <ctime>
#include <cstdlib>
#include <iomanip>

using namespace std;

constexpr double pi(3.14159265358979323846);
constexpr double epsilon = numeric_limits<double>::epsilon();

enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

template <typename Type>
bool typeEqual(Type x, Type y)
{
    if((x<0 && y>0) || (x>0 && y<0))
        return false;

    Type eps = 10 * numeric_limits<Type>::epsilon() * (fabs(x) + fabs(y));
    return fabs(x-y)<=eps;
}

template <typename Type>
int sgn(Type x)
{
    if (typeEqual(x, 0.))
        return 0;
    return (0 < x) - (x < 0);
}

template <typename Type>
bool is_zero(Type x, double eps = 1e-12)
{
    return abs(x) <= eps;
}

double randomReal(double min, double max)
{
    double x = (double) rand() / RAND_MAX;
    return min + x * (max - min);
}

complex<double> randomComplex(double rmin, double rmax, double imin, double imax)
{
    return {randomReal(rmin, rmax), randomReal(imin, imax)};
}

bool vectorEqual(vector<complex<double>> v1, vector<complex<double>> v2)
{
    if (v1.size() != v2.size())
        return false;

    for (int i(0); i<v1.size(); i++)
        if (abs(v1[i] - v2[i]) > 0.00001)
            return false;

    return true;
}

pair<complex<double>, bool> Laguerre(vector<complex<double>> p, int n, complex<double> x, double eps, int maxiter)
{
    complex<double> dx = numeric_limits<double>::infinity();
    int k(1);

    while (abs(dx) > eps && k < maxiter) {
        complex<double> f = p.at(n);
        complex<double> d(0), s(0);

        for (int i(n-1); i>=0; i--) {
            s = s*x + 2.0 *d;
            d = d*x + f;
            f = f*x + p.at(i);
        }

        if (abs(f)<=eps || typeEqual(abs(f), eps))
            return {x, true};

        complex<double> r(sqrt(double(n - 1) * (double(n - 1) * d*d - double(n) * f * s)));

        if (abs(d + r) > abs(d - r))
            dx = double(n)*f / (d + r);
        else
            dx = double(n)*f / (d - r);

        x-=dx;
        k++;
    }

    if (abs(dx)<=eps || typeEqual(abs(dx),eps))
        return {x, true};

    return {x, false};
}

pair<complex<double>, bool> Laguerre(vector<double> p, int n, complex<double> x, double eps, int maxiter)
{
    vector<complex<double>> pComplex;

    for(double &d : p)
        pComplex.push_back({d,0});

    return Laguerre(pComplex, n, x, eps, maxiter);
}

template <typename FunType>
bool BracketRoot(FunType f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4)
{
    if (hinit <= 0 || hmax <= 0 || lambda <= 0)
        throw domain_error("Invalid parameters");

    bool found(true);
    double aTemp(x0), bTemp;
    double h = hinit;
    double f1(f(aTemp)), f2(0);

    while (fabs(h) < hmax) {
        bTemp = aTemp + h;
        f2 = f(bTemp);

        while (!isfinite(f2)) {
            h = h / (2 * ( 1 + lambda));

            if (fabs(h) <= fabs(aTemp) * numeric_limits<double>::epsilon()) {
                found = false;
                break;
            }

            bTemp = aTemp + h;
            f2 = f(bTemp);
        }

        if (!found)
            break;


        if (f1*f2 <=00 || typeEqual(f1*f2, 0.)) {
            if (aTemp > bTemp)
                swap(aTemp, bTemp);

            a = aTemp;
            b = bTemp;

            return true;
        }

        h = lambda * h;
        aTemp = bTemp;
        f1 = f2;
    }


    aTemp = x0;
    h = -hinit;
    f1 = f(aTemp);

    while (fabs(h) < hmax) {
        bTemp = aTemp + h;
        f2 = f(bTemp);

        while (!isfinite(f2)) {
            h = h / (2 * ( 1 + lambda));

            if (fabs(h) <= fabs(aTemp) * numeric_limits<double>::epsilon())
                return false;

            bTemp = aTemp + h;
            f2 = f(bTemp);
        }


        if (f1*f2 < 0 || is_zero(f1*f2)) {
            if (aTemp > bTemp)
                swap(aTemp, bTemp);

            a = aTemp;
            b = bTemp;

            return true;
        }

        h = lambda * h;
        aTemp = bTemp;
        f1 = f2;
    }

    return false;
}

template <typename FunType>
double RegulaFalsiSolve(FunType f, double a, double b, RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100)
{
    if (f(a) * f(b) > 0)
        throw range_error("Root must be bracketed");

    if (maxiter <= 0 || eps <= 0)
        throw domain_error("Invalid parameters");

    if (is_zero(f(a), 0.))
        return a;
    if (is_zero(f(b), 0.))
        return b;

    auto varphi = [](double x) {
        return x/(1 + fabs(x));
    };

    bool slavicMode = true;
    double f1, f2, f3(0);

    if (mode == Unmodified || mode == Illinois)
        slavicMode = false;

    f1 = slavicMode?(varphi(f(a))):f(a);
    f2 = slavicMode?(varphi(f(b))):f(b);

    double c(a), cOld(b);

    for (int i(0); i <= maxiter && fabs(c - cOld) > eps; i++) {
        if (i == maxiter)
            throw logic_error("Given accuracy has not achieved");

        cOld = c;
        c = (a*f2 - b*f1) / (f2 - f1);

        f3 = slavicMode?(varphi(f(c))):f(c);

        if (typeEqual(f3, 0.))
            return c;

        if (f1*f3 < 0) {
            b = a;
            f2 = f1;
        } else if (mode == Illinois || mode == IllinoisSlavic)
            f2 = f1*f2 / (f1 + f3); // Pegasus algoritam, nerijetko daje bolje performanse

        a = c;
        f1 = f3;
    }

    return c;
}

template <typename FunType>
double RiddersSolve(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100)
{
    if (f(a) * f(b) > 0)
        throw range_error("Root must be bracketed");

    if (maxiter <= 0 || eps <= 0)
        throw domain_error("Invalid parameters");

    if (is_zero(f(a), 0.))
        return a;

    if (is_zero(f(b), 0.))
        return b;


    double f1(f(a)), f2(f(b));

    for (int i(0); i<=maxiter && fabs(b - a) > eps; i++) {

        if (i == maxiter)
            throw logic_error("Given accuracy has not achieved");

        double c = (a + b)/2;
        double f3 = f(c);

        if (typeEqual(f3, 0.))
            return c;

        double d = c + (f3 * (c - a)*sgn(f1 - f2)) / sqrt(f3*f3 - f1*f2);

        double f4(f(d));

        if (typeEqual(f4, 0.))
            return d;

        if (f3*f4 <= 0 || typeEqual(f3*f4, 0.)) {
            a = c;
            b = d;
            f1 = f3;
            f2 = f4;
        } else if (f1*f4 <= 0 || typeEqual(f1*f4, 0.)) {
            b = d;
            f2 = f4;
        } else {
            a = d;
            f1 = f4;
        }
    }

    return (a + b) / 2;
}

template <typename FunTip1, typename FunTip2>
double NewtonRaphsonSolve(FunTip1 f, FunTip2 fprim, double x0, double eps = 1e-10, double damping = 0, int maxiter = 100)
{
    if (maxiter<=0 || eps<=0 || damping<0 || damping>=1)
        throw domain_error("Invalid parameters");


    double dx = numeric_limits<double>::infinity();
    double x = x0;
    double v(f(x)), d(fprim(x));

    for (int i(0); i<=maxiter && fabs(dx) > eps; i++) {

        if (i == maxiter)
            throw logic_error("Convergence has not achieved");

        if (fabs(v)<=eps || typeEqual(v, eps))
            return x;

        dx = v/d;
        double w(v);
        v = f(x-dx);
        d = fprim(x-dx);

        while (fabs(v) > fabs(w) || !isfinite(v) || typeEqual(d, 0.)) {
            if (typeEqual(damping, 0.))
                throw logic_error("Convergence has not achieved");

            dx = damping*dx;
            v = f(x - dx);
            d = fprim(x - dx);
        }

        x -= dx;
    }

    return x;
}

vector<complex<double>> PolyRoots(vector<complex<double>> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10)
{
    if (eps<=0 || maxiters<=0 || maxtrials<=0)
        throw domain_error("Invalid parameters");

    vector<complex<double>> initialCoefficients(coefficients);
    vector<complex<double>> roots;

    complex<double> x, w, v;

    for (int i(coefficients.size()-1); i>=1; i--) {
        int t(1);
        bool c(false);

        while (!c &&  t < maxtrials) {
            x = randomComplex(-10, 10, -10, 10);
            pair<complex<double>,bool> lc = Laguerre(coefficients, i, x, eps, maxiters);
            x = lc.first;
            c = lc.second;
            t++;
        }

        if (!c)
            throw logic_error("Convergence has not achieved");

        // poliranje //
        // P O L I R A N J E //
        pair<complex<double>,bool> lc = Laguerre(initialCoefficients, initialCoefficients.size()-1, x, eps, maxiters);
        if (lc.second)
            x = lc.first;

        if (fabs(x.imag()) <= eps || typeEqual(fabs(x.imag()), eps))
            x = x.real();

        if (fabs(x.real()) <= eps || typeEqual(fabs(x.real()),eps))
            x = {0, x.imag()};

        roots.push_back(x);

        v = coefficients.at(i);
        for (int j(i-1); j>=0; j--) {
            w = coefficients.at(j);
            coefficients.at(j) = v;
            v = w+x*v;
        }
    }

    sort(roots.begin(), roots.end(), [](complex<double> c1, complex<double> c2) {
        return (c1.real()<c2.real() || (typeEqual(c1.real(), c2.real()) && c1.imag()<c2.imag() && !typeEqual(c1.imag(), c2.imag())));
    });

    return roots;
}

vector<complex<double>> PolyRoots(vector<double> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10)
{
    if (eps<=0 || maxiters<=0 || maxtrials<=0)
        throw domain_error("Invalid parameters");

    vector<double> initialCoefficients(coefficients);
    vector<complex<double>> roots(coefficients.size());
    complex<double> x;

    int i(coefficients.size()-1);
    double v(0), a(0), b(0), u(0), w(0);

    while (i >= 1) {
        int t(1);
        bool c(false);

        while (!c &&  t < maxtrials) {
            x = randomComplex(-10, 10, -10, 10);
            pair<complex<double>,bool> lc = Laguerre(coefficients, i, x, eps, maxiters);
            x = lc.first;
            c = lc.second;
            t++;
        }

        if (!c)
            throw logic_error("Convergence has not achieved");

        // poliranje //
        // P O L I R A N J E //
        pair<complex<double>, bool> lc = Laguerre(initialCoefficients, initialCoefficients.size()-1, x, eps, maxiters);
        if (lc.second)
            x = lc.first;


        if (fabs(x.imag())<=eps || typeEqual(fabs(x.imag()),eps)) {
            x = x.real();
            roots.at(i) = x;
            v = coefficients.at(i);

            for (int j(i-1); j>=0; j--) {
                w = coefficients.at(j);
                coefficients.at(j) = v;
                v = w + x.real() * v;
            }

            i--;
        }

        else {
            roots.at(i) = x;
            roots.at(i-1) = conj(x);
            a = 2 * x.real();
            b = abs(x) * abs(x);
            u = coefficients.at(i);
            v = coefficients.at(i-1) + a*u;
            for (int j(i-2); j>=0; j--) {
                w = coefficients.at(j);
                coefficients.at(j) = u;
                u = v;
                v = w + a*v - b*coefficients.at(j);
            }
            i-=2;
        }

    }

    // vracaju se vrijednosti od i = 1...n
    roots.erase(roots.begin());

    sort(roots.begin(), roots.end(), [](complex<double> c1, complex<double> c2) {
        return (c1.real()<c2.real() || (typeEqual(c1.real(), c2.real()) && c1.imag()<c2.imag() && !typeEqual(c1.imag(), c2.imag())));
    });

    return roots;
}

template <typename FunTip>
void BracketMinimum(FunTip f, double &x, double &h, double hmax, double lambda, double  &a, double &b)
{

    while (fabs(h) < hmax) {

        if (f(x + h) < f(x)) {
            b = x + h;
            a = x - h;
        } else if (f(x - h) < f(x)) {
            b = x - h;
            a = b - h;
        } else {
            a = x - h;
            b = x + h;
            break;
        }
        x = b;
        h *= lambda;
    }

    if (fabs(h) >= hmax)
        throw logic_error("Minimum has not found");
}

template <typename FunTip>
double FindMinimum(FunTip f, double x0, double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4)
{
    if (eps<=0 || hinit<=0 || hmax<=0 || lambda<=0)
        throw domain_error("Invalid parameters");

    double x(x0), h(hinit);
    double a(x-h), b(x+h);

    BracketMinimum(f, x, h, hmax, lambda, a, b);

    // Uklanjanje ograničenja na unimodalnost
    double phi((1 + sqrt(5))/2), c(x); // Zlatna sredina
    double d, u, v;

    if (fabs(c - a) < fabs(b - c))
        d = b - (b - c) / phi;

    else {
        d = c;
        c = a + (c - a) / phi;
    }

    u = f(c);
    v = f(d);

    while (fabs(b - a) > eps) {
        if (u < v) {
            b = d;
            d = c;
            c = a + (c - a) / phi;
            v = u;
            u = f(c);
        } else {
            a = c;
            c = d;
            d = b -(b - d) / phi;
            u = v;
            v = f(d);
        }
    }

    return (a + b) / 2;
}

template <typename FunTip>
double RK4Step(FunTip f, double x, double y, double h)
{
    double K1 = f(x, y);
    double K2 = f(x + h/2, y + (h*K1)/2);
    double K3 = f(x + h/2, y + (h*K2)/2);
    double K4 = f(x + h, y + h*K3);

    return y + (h*(K1 + 2*K2 + 2*K3 + K4))/6;
}

template <typename FunTip>
vector<pair<double, double>> RK4Integrator(FunTip f, double x0, double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false)
{
    vector<pair<double, double>> points;

    double x(x0), y(y0);
    if (h < 0)
        xmax = -1*xmax;

    points.push_back({x, y});
    bool x_max(false);

    while (fabs(x) <= (xmax + eps) && (!adaptive || !x_max)) {

        if (!adaptive) {
            double step = RK4Step(f, x, y, h);
            y = step;
            x += h;
            if (fabs(x) <= (xmax + eps))
                points.push_back({x, y});
        }

        else {
            double u = RK4Step(f, x, y, h/2);
            double v = RK4Step(f, x+h/2, u, h/2);
            double w = RK4Step(f, x, y, h);
            double delta = fabs(w - v) / fabs(h);

            if (delta <= eps) {
                x += h;
                y = v;
                points.push_back({x, y});
            }
            h *= min(5.0, 0.9*pow((eps/delta), 1./4));
        }
    }

    if (h < 0)
        xmax = -1*xmax;

    if ((h > 0) && (points[points.size()-1].first > xmax) && adaptive) {
        h = xmax - points[points.size()-1].first;
        double u = RK4Step(f, x, y, h/2);
        double v = RK4Step(f, x+h/2, u, h/2);
        points[points.size()-1] = {xmax, v};
    } else if ((h < 0) && (points[points.size()-1].first < xmax) && adaptive) {
        h = xmax - points[points.size()-1].first;
        double u = RK4Step(f, x, y, h/2);
        double v = RK4Step(f, x+h/2, u, h/2);
        points[points.size()-1] = {xmax, v};
    }

    if (adaptive)
        points.erase(points.begin());

    return points;
}

void FFT(vector<double> &x, vector<complex<double>> &xt, int N, int s=0, int d=0, int t=1)
{
    if (N==1)
        xt[d] = x[s];
    else {
        FFT(x, xt, N/2, s, d, 2*t);
        FFT(x, xt, N/2, s + t, d + N/2, 2*t);
        complex<double> mi, w, u, v;

        mi = {1, 0};
        w = exp(complex<double>(0, (-2 * pi)/N));

        for(int k=d; k<=(d+N/2-1); k++) {
            u = xt[k];
            v = mi * xt[k + N/2];
            xt[k] = u + v;
            xt[k + N/2] = u - v;
            mi *= w;
        }
    }
}

void invFFT(vector<complex<double>> &xt, vector<complex<double>> &x, int N, int s=0, int d=0, int t=1)
{
    if (N==1)
        x[d] = xt[s];

    else {
        invFFT(xt, x, N/2, s, d, 2*t);
        invFFT(xt, x, N/2, s + t, d + N/2, 2*t);
        complex<double> mi, w, u, v;
        mi = {1, 0};
        w = exp(complex<double>(0, (2*pi)/N));

        for(int k(d); k<=(d + N/2 -1); k++) {
            u = x[k];
            v = mi * x[k + N/2];
            x[k] = (u + v) / 2.;
            x[k + N/2] = (u - v) / 2.;

            mi *= w;
        }
    }
}

vector<double> LossyCompress(vector<double> data, int new_size)
{
    if(new_size<=1 || new_size>data.size())
        throw range_error("Bad new size");

    int N(data.size());
    if(N<=0 || (N & (N - 1)) != 0)
        throw range_error("Data size must be a power of two");

    vector<double> y(data.size());
    for (int i(0); i<data.size()/2; i++)
        y.at(i) = data.at(2 * i);

    for (int i(data.size()/2); i<data.size(); i++)
        y.at(i) = data.at(2 * (data.size() - i) -1);

    vector<complex<double>> y_DFT(y.size());

    FFT(y, y_DFT, y.size());

    vector<double> x_DFT(y_DFT.size());
    for(int i=0; i < y_DFT.size(); i++)
        x_DFT[i] = (exp(complex<double>(0, (-i*pi)/(2 *N))) * y_DFT[i]).real();

    x_DFT.resize(new_size);
    x_DFT[new_size - 1] = data.size();

    return x_DFT;
}

vector<double> LossyDecompress(vector<double> compressed)
{
    int N(compressed[compressed.size()-1]);

    if (N<1 || N<compressed.size())
        throw logic_error("Bad compressed sequence");

    if(N<=0 || (N & (N - 1)) != 0)
        throw range_error("Data size must be a power of two");


    compressed.at(compressed.size()-1) = 0;
    for (int i(compressed.size()); i<N; i++)
        compressed.push_back(0);

    vector<complex<double>> y(N);
    vector<complex<double>> y_temp(N); // pomocna

    y_temp[0] = compressed[0];

    for (int i(1); i<N; i++)
        y_temp.at(i) = 2. * exp(complex<double>(0, (pi*i)/(2*N))) * compressed[i];

    invFFT(y_temp, y, N);

    for (int i(0); i<N; i++) {
        if (i%2 == 0)
            compressed.at(i) = y[i/2].real();
        else
            compressed.at(i) = y[N - (i + 1)/2].real();
    }

    return compressed;
}


    // T E S T //
// F U N C T I O N S //

void BracketRootTest()
{
    bool result;
    double a(0), b(0);

    auto f1 = [](double x) {
        return x * (x+1);
    };
    auto f2 = [](double x) {
        return (x+3)*(x+3);
    };
    auto f3 = [](double x) {
        return x*x + 3*x - 15;
    };


    try {
        result = BracketRoot(f1, 3, a, b, -1);
        cout << "BracketRoot Test 1: Error" << endl;
    } catch(domain_error) {
        cout << "BracketRoot Test 1: Good" << endl;
    } catch(...) {
    }


    result = BracketRoot(f1, 2, a, b);
    if (result)
        cout << "BracketRoot Test 2: Good" << endl;
    else
        cout << "BracketRoot Test 2: Error" << endl;

    result = BracketRoot(f2, -6, a, b);
    if (!result)
        cout << "BracketRoot Test 3: Good" << endl;
    else
        cout << "BracketRoot Test 3: Error" << endl;

    result = BracketRoot(f3, -1, a, b);
    if (result)
        cout << "BracketRoot Test 4: Good" << endl;
    else
        cout << "BracketRoot Test 5: Error" << endl;
}

void RegulaFalsiSolveTest()
{
    double a(0), b(1), result(0);

    auto f2 = [](double x) {
        return (x+2)*(x-8);
    };
    auto f3 = [](double x) {
        return sin(x);
    };

    try {
        RegulaFalsiSolve(f2, a, b);
        cout << "RegulaFalsiSolve Test 1: Error" << endl;
    } catch(range_error) {
        cout << "RegulaFalsiSolve Test 1: Good" << endl;

    } catch(...) {}


    try {
        a=-3;
        b=-1;
        RegulaFalsiSolve(f2, a, b, IllinoisSlavic, -1);
        cout << "RegulaFalsiSolve Test 2: Error" << endl;
    } catch(...) {
        cout << "RegulaFalsiSolve Test 2: Good" << endl;
    }

    try {
        a = -2.3;
        b = -1.9;
        result = RegulaFalsiSolve(f2, a, b, Unmodified);
        if (fabs(result + 2) < 0.0001)
            cout << "RegulaFalsiSolve Test 3: Good" << endl;
        else
            throw range_error("Greska!");
    } catch(range_error) {
        cout << "RegulaFalsiSolve Test 3: Error" << endl;

    } catch(...) {}


    a=-2.3;
    b=-1.9;
    result=RegulaFalsiSolve(f2, a, b, Slavic);
    if (typeEqual(-2.0, result))
        cout << "RegulaFalsiSolve Test 4: Good" << endl;
    else
        cout << "RegulaFalsiSolve Test 4: Error" << endl;

    a=7.95;
    b=8.05;
    result=RegulaFalsiSolve(f2, a, b, Illinois);
    if (typeEqual(8.0, result))
        cout << "RegulaFalsiSolve Test 5: Good" << endl;
    else
        cout << "RegulaFalsiSolve Test 5: Error" << endl;


    result=RegulaFalsiSolve(f3, 0, 0.01, IllinoisSlavic);
    if (typeEqual(0.0, result))
        cout << "RegulaFalsiSolve Test 6: Good" << endl;
    else
        cout << "RegulaFalsiSolve Test 6: Error" << endl;
}

void PolyRootsTest()
{
    vector<complex<double>> v1({{169, 0}, {10, 0}, {1, 0}});
    vector<complex<double>> nule1({{-5, -12}, {-5, 12}});

    vector<double> v2({153, -145, -9, 1});
    vector<complex<double>> nule2({{-9, 0}, {1, 0}, {17, 0}});

    vector<complex<double>> v3({{-42, -9}, {29, -6}, {-8, 3}, 1});
    vector<complex<double>> nule3({{1, 2}, 3, {4, -5}});

    try {
        PolyRoots(v1, -1);
        cout << "PolyRoots Test 1: Error" << endl;
    } catch(domain_error) {
        cout << "PolyRoots Test 1: Good" << endl;
    } catch(...) {}

    try {
        PolyRoots(v2, -1);
        cout << "PolyRoots Test 1: Error" << endl;
    } catch(domain_error) {
        cout << "PolyRoots Test 2: Good" << endl;
    } catch(...) {}

    try {
        auto rez1 = PolyRoots(v1);
        auto rez2 = PolyRoots(v2);

        if (vectorEqual(rez1, nule1))
            cout << "PolyRoots Test 3: Good" << endl;
        else
            cout << "PolyRoots Test 3: Error" << endl;

        if (vectorEqual(rez2, nule2))
            cout << "PolyRoots Test 4: Good" << endl;
        else
            cout << "PolyRoots Test 4: Error" << endl;
    } catch(...) {}

}

void RiddersSolveTest()
{
    double a(0), b(1), result(0);

    auto f2 = [](double x) {
        return (x+2)*(x-8);
    };
    auto f3 = [](double x) {
        return sin(x);
    };

    try {
        RiddersSolve(f2, a, b);
        cout << "RiddersSolve Test 1: Error" << endl;
    } catch(...) {
        cout << "RiddersSolve Test 1: Good" << endl;
    }

    try {
        a=-3;
        b=-1;
        RiddersSolve(f2, a, b, 1e-10, -1);
        cout << "RiddersSolve Test 2: Error" << endl;
    } catch(...) {
        cout << "RiddersSolve Test 2: Good" << endl;
    }

    try {
        a=-2.2;
        b=-1.8;
        result=RiddersSolve(f2, a, b, IllinoisSlavic);
        if (fabs(result+2)<0.00001)
            cout << "RiddersSolve Test 2: Good" << endl;
        else
            cout << "RiddersSolve Test 2: Error" << endl;

        a=1.3;
        b=4.99;
        result=RiddersSolve(f3, a, b, IllinoisSlavic);
        if (fabs(result-pi)<0.001)
            cout << "RiddersSolve Test 3: Good" << endl;
        else
            cout << "RiddersSolve Test 3: Error" << endl;
    } catch(...) {}
}

void NewtonRaphsonSolveTest()
{
    double result(0);

    auto f1 = [](double x) {
        return log(x);
    };
    auto f1prim = [](double x) {
        return 1./x;
    };

    auto f2 = [](double x) {
        return (x+2)*(x-8);
    };
    auto f2prim=[](double x) {
        return 2*x-6;
    };

    auto f3=[](double x) {
        return sin(x);
    };
    auto f3prim=[](double x) {
        return cos(x);
    };

    try {
        NewtonRaphsonSolve(f2, f2prim, 3, -10);
        cout << "NewtonRaphsonSolve Test 1: Error" << endl;
    } catch(...) {
        cout << "NewtonRaphsonSolve Test 1: Good" << endl;
    }

    try {
        result=NewtonRaphsonSolve(f1, f1prim, 3, 1e-10, 0.5);
        if (fabs(result-1)<0.00001)
            cout << "NewtonRaphsonSolve Test 2: Good" << endl;
        else
            cout << "NewtonRaphsonSolve Test 2: Error" << endl;

        result=NewtonRaphsonSolve(f2, f2prim, 6, 1e-10, 0.5);
        if (fabs(result+2)<0.00001 || fabs(result-8)<0.00001)
            cout << "NewtonRaphsonSolve Test 3: Good" << endl;
        else
            cout << "NewtonRaphsonSolve Test 3: Error" << endl;  ;

        result=NewtonRaphsonSolve(f3, f3prim, 2.6);
        if (fabs(result/pi)-1 < 0.00001)
            cout << "NewtonRaphsonSolve Test 3: Good" << endl;
        else
            cout << "NewtonRaphsonSolve Test 3: Error" << endl;
    } catch(...) {}
}

void FindMinimumTest()
{
    auto f=[](double x) {
        return x*x + 2*x + 1;
    };

    try {
        FindMinimum(f, 1, -1);
        cout << "FindMinimum Test 1: Error" << endl;
    } catch(...) {
        cout << "FindMinimum Test 1: Good" << endl;
    }

    double min = FindMinimum(f, -1.5);
    if (fabs(min + 1) < 0.00001)
        cout << "FindMinimum Test 2: Good" << endl;
    else
        cout << "FindMinimum Test 2: Error" << endl;

}

void RK4IntegratorTest()
{
    int counter(0);

    auto equation=[](double x, double y) {
        return 2*x - y + 7;
    };
    // y(0)=1
    auto solution=[](double x) {
        return 2*x - 4*exp(-x) + 5;
    };

    vector<pair<double, double>> v = RK4Integrator(equation, 0, 1, 1.5, 0.1);

    int n1(v.size());
    for (int i(0); i<n1; i++)
        if (fabs(v[i].second - solution(v[i].first)) < 0.00001)
            counter++;


    v = RK4Integrator(equation, 0, 1, 1.5, 0.1, 1e-8, true);

    int n2(v.size());
    for (int i=0; i<n2; i++)
        if (fabs(v[i].second - solution(v[i].first)) < 0.00001)
            counter++;

    if (counter == n1 + n2)
        cout << "RK4Integrator General Test: Good" << endl;
    else
        cout << "RK4Integrator General Test: Error" << endl;
}

void LossyCompressDecompressTest()
{
    vector<double> seq(32), comp, decomp;
    for (int i(0); i<seq.size(); i++) seq[i] = rand()%10;

    try {
        LossyCompress(seq, 20);
        cout << "LossyCompressDecompress Test 1: Good" << endl;
    } catch(...) {
        cout << "LossyCompressDecompress Test 1: Error" << endl;
    }

    comp = LossyCompress(seq, 32);
    decomp = LossyDecompress(comp);

    bool error(false);
    double MaxError = 1;
    for (int i(0); i<seq.size(); i++)
        if (fabs(seq[i] - decomp[i]) > MaxError)
            error = true;

    if (error)
        cout << "LossyCompressDecompress Test 2: Error" << endl;
    else
        cout << "LossyCompressDecompress Test 2: Good" << endl;
}



int main (){
    BracketRootTest();
    cout << endl;
    RegulaFalsiSolveTest();
    cout << endl;
    RiddersSolveTest();
    cout << endl;
    NewtonRaphsonSolveTest();
    cout << endl;
    PolyRootsTest();
    cout << endl;
    FindMinimumTest();
    cout << endl;
    RK4IntegratorTest();
    cout << endl;
    LossyCompressDecompressTest();

    return 0;
}
