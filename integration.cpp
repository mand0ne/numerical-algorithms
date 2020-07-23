#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::domain_error;


constexpr double pi(3.14159265358979323846);
constexpr double epsilon = std::numeric_limits<double>::epsilon();

template <typename Type>
bool typeEqual(Type x, Type y, double eps = epsilon){
    if((x<0 && y>0) || (x>0 && y<0))
        return false;
    Type eps1 = 10 * eps * (fabs(x) + fabs(y));
    return fabs(x - y) <= eps1;
}

class ChebyshevApproximation {
private:
    vector<double> c;
    int m;
    double min, max;

    ChebyshevApproximation(const vector<double> &coefficients, double xmin, double xmax) : c(coefficients),
                                                                                           m(coefficients.size() - 1),
                                                                                           min(xmin), max(xmax) {}
public:
    template<typename FunctionType>
    ChebyshevApproximation(FunctionType f, double xmin, double xmax, int n){
        if (xmax <= xmin || n < 1)
            throw domain_error("Bad parameters");

        vector<double> w((unsigned int) (n + 2)), v((unsigned int) (n + 1));
        c.resize((unsigned int) (n + 1));
        m = n;
        min = xmin;
        max = xmax;

        for (int i(0); i <= n + 1; i++)
            w[i] = std::cos(pi * i / (2 * n + 2));

        for (int i(0); i <= (int) floor(n / 2); i++)
            v[i] = f((xmin + xmax + (xmax - xmin) * w[2 * i + 1]) / 2);

        for (int i((int) floor(n / 2) + 1); i <= n; i++)
            v[i] = f((xmin + xmax - (xmax - xmin) * w[2 * n + 1 - 2 * i]) / 2);

        for (int k(0); k <= n; k++) {
            double s(0);
            for (int i(0); i <= n; i++) {
                double p = (k * (2 * i + 1)) % (4 * n + 4);
                if (p > 2 * n + 2)
                    p = 4 * n + 4 - p;
                if (p > n + 1)
                    s = s - v[i] * w[static_cast<unsigned int>(2 * n + 2 - p)];
                else
                    s = s + v[i] * w[(unsigned int) p];
            }
            c[k] = 2 * s / (n + 1);
        }
    }

    void set_m(int m) {
        if (m > c.size() - 1 || m <= 1)
            throw domain_error("Bad order");

        this->m = m;
    }

    void trunc(double eps) {
        if (eps < 0)
            throw domain_error("Bad tolerance");

        for (int i(m); i >= 0; i--) {
            if (std::abs(c[i]) > eps) {
                m = i;
                return;
            }
        }

        throw domain_error("Bad tolerance");
    }

    double operator()(double x) const {
        if (x < min || x > max)
            throw domain_error("Bad argument");

        double t = (2 * x - min - max) / (max - min), p = 0, q = c[m];

        for (int k(m - 1); k >= 1; k--) {
            double r(2 * t * q - p + c[k]);
            p = q;
            q = r;
        }

        return t * q - p + c[0] / 2;        
    }

    double derivative(double x) const {
        if (x < min || x > max)
        throw domain_error("Bad argument");

        ChebyshevApproximation derivativeApproximation = this->derivative();
        return derivativeApproximation(x); 
    }

    ChebyshevApproximation derivative() const {
        double u(4. / (max - min));
        
        vector<double> cPrim(c.size());
        cPrim[m - 1] = u * m * c[m];
        cPrim[m - 2] = u * (m - 1) * c[m - 1];
        
        for (int k(m - 3); k >= 0; k--)
            cPrim[k] = cPrim[k + 2] + u * (k + 1) * c[k + 1];

        return ChebyshevApproximation(cPrim, min, max);         
    }

    ChebyshevApproximation antiderivative() const {
        double u((max - min) / 4);
        
        vector<double> temp((unsigned int) (m + 2));
        temp[0] = 0;
        for (int i(1); i <= m - 1; i++)
            temp[i] = (u * (c[i - 1] - c[i + 1])) / i;

        temp[m] = (u * c[m - 1]) / m;
        temp[m + 1] = (u * c[m]) / (m + 1);

        return ChebyshevApproximation(temp, min, max);
    }
    
    double integrate(double a, double b) const {
        if(a > max || a < min || b > max || b < min) 
            throw domain_error("Bad interval");
        
        ChebyshevApproximation F(this->antiderivative());
        return F(b)-F(a);
    }
    
    double integrate() const {
        double s(0);
        
        for(int k(1); k <= (m-1)/2+1; k++)
            s += 2 * c[2*k] / (1 - 4*k*k);
            
        s *= (max - min) / 2;
        s += c[0]*(max - min) / 2;
        return s;
    }
};


template <typename FunctionType>
pair<double,bool> RombergIntegration(FunctionType f, double a, double b, double eps=1e-8, int nMax=1000000, int nMin=50) {
    if(eps<0 || nMin<0 || nMax<0 || nMax<nMin) 
        throw domain_error("Bad parameter");
    
    int N=2;
    double h = (b - a) / N;
    double s = (f(a) + f(b)) / 2;
    double Iold=s;
    vector<double> I;
    for(int i(1); N<=nMax; i++) {
        
        for(int j(1); j<=N/2; j++)
            s += f(a + (2*j - 1) * h);
        
        I.push_back(h*s);
        double p(4);
        for(int k(I.size() - 2); k>=0; k--) {
            I[k] = (p * I[k + 1] - I[k]) / (p - 1);
            p *= 4;
        }
        
        if(N>=nMin && std::abs(I[0]-Iold) <= eps) 
            return {I[0], true};
        
        Iold = I[0];
        h /= 2;
        N *= 2;
    }
    
    return {Iold, false};
}

template <typename FunctionType>
pair<double, bool> TanhSinhIntegration(FunctionType f, double a, double b, double eps = 1e-8, int nMax = 1000000, int nMin = 20, double range = 3.5) {
    if(nMax<nMin || nMin<0 || eps<0 || range<0) 
        throw domain_error("Bad parameter");
        
    int sign(1);
    
    auto fp = [&](double t) { 
                    double t1 = std::cosh(pi * std::sinh(t) / 2);
                    double t2 = (b + a)/2 + (b - a)/2 * std::tanh(pi * std::sinh(t) / 2);
                    double d  = std::cosh(t) * f(t2) / (t1 * t1); 
        
                    if(!std::isfinite(d))
                        return 0.; 
            
                    return (b - a) * pi * d / 4; 
                };
                
    double N(2), h = (2*range)/N, s = (fp(-range) + fp(range))/2, Iold=s;
    
    while(N < nMax) {
        for(int i(1); i<=N/2; i++)
            s += fp(-range + (2*i - 1)*h);
            
        double I(h*s);
        if(N > nMin && std::fabs(I - Iold) <= eps) 
            return {sign*I,true};
        
        Iold = I;
        N *= 2;
        h /= 2;
    }
    
    return {Iold*sign,false};
}

pair<double, bool> operator +(const pair<double,bool> &p1, const pair<double,bool> &p2) { 
    return {p1.first + p2.first, p1.second && p2.second}; 
}

pair<double, bool> operator +=(pair<double,bool> &p1, const pair<double,bool> &p2){
    p1 = p1 + p2;
    return p1;
}

template <typename FunctionType>
pair<double,bool> AdaptiveAux(FunctionType f, double a, double b, double eps, double f1, double f2, double f3, int maxRecursionDepth) {
    if(!std::isfinite(f1)) 
        f1=0;
    
    if(!std::isfinite(f2)) 
        f2=0;
    
    if(!std::isfinite(f3)) 
        f3=0;
        
    double c = (a + b) / 2;
    double I1 = (b - a) * (f1 + 4*f3+ f2) / 6;
    double f4 = f((a + c) / 2); 
    double f5 = f((c + b) / 2);
    
    if(!std::isfinite(f4)) 
        f4=0;
    
    if(!std::isfinite(f5)) 
        f5=0;
    
    double I2 = (b - a) * (f1 + 4*f4 + 2*f3 + 4*f5 + f2) / 12;
    
    if(std::abs(I1-I2)<=eps) 
        return {I2,true};
    
    if(maxRecursionDepth<=0) 
        return {I2,false};
    
    return AdaptiveAux(f, a, c, eps, f1, f3, f4, maxRecursionDepth - 1)  +  AdaptiveAux(f, c, b, eps, f3, f2, f5, maxRecursionDepth - 1);
}

template <typename FunctionType>
pair<double, bool> AdaptiveIntegration(FunctionType f, double a, double b, double eps = 1e-10, int maxRecursionDepth = 30, int nMin = 1) {
    if(nMin<0 || eps<0 || maxRecursionDepth<0) 
        throw domain_error("Bad parameter");
    
    pair<double,bool> s{0,true};

    double h = (b - a) / nMin;

    for(int i(1); i<=nMin; i++) {
        s += AdaptiveAux(f, a, a + h, eps, f(a), f(a + h), f(a + h/2), maxRecursionDepth);
        a += h;
    }

    return s;
}

void testChebyshevApproximation() {
    int counter(0);
    
    auto f = [](double x){ return 25*x*x*x + 31*x*x + 10; };
    
    try {
        ChebyshevApproximation(f, 20, 10, 10);
    }
    catch(domain_error) { 
        counter++; 
    }
    catch(...) {}
    
    ChebyshevApproximation cheby = ChebyshevApproximation(f, 0, 10, 20);
    cheby.trunc(0.00001);
    
    try {
        cheby(11);
    }
    catch(domain_error) { 
        counter++; 
    }
    catch(...) {}
    

    if (typeEqual(2106., cheby(4)))
        counter++;
    if (typeEqual(cheby.derivative(4), 1448.))
        counter++;
    if (typeEqual(cheby.integrate(2, 4), 6926/3., 10))
        counter++;
        
    auto f2 = [](double x){ return cos(x);};
    ChebyshevApproximation cheby2 = ChebyshevApproximation(f2, 0, 15, 3);
    if(typeEqual(cheby2.derivative(2), -sin(2), 10))
        counter++;
    
    if (counter==6) cout << "ChebyshevApproximation Test: Good" << endl;
    else cout << "ChebyshevApproximation Test: Error" << endl;
}

void testRombergIntegration() {
    int counter(0);
    
    auto f=[](double x){ return std::sin(x); };
    auto f3=[](double x){ return exp(x); };
    
    try {
        RombergIntegration(f, 0, pi, -0.001);
    }
    catch(std::domain_error) { 
        counter++; }
    catch(...) {}
    
    
    auto rez = RombergIntegration(f, 0, pi);
    
    if (typeEqual(rez.first, 2., 0.0000000001)) 
        counter++;
    if (rez.second) 
        counter++;
        
    rez = RombergIntegration(f3, 0,3);
    if (typeEqual(rez.first, exp(3)-1, 0.0000000001))
        counter++;
    
    if (counter==4) cout << "RombergIntegration Test: Good" << endl;
    else cout << "RombergIntegration Test: Error" << endl;
}

void testTanhSinhIntegration() { 
    int counter(0);
    
    auto f=[](double x){ return 1/sqrt(1-x*x);};
    
    //Bad parameter
    try {
        RombergIntegration(f, 0, pi, -0.1);
    }
    catch(std::domain_error) { counter++; }
    catch(...) {}
    
    auto rez = RombergIntegration(f, 0.01, 0.99);
    
    if (typeEqual(rez.first, asin(0.99) - asin(0.01), 0.0000000001))
        counter++;
    if (rez.second) 
        counter++;
    
    if (counter==3) std::cout << "TanhSinhIntegration Test: Good" << std::endl;
    else std::cout << "TanhSinhIntegration Test: Error" << std::endl;
}

void testAdaptiveIntegration() {
    int counter(0);
    
    auto f=[](double x){ return 1./std::sqrt(x); };
    
    //Bad parameter
    try {
        AdaptiveIntegration(f, 0, 1, -0.1);
    }
    catch(std::domain_error) { counter++; }
    catch(...) {}
    
    auto rez = AdaptiveIntegration(f, 0, 1);
    
    if (typeEqual(rez.first, 2., 0.0001)) 
        counter++;
    if (!rez.second) 
        counter++;
    
    if (counter==3) std::cout << "AdaptiveIntegration Test: Good" << std::endl;
    else std::cout << "AdaptiveIntegration Test: Error" << std::endl;
}

int main (){
    testChebyshevApproximation();
    testRombergIntegration();
    testTanhSinhIntegration();
    testAdaptiveIntegration();
    return 0;
}