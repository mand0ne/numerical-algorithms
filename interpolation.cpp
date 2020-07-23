#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
#include <exception>
#include <cmath>

constexpr double PI=4*atan(1.);
using namespace std;

template <typename Type>
bool typeEqual(Type x, Type y){
    if((x<0 && y>0) || (x>0 && y<0))
        return false;
    Type eps = 10 * numeric_limits<Type>::epsilon()* (fabs(x) + fabs(y));
    return fabs(x - y) <= eps;
}

template <typename Type>
bool between(Type x1, Type x2, Type x3){
    return (((x1 < x2) && !typeEqual(x1,x2)) && ((x2 < x3) || typeEqual(x1, x2)));
}

class Comparator{
public:
    bool operator()(const pair<double, double> &p1, const pair<double, double> &p2){
        return p1.first < p2.first;
    }
};

class AbstractInterpolator{
    protected:
        vector<pair<double,double>> data;
        mutable int cachedIndex;
        
        int Locate(double x) const {
            if(cachedIndex != -1) {
                if(between(data.at(cachedIndex).first, x, data.at(cachedIndex + 1).first))
                    return cachedIndex + 1;

                if(cachedIndex >=1 && between(data.at(cachedIndex - 1).first, x, data.at(cachedIndex).first)){
                    cachedIndex--;
                    return cachedIndex + 1;
                }

                if(cachedIndex < data.size() - 2 && between(data.at(cachedIndex + 1).first, x, data.at(cachedIndex + 2).first)){
                    cachedIndex++;
                    return cachedIndex + 1;
                }
            }

            if(x<data.at(0).first || typeEqual(x,data.at(0).first))
                return 0;
            if(x>data.at(data.size() - 1).first && !typeEqual(x,data.at(data.size() - 1).first))
                return data.size();
    
            cachedIndex=distance(data.begin(),lower_bound(data.begin(),data.end(),make_pair(x, numeric_limits<double>::min()))) - 1;
            return cachedIndex + 1;
        }
    
    public:
        AbstractInterpolator(const vector<pair<double,double>> &_data) : data(_data), cachedIndex(-1){
            sort(data.begin(),data.end(),Comparator());
            if(adjacent_find(data.begin(),data.end(),[](const pair<double,double> &p1, const pair<double,double>p2) 
            -> bool {return typeEqual(p1.first,p2.first);})!=data.end())
                throw domain_error("Invalid data set");
        }
    
        virtual double operator()(double x) const =0;
};

class LinearInterpolator : public AbstractInterpolator {
    private:
    static double lineThroughTwoDots(const pair<double,double> &p1, const pair<double,double> &p2, double x) {
        return ((p2.first-x)*p1.second)/(p2.first-p1.first) + ((x-p1.first)*p2.second)/(p2.first-p1.first);    
    }
    
    public:
        LinearInterpolator(const vector<pair<double,double>> &_data) : AbstractInterpolator(_data) {}
        
        double  operator()(double x) const override {
            int i = Locate(x) - 1;
            if(i <= 0) 
                return lineThroughTwoDots(data.at(0),data.at(1),x);
            else if(i>=data.size() - 1) 
                return lineThroughTwoDots(data.at(data.size() - 2),data.at(data.size() - 1),x);
            else 
                return lineThroughTwoDots(data.at(i),data.at(i + 1),x);
        }
};

class PolynomialInterpolator : public AbstractInterpolator {
    vector<double> yn;
public:
    PolynomialInterpolator(const vector<pair<double,double>> &_data) : AbstractInterpolator(_data) {
        yn.resize(data.size());
        yn[0]=data[data.size() - 1].second;
        for(int j(1); j<=data.size() - 1; j++) {
            for(int i(data.size()); i >= j + 1; i--) {
                if(typeEqual(data[i - 1].first, data[i - j - 1].first))
                    throw domain_error("Invalid data set");
                data[i - 1].second=(data[i - 1].second-data[i - 2].second)/(data[i - 1].first-data[i - j - 1].first);
            }
            yn.at(j)=data[data.size() - 1].second;
        }
    }
    
    double operator()(double x) const override{
        double s(data[0].second), p(1);
        for(int j(1); j <= data.size() - 1; j++) {
            p *= (x - data[j - 1].first);
            s += data[j].second*p;
        }
        return s;
    }
    
    void AddPoint(const pair<double,double> &p){
        if(any_of(data.begin(),data.end(),[p](const pair<double,double> &p1)->bool {return typeEqual(p.first,p1.first);})) 
            throw domain_error("Invalid point");
        
        data.push_back(p);
        int n(data.size());
        yn.resize(n);
    
        for(int i(1); i<=n-1; i++) {
            double ynTemp(yn[i - 1]);
            yn[i - 1] = data[n-1].second;
            data[n - 1].second = (data[n - 1].second-ynTemp)/(data[n - 1].first-data[n - i - 1].first);
        }
        yn[n - 1] = data[n - 1].second;
    }
    
    vector<double> GetCoefficients() const{
        int n(data.size());
        vector<double> w(n + 1),p(n);
        w[0]=1;
        for(int i(1);i<=n;i++) {
            for(int j(0);j<i;j++)
                p[j] += data[i - 1].second*w[j];
            
            w[i]=w[i - 1];
            for(int j(i - 1); j>=1; j--)
                w[j]=w[j - 1] - data[i - 1].first*w[j];
            
            w[0] *= -data[i - 1].first;
        }
        return p; 
    }
};

class PiecewisePolynomialInterpolator : public AbstractInterpolator {
    private:
        int k;
    public:
        PiecewisePolynomialInterpolator(const vector<pair<double, double>> &_data, int red) : AbstractInterpolator(_data){
            if(red<1 || red>data.size()) 
                throw domain_error("Invalid order");
            k=red;
        }
        
        double operator()(double x) const override{
            int i(Locate(x)),j,n;
    
            j = (k%2==0)?(i-k/2-1):(i-(k-1)/2-1);
            n = (k%2==0)?(i+k/2):(i+(k+1)/2);
    
            if(j<=0) {
                j=0;
                n=k+1;
            }
    
            if(n>=data.size()) {
                j=data.size()-k-1;
                n=data.size();
            }
    
            double s(0);
            for(int l(j); l<n; l++) {
                double p=data[l].second;
                for(int m(j); m<n; m++)
                    if(l!=m)
                        p *= (x - data[m].first) / (data[l].first - data[m].first);
                s += p;
            }
            return s;
        }
};

class SplineInterpolator : public AbstractInterpolator {
    private:
        std::vector<double> r,q,s;
    public:
        SplineInterpolator(const vector<pair<double,double>> &_data) : AbstractInterpolator(_data) {
            r.resize(data.size());
            s.resize(data.size());
            q.resize(data.size());
            r[0]=0;
            r[data.size()-1]=0;
            
            for(int i(1); i<data.size()-1; i++) {
                s[i] = 2 * (data[i + 1].first - data[i - 1].first);
                r[i] = 3 * ((data[i + 1].second - data[i].second)/(data[i + 1].first - data[i].first) - (data[i].second-data[i - 1].second)/(data[i].first - data[i - 1].first));
            }
            
            for(int i(1); i<data.size()-2; i++) {
                double u = (data[i + 1].first - data[i].first) / s[i];
                s[i + 1] -= u * (data[i + 1].first - data[i].first);
                r[i + 1] -= u * r[i];
            }
        
            r[data.size()-2] /= s[data.size()-2];
            
            for(int i(data.size()-3); i>=1; i--)
            r[i] = (r[i] - (data[i + 1].first - data[i].first) * r[i + 1]) / s[i];
            
            for(int i(0); i<data.size()-1; i++) {
                double x = data[i + 1].first - data[i].first;
                s[i] = (r[i + 1] - r[i]) / (3 * x);
                q[i] = (data[i + 1].second - data[i].second)/x - x*(r[i + 1] + 2 * r[i])/3;
            }
        } 
    
        double operator()(double x) const override {
            int i(Locate(x)-1);
            if(i <= -1) 
                i=0;
            if(i >= data.size()-1) 
                i = data.size()-2;
            double t(x - data[i].first);
            return data[i].second +  t*(q[i] + t*(r[i] + s[i]*t));
        }
};

class BarycentricInterpolator : public AbstractInterpolator {
    private:
        int r;
        vector<double> w;
    public:
        BarycentricInterpolator(const vector<pair<double, double>> &_data, int red) : AbstractInterpolator(_data) {
            if(red<0 || red>data.size()) throw domain_error("Invalid order");
            r=red;
            w.resize(data.size());
            
            for(int i(0); i<data.size(); i++) {
                w[i]=0;
                double p{};
                int m,l;
                
                m=max(1,i-r);
                l=min(i,int(data.size())-r);
                    
                for(int k(m-1); k<l+1; k++) {
                    p=1;
                    for(int j(k); j<k+r; j++)
                        if(j!=i)
                            p = p / (data[i].first-data[j].first);
                            
                    if(k%2==1) p=-p;
                }
                w[i]+=p;
            }
        }
        
        double operator()(double x) const override {
            double p(0),q(0);
            for(int i(0); i<data.size(); i++) {
                if(typeEqual(x,data[i].first)) 
                    return data[i].second;
                double u = w[i] / (x - data[i].first);
                p += u * data[i].second;
                q += u;
            }
            return p/q;
        }
        
        std::vector<double> GetWeights() const { 
            return w; 
        }
};

template <typename FunType>
pair<double,bool> Limit(FunType f, double x0, double h=0, double eps=1e-8, double nMax=20) {
    if(eps<0 || typeEqual(eps,0.) || nMax<3 || nMax>30) 
        throw std::domain_error("Invalid parameters");
    
    if(!isinf(x0) && typeEqual(h,0.))
        h = max(fabs(x0),1.)*0.001;
    
    vector<double> y(nMax,0);
    double yOld(numeric_limits<double>::infinity());
    pair<double,bool> lim(numeric_limits<double>::min(),false);
    
    for(int i(0); i<nMax; i++) {
        if(isinf(x0) && x0<0) y[i] = f(1./(-numeric_limits<double>::min()+h));
        else if(isinf(x0) && x0>0) y[i] = f(1./(0+h));
        else y[i]=f(x0+h);
        double p=2;
        
        for(int k(i-1); k>=0; k--) {
            y[k] = (p * y[k+1] - y[k]) / (p-1);
            p *= 2;
        }
        
        if(fabs(y[0]-yOld)<eps){ 
            lim.second=true;
            break;
        }
        
        yOld = y[0];
        h /= 2;
    }
    
    lim.first=y[0];
    return lim;
}



void classLinearInterpolatorTest() {
        LinearInterpolator lin({{4,16},{2,4},{5,25},{6,36},{1,1},{2.5,6.25},{2.85,8.1225},{3.15,9.9225}});
        if(!typeEqual(lin(2.7),7.32) && !typeEqual(lin(5.8),33.8) && !typeEqual(lin(0.1),-1.7) && !typeEqual(lin(3),9.0225) && !typeEqual(lin(8),58.))
        cout << "LinearInterpolator Test 1 : Error" <<endl;
        else cout << "LinearInterpolator Test 1: Good" << endl;
        
        try{
        LinearInterpolator lin2({{1,13},{2,16},{2,12},{4,20}});
        cout << "LinearInterpolator Test 2: Error" << endl;
        }
        catch(...){
            cout << "LinearInterpolator Test 2: Good" << endl;
        }
}
void classPolynomialInterpolatorTest(){
    PolynomialInterpolator pol1({{-1,0},{-0.5,0.875},{0,1},{0.5,1.125},{1,2}});
    if(typeEqual(pol1(2),9.)) cout << "PolynomialInterpolator Test 1: Good" << endl;
    else cout << "PolynomialInterpolator Test 1: Error" << endl;
    
    PolynomialInterpolator pol( {{-0.439,0} , {-0.409, 0.312} , {-0.418,0.221} , {-0.36,0.756} , {-0.239,1.521} , {0,2} , {0.17,1.876} , {0.29,1.739},{0.42,1.6857},{0.45,1.6963}});
    double fx = 10*0.44*0.44*0.44 -6*0.44*0.44 +2;
    if(fabs(pol(0.44)-fx)<0.005) cout <<"PolynomialInterpolator Test 2: Good" << endl;
    else cout <<"PolynomialInterpolator Test 2: Error" << endl;
    
    try {
        pol1.AddPoint({-1,2});
        cout << "PolynomialInterpolator Test 3: Error" << endl;
    }
    catch (...) {
        cout << "PolynomialInterpolator Test 3: Good" << endl;
    }
    fx = 10*0.5*0.5*0.5 -6*0.5*0.5 +2;
    pol.AddPoint({0.6,2});
    bool good=false;
    if(fabs(pol(0.5)-fx)>0.01) good=true;
    pol.AddPoint({0.58,1.9327});
    if(fabs(pol(0.5)-fx)<0.01 && good) cout << "PolynomialInterpolator Test 4: Good" << endl;
    else cout << "PolynomialInterpolator Test 4: Error" << endl;
}
void classPiecewisePolynomialInterpolatorTest() {
        PiecewisePolynomialInterpolator ppol({{-0.439,0} , {-0.409, 0.312} , {-0.418,0.221} , {-0.36,0.756} , {-0.239,1.521} , {0,2} , {0.17,1.876} , {0.29,1.739},{0.42,1.6857},{0.45,1.6963}},3);
        double fx = 10*0.44*0.44*0.44 -6*0.44*0.44 +2;
        if(fabs(ppol(0.44)-fx)<0.002) cout << "PiecewisePolynomialInterpolator Test 1: Good" << endl;
        else cout <<"PolynomialInterpolator Test 1: Error" << endl;
        
        PiecewisePolynomialInterpolator ppol2({{-1,0},{-0.5,0.875},{0,1},{0.5,1.125},{1,2}},3);
        if(typeEqual(ppol2(2),9.)) cout <<"PiecewisePolynomialInterpolator Test 2: Good" <<endl;
        else cout << "PiecewisePolynomialInterpolator Test 2: Error" << endl;

        try {
            PiecewisePolynomialInterpolator ppol3({{1,320},{2,420},{3,520}},0);
            cout << "PiecewisePolynomialInterpolator Test 3: Error" << endl;
        }
        catch (...) {
           cout << "PiecewisePolynomialInterpolator Test 3: Good" <<endl;
        }
}
void classSplineInterpolatorTest() {
    try{
        vector<pair<double,double>> data;
        SplineInterpolator spl({ {2*PI,0} , {2*PI-PI/8,sin(2*PI-PI/8)} , {2*PI-PI/4,sin(2*PI-PI/4)} , {2*PI-3*PI/8,sin(2*PI-3*PI/8)}, 
                                {2*PI-4*PI/8,sin(2*PI-4*PI/8)}, {2*PI-5*PI/8,sin(2*PI-5*PI/8)}, {2*PI-6*PI/8,sin(2*PI-6*PI/8)},
                                {2*PI-7*PI/8,sin(2*PI-7*PI/8)}, {PI,sin(PI)}, {2*PI-9*PI/8,sin(2*PI-9*PI/8)}, 
                                {2*PI-10*PI/8,sin(2*PI-10*PI/8)}, {2*PI-11*PI/8,sin(2*PI-11*PI/8)}, {2*PI-12*PI/8,sin(2*PI-12*PI/8)},
                                {2*PI-13*PI/8,sin(2*PI-13*PI/8)}, {2*PI-14*PI/8,sin(2*PI-14*PI/8)}, {2*PI-15*PI/8,sin(2*PI-15*PI/8)},
                                {2*PI-16*PI/8,sin(2*PI-16*PI/8)}});
       
       if(fabs(spl(2.2)-sin(2.2))<0.001) cout << "SplineInterpolator Test 1: Good" << endl;
       else cout << "SplineInterpolator Test 1: Error" << endl;
       
       if(!typeEqual((spl(3*PI)-sin(3*PI)),0.)) cout << "SplineInterpolator Test 2: Good" << endl;
       else cout << "SplineInterpolator Test 2: Error " << endl;
    }
    catch(std::domain_error e) {
        std::cout << e.what() << std::endl;
    }
    catch(std::logic_error e) {
        std::cout << e.what() << std::endl;
    }
    catch(...) {
        std::cout << "Nepoznat izuzetak" << std::endl;
    }
}
void classBarycentricInterpolatorTest() {
        vector<pair<double,double>> data;
        BarycentricInterpolator bin({ {2*PI,0} , {2*PI-PI/8,sin(2*PI-PI/8)} , {2*PI-PI/4,sin(2*PI-PI/4)} , {2*PI-3*PI/8,sin(2*PI-3*PI/8)}, 
                                {2*PI-4*PI/8,sin(2*PI-4*PI/8)}, {2*PI-5*PI/8,sin(2*PI-5*PI/8)}, {2*PI-6*PI/8,sin(2*PI-6*PI/8)},
                                {2*PI-7*PI/8,sin(2*PI-7*PI/8)}, {PI,sin(PI)}, {2*PI-9*PI/8,sin(2*PI-9*PI/8)}, 
                                {2*PI-10*PI/8,sin(2*PI-10*PI/8)}, {2*PI-11*PI/8,sin(2*PI-11*PI/8)}, {2*PI-12*PI/8,sin(2*PI-12*PI/8)},
                                {2*PI-13*PI/8,sin(2*PI-13*PI/8)}, {2*PI-14*PI/8,sin(2*PI-14*PI/8)}, {2*PI-15*PI/8,sin(2*PI-15*PI/8)},
                                {2*PI-16*PI/8,sin(2*PI-16*PI/8)}},2);
       
       if(fabs(bin(2.2)-sin(2.2))<0.1) cout << "BarycentricInterpolator Test 1: Good" << endl;
       else cout << "BarycentricInterpolator Test 1: Error" << endl;
      
       if(bin.GetWeights().size()!=0) cout << "BarycentricInterpolator Test 2: Good" << endl;
       else cout << "BarycentricInterpolator Test 2: Error" <<endl;

        BarycentricInterpolator ppol({{-0.439,0} , {-0.409, 0.312} , {-0.418,0.221} , {-0.36,0.756} , {-0.239,1.521} , {0,2} , {0.17,1.876} , {0.29,1.739},{0.42,1.6857},{0.45,1.6963}},3);
        double fx = 10*0.44*0.44*0.44 -6*0.44*0.44 +2;
        if(fabs(ppol(0.44)-fx)<0.004) cout << "BarycentricInterpolator Test 3: Good" << endl;
        else cout <<"BarycentricInterpolator Test 3: Error" << endl;

        try {
            BarycentricInterpolator ppol3({{1,320},{2,420},{3,520}},-1);
            cout << "BarycentricInterpolator Test 4: Error" << endl;
        }
        catch (...) {
           cout << "BarycentricInterpolator Test 4: Good" <<endl;
        }
}
void functionLimitTest() {
    try{
        Limit([](double x){return sin(x)/x;},0,0,-0.0001,20);
        cout << "Limit Test 1: Error" << endl;
    }
    catch(...){
        cout << "Limit Test 1: Good " <<endl;
    }
    
    pair<double,bool> lim(Limit([](double x){return sin(x)/x;},0));
    if(typeEqual(lim.first,1.) && lim.second==true) cout << "Limit Test 2: Good" << endl;
    else cout << "Limit Test 2: Error" << endl;
    
    double inf = std::numeric_limits<double>::infinity();
    lim=Limit([](double x){return (6*x*x + x + 2)/(3*x*x + 1);},inf,1);
    if(fabs(lim.first-2)<0.0000001) cout <<"Limit Test 3: Good" << endl;
    else cout << "Limit Test 3: Error" << endl;
    
    lim=Limit([](double x){return pow((1+1./x),x);},-inf,1);
    if(fabs(lim.first-exp(1))<0.0001) cout <<"Limit Test 4: Good" << endl;
    else cout << "Limit Test 5: Error" << endl;

    lim=Limit([](double x){return (((exp(x) + exp(-x))/2)  - 1)/(x*x);},0,1);
    if(lim.first<0.5001) cout <<"Limit Test 6: Good" << endl;
    else cout << "Limit Test 6: Error" << endl;   
}

int main(){
    classLinearInterpolatorTest();
    classPolynomialInterpolatorTest();
    classPiecewisePolynomialInterpolatorTest();
    classSplineInterpolatorTest();
    classBarycentricInterpolatorTest();
    functionLimitTest();
    return 0;
}
