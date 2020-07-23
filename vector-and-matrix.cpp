#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <exception>
#include <iomanip>
#include <ios>
#include <initializer_list>
#include <stdlib.h>
#include <time.h>

using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::initializer_list;

/////////////////////
//// V E C T O R ///
////////////////////
class Vector
{
    vector<double> v;
public:
    explicit Vector(int n)
    {
        if(n<=0) throw std::range_error("Bad dimension");
        v.resize(n,0);
    }

    Vector(std::initializer_list<double> l)
    {
        if(l.size()<=0) throw std::range_error("Bad dimension");
        for(double i: l) v.push_back(i);
    }

    int NElems() const
    {
        return v.size();
    }

    double &operator[](int i)
    {
        return v[i];
    }

    double operator[](int i) const
    {
        return v[i];
    }

    double &operator()(int i)
    {
        if(i>v.size() || i<1) throw std::range_error("Invalid index");
        return v.at(i-1);
    }

    double operator()(int i) const
    {
        if(i>v.size() || i<1) throw std::range_error("Invalid index");
        return v.at(i-1);
    }

    double Norm() const
    {
        double sum{};
        for(int i=0; i<v.size(); i++)
            sum+= v.at(i)*v.at(i);
        return sqrt(sum);
    }

    friend double VectorNorm(const Vector &v);

    double GetEpsilon() const
    {
        return std::numeric_limits<double>::epsilon()*10*this->Norm();
    }

    void Print(char separator = '\n', double eps = -1) const
    {
        if(eps<0) eps=this->GetEpsilon();
        for(int i(0); i<v.size(); i++) {
            if(fabs(v.at(i))<eps) cout << 0;
            else cout << v.at(i);
            if(i!=v.size()-1)
                cout << separator;
        }
        if(separator=='\n') cout << separator;
    }

    friend void PrintVector(const Vector &v, char separator = '\n',double eps = -1);

    friend Vector operator +(const Vector &v1, const Vector &v2);
    Vector &operator +=(const Vector &v)
    {
        *this=*this+v;
        return *this;
    }

    friend Vector operator -(const Vector &v1, const Vector &v2);
    Vector &operator -=(const Vector &v)
    {
        *this=*this-v;
        return *this;
    }

    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    Vector &operator *=(double s)
    {
        *this=*this*s;
        return *this;
    }

    friend double operator *(const Vector &v1, const Vector &v2);

    friend Vector operator /(const Vector &v, double s);
    Vector &operator /=(double s)
    {
        *this=*this/s;
        return *this;
    }

    auto begin() ->decltype(v.begin())
    {
        return v.begin();
    }

    auto end() -> decltype(v.end())
    {
        return v.end();
    }
};

double VectorNorm(const Vector &v)
{
    return v.Norm();
}
void PrintVector(const Vector &v, char separator,double eps)
{
    v.Print(separator,eps);
}
Vector operator +(const Vector &v1, const Vector &v2)
{
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i=1; i<=v3.NElems(); i++)
        v3(i)=v1(i)+v2(i);
    return v3;
}
Vector operator -(const Vector &v1, const Vector &v2)
{
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i=1; i<=v3.NElems(); i++)
        v3(i)=v1(i)-v2(i);
    return v3;
}
Vector operator *(const Vector &v, double s)
{
    Vector v1(v.NElems());
    for(int i=1; i<=v.NElems(); i++)
        v1(i)=v(i)*s;
    return v1;
}
Vector operator *(double s, const Vector &v)
{
    return v*s;
}
double operator *(const Vector &v1, const Vector &v2)
{
    double product {};
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    for(int i=1; i<=v1.NElems(); i++)
        product+=v1(i)*v2(i);
    return product;
}
Vector operator /(const Vector &v, double s)
{
    if(s<std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
    Vector v1(v.NElems());
    for(int i=1; i<=v.NElems(); i++)
        v1(i)=v(i)/s;
    return v1;
}

/////////////////////
//// M A T R I X ///
////////////////////
class Matrix
{
    vector<vector<double>> matrix;

public:
    Matrix(int m, int n)
    {
        if(m<=0 || n<=0) throw std::range_error("Bad dimension");
        matrix.resize(m);
        for(int i=0; i<m; i++) matrix.at(i).resize(n);
    }
    Matrix(const Vector &v)
    {
        if(v.NElems()==0) throw std::range_error("Bad dimension");
        matrix.resize(v.NElems());
        for(int i=0; i<v.NElems(); i++)
            matrix.at(i).at(0)=v(i);
    }
    Matrix(std::initializer_list<std::vector<double>> l)
    {
        if(l.size()==0) throw std::range_error("Bad dimension");
        int n=l.begin()->size();
        for(auto i: l) if(i.size()==0) throw std::range_error("Bad dimension");
            else if(i.size()!=n) throw std::logic_error("Bad matrix");
        matrix.resize(l.size());
        for(int i=0; i<matrix.size(); i++) matrix.at(i).resize(n);
        int j=0;
        for(const vector<double> &i: l) {
            for(int k=0; k<i.size(); k++) {
                matrix.at(j).at(k)=i.at(k);
            }
            j++;
        }
    }

    int NRows() const
    {
        return matrix.size();
    }
    int NCols() const
    {
        return matrix.at(0).size();
    }

    double *operator[](int i)
    {
        return &(matrix.at(i).at(0));
    }
    const double *operator[](int i) const
    {
        return &(matrix.at(i).at(0));
    }

    double &operator()(int i, int j)
    {
        if(i>this->NRows() || i<1 || j<1 || j>this->NCols()) throw std::range_error("Invalid index");
        return matrix.at(i-1).at(j-1);
    }

    double operator()(int i, int j) const
    {
        if(i>this->NRows() || i<1 || j<1 || j>this->NCols()) throw std::range_error("Invalid index");
        return matrix.at(i-1).at(j-1);
    }

    double Norm() const
    {
        double sum{};
        for(auto r: matrix)
            for(auto k: r) sum+=k*k;
        return sqrt(sum);
    }
    friend double MatrixNorm(const Matrix &m);

    double GetEpsilon() const
    {
        return 10*this->Norm()*std::numeric_limits<double>::epsilon();
    }

    void Print(int width = 10, double eps = -1) const
    {
        if(eps<0) eps=this->GetEpsilon();
        for(int i=0; i<matrix.size(); i++) {
            for(int j=0; j<matrix.at(i).size(); j++) {
                std::cout.width(width);
                if(fabs(matrix.at(i).at(j))<eps) cout << 0;
                else cout << matrix.at(i).at(j);
            }
            cout << endl;
        }
    }
    friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1);

    friend Matrix operator +(const Matrix &m1, const Matrix &m2);
    Matrix &operator +=(const Matrix &m)
    {
        *this=*this+m;
        return *this;
    }

    friend Matrix operator -(const Matrix &m1, const Matrix &m2);
    Matrix &operator -=(const Matrix &m)
    {
        *this=*this-m;
        return *this;
    }

    friend Matrix operator *(const Matrix &m, double s);
    friend Matrix operator *(double s, const Matrix &m);
    Matrix &operator *=(double s)
    {
        *this=*this*s;
        return *this;
    }

    friend Matrix operator *(const Matrix &m1, const Matrix &m2);
    Matrix &operator *=(const Matrix &m)
    {
        *this=*this*m;
        return *this;
    }

    friend Vector operator *(const Matrix &m, const Vector &v);

    friend Matrix Transpose(const Matrix &m);
    void Transpose()
    {
        if(this->NRows()!=this->NCols()) {
            Matrix m1(this->NCols(),this->NRows());
            for(int i=0; i<this->NRows(); i++)
                for(int j=0; j<this->NCols(); j++)
                    m1[j][i]=(*this)[i][j];
            *this=m1;
        }

        else {
            int j=1;
            for(int i=0; i<this->NRows(); i++) {
                for(int k=j; k<this->NRows(); k++) {
                    std::swap((*this)[k][i],(*this)[i][k]);
                }
                j++;
            }
        }
    }
};

Matrix Transpose(const Matrix &m)
{
    Matrix m1=m;
    m1.Transpose();
    return m1;
}
double MatrixNorm(const Matrix &m)
{
    return m.Norm();
}
void PrintMatrix(const Matrix &m, int width, double eps)
{
    m.Print(width,eps);
}
Matrix operator +(const Matrix &m1, const Matrix &m2)
{
    if(m1.NRows()!=m2.NRows() || m1.NCols()!=m2.NCols()) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(),m1.NCols());
    for(int i=0; i<m3.NRows(); i++)
        for(int j=0; j<m3.NCols(); j++)
            m3[i][j]=m1[i][j]+m2[i][j];
    return m3;
}
Matrix operator -(const Matrix &m1, const Matrix &m2)
{
    if(m1.NRows()!=m2.NRows() || m1.NCols()!=m2.NCols()) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(),m1.NCols());
    for(int i=0; i<m3.NRows(); i++)
        for(int j=0; j<m3.NCols(); j++)
            m3[i][j]=m1[i][j]-m2[i][j];
    return m3;
}
Matrix operator *(const Matrix &m, double s)
{
    Matrix m1(m.NRows(),m.NCols());
    for(int i=0; i<m1.NRows(); i++)
        for(int j=0; j<m1.NCols(); j++)
            m1[i][j]=m[i][j]*s;
    return m1;
}
Matrix operator *(double s, const Matrix &m)
{
    return m*s;
}
Matrix operator *(const Matrix &m1, const Matrix &m2)
{
    if(m1.NCols()!=m2.NRows()) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(),m2.NCols());
    for(int i = 0; i < m3.NRows(); i++)
        for(int j = 0; j < m3.NCols(); j++) {
            for(int k = 0; k < m1.NCols(); k++) {
                m3[i][j] += m1[i][k] * m2[k][j];
            }
        }

    return m3;
}
Vector operator *(const Matrix &m, const Vector &v)
{
    if(m.NCols()!=v.NElems()) throw std::domain_error("Incompatible formats");
    Vector v1(m.NRows());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < 1; j++) {
            for(int k = 0; k < m.NCols(); k++) {
                v1[i] += m[i][k] * v[k];
            }
        }
    return v1;
}


int main()
{
    {
        // Testiranje nula vektor
        Vector test(5);
        for(auto i : test) if(!(i<=std::numeric_limits<double>::epsilon()*10)) {
                cout << "Greska!";
                return 1;
            }
    }
    {
        try {
            Vector v2(-3);
        } catch(std::exception &ex) {
            cout <<ex.what() << endl;
        }
    }

    {
        try {
            Vector v3 {1,2,3}, v4{1,2}, v(2);
            v=v3+v4;
        } catch (std::domain_error e) {
            std::cout << "'" << e.what() << "'";
        }
    }

    cout << endl;
    Vector v1{1,2,3,4,5,6};
    Vector v2(5);
    for(int i=0; i<5; i++)
        v2[i]=i*log(i+1);

    Vector v3{0.132,-1.2345,3.14159265,991.113,0.0003};
    Vector v4{11./8, 33./241,33.11,42.42,0.01, 874.423};

    cout << "Vektor v1: " << endl;
    v1.Print(' ');
    cout << endl;

    cout << "Vektor v2: " << endl;
    v2.Print(' ');
    cout << endl;

    cout << "Vektor v3: " << endl;
    v3.Print(' ');
    cout << endl;

    cout << "Vektor v4: " << endl;
    v4.Print(' ');
    cout << endl;

    ////
    cout << "Broj elemenata vektora (v1,v2,v3,v4) : " <<endl;
    cout << v1.NElems() <<" " << v2.NElems()<<" " << v3.NElems()<<" " << v4.NElems() << endl;

    // Losa indeksacija
    try {
        cout << v1(-1);
    } catch(std::exception &ex) {
        cout << ex.what() << endl;
    }
    cout << endl << v1(1) << endl;

    try {
        cout << v2(0);
    } catch(std::exception &ex) {
        cout << ex.what() << endl;
    }
    cout << "v2(1) " << v2(1) << endl << endl;


    cout << "Norme vektora (v1,v2,v3,v4) : " <<endl;
    cout << v1.Norm() << " " << v2.Norm() << " " << v3.Norm() << " "<< v4.Norm() << endl;

    {
        cout << "Norma novokreiranog: " << VectorNorm({3.4,21,32.09,99}) << endl;
    }

    cout << "Epsiloni vektora (v1,v2,v3,v4) : " <<endl;
    cout << v1.GetEpsilon() << " " << v2.GetEpsilon() << " " << v3.GetEpsilon() << " "<< v4.GetEpsilon() << endl;

    cout << endl;
    cout << "v1+v4 \n";
    (v1+v4).Print(' ');
    cout << endl;

    cout << endl;
    cout << "v4-v1 \n";
    (v4-v1).Print(' ');
    cout << endl;

    cout << endl;
    cout << "v3+=v2 \n";
    v3+=v2;
    v3.Print(' ');

    cout << endl;
    cout << "v1*v4 \n";
    cout << v1*v4;

    cout << endl;
    cout << "v1*3 \n";
    (v1*3.).Print(' ');

    cout << endl;
    try {
        cout << "v2/0 \n";
        (v1/0).Print(' ');
    } catch(std::exception &ex) {
        cout << ex.what() << endl;
    }

    cout << endl;
    cout << "v2*=3 \n";
    v2*=3;
    v2.Print(' ');

    cout << endl;
    cout << "v2/=3 \n";
    v2/=3;
    v2.Print(' ');

    /////////////////////////////
    // TESTIRANJE KLASE MATRIX //
    ////////////////////////////

    {
        try {
            Matrix m(-1,4);
        } catch(std::exception &ex) {
            cout << ex.what() <<endl;
        }
    }

    cout << endl << endl;
    cout << "// TESTIRANJE KLASE MATRIX //" << endl;
    cout << endl << endl;
    Matrix m1(3,3);
    Matrix m2(2,4);
    Matrix m3(4,2);
    Matrix m4(3,3);
    {
        try {
            Matrix m5{{1,2,3,4},{3,2,1}};
        } catch(std::exception &ex) {
            cout << ex.what() <<endl;
        }
    }

    cout << endl;
    Matrix m5{{3.14159265,2.71828,sqrt(2),0.57721},{0.03213,1.24123,99.321,8.54}};
    m5.Print();

    cout <<endl;
    cout << "Broj redova/kolona m1,m2,m3,m4,m5: " << endl;
    cout << "[" << m1.NRows() << "," << m1.NCols() << "]  "
         << "[" << m2.NRows() << "," << m2.NCols() << "]  "
         << "[" << m3.NRows() << "," << m3.NCols() << "]  "
         << "[" << m4.NRows() << "," << m4.NCols() << "]  "
         << "[" << m5.NRows() << "," << m5.NCols() << "]  " << endl << endl;



    srand (time(NULL));
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++) {
            m1[i][j]=log((rand() %999+1)*26.08);
            m4[i][j]=log((rand() %9999+1)*66.99);
        }

    for(int i=0; i<2; i++)
        for(int j=0; j<4; j++) {
            m2[i][j]=((rand() %999+1)*3.14159)/26.08;
            m3[j][i]=((rand() %9999+1)*0.57*log(2))/66.99;
        }

    cout << "Norme  m1,m2,m3,m4,m5: " << endl;
    cout << m1.Norm() << " " << m2.Norm() << " " << m3.Norm()
         << " " << m4.Norm() << " " << m5.Norm() << endl << endl;

    cout << "m1: \n";
    m1.Print();
    cout << endl << endl;

    cout << "m4: \n";
    m4.Print();
    cout << endl << endl;

    cout <<"m1+m4 \n";
    (m1+m4).Print();
    cout << endl << endl;

    cout << "m2*m3: \n";
    (m2*m3).Print();
    cout << endl << endl;

    cout << "m5: \n";
    m5.Print();
    cout << endl << endl;

    cout << "m2: \n";
    m2.Print();
    cout << endl << endl;

    cout << "m5+=m2: \n";
    m5+=m2;
    m5.Print();
    return 0;
}
