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

const double epsilon=std::numeric_limits<double>::epsilon();

/////////////////////
//// V E C T O R ///
////////////////////
class Vector{
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
        for(int i=0; i<v.size(); i++) sum+= v.at(i)*v.at(i);
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
        for(int i=0; i<v.size(); i++) {
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

    void Chop(double eps=-1){
    if(eps<0) eps=this->GetEpsilon();
    for(int i=0;i<NElems();i++)
    if(fabs(v.at(i))<eps) v.at(i)=0;}

    bool EqualTo(const Vector &v2,double eps=-1) const
    {
        if(this->v.size()!=v2.NElems()) return false;
        if(eps<0) eps=this->GetEpsilon();
        for(int i=0; i<v.size(); i++)
            if(fabs(v.at(i)-v2[i])>=eps) return false;

        return true;
    }
};
double VectorNorm(const Vector &v){
    return v.Norm();
}
void PrintVector(const Vector &v, char separator,double eps){
    v.Print(separator,eps);
}
Vector operator +(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i=1; i<=v3.NElems(); i++)
        v3(i)=v1(i)+v2(i);
    return v3;
}
Vector operator -(const Vector &v1, const Vector &v2){
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i=1; i<=v3.NElems(); i++)
        v3(i)=v1(i)-v2(i);
    return v3;
}
Vector operator *(const Vector &v, double s){
    Vector v1(v.NElems());
    for(int i=1; i<=v.NElems(); i++)
        v1(i)=v(i)*s;
    return v1;
}
Vector operator *(double s, const Vector &v){
    return v*s;
}
double operator *(const Vector &v1, const Vector &v2){
    double product {};
    if(v1.NElems()!=v2.NElems()) throw std::domain_error("Incompatible formats");
    for(int i=1; i<=v1.NElems(); i++)
        product+=v1(i)*v2(i);
    return product;
}
Vector operator /(const Vector &v, double s){
    if(s<std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
    Vector v1(v.NElems());
    for(int i=1; i<=v.NElems(); i++)
        v1(i)=v(i)/s;
    return v1;
}

/////////////////////
//// M A T R I X ///
////////////////////
class Matrix {
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
    for(int i=0;i<v.NElems();i++) {
        matrix[i].resize(1);
        matrix[i][0]=v[i];
    }}
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
    void Transpose() {
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

    void ZamijeniKolone(int i,int j){
        double temp;
    for(int k=0;k<NRows();k++) {
        temp=matrix[k][i];
        matrix[k][i]=matrix[k][j];
        matrix[k][j]=temp;
    }
    }
    void ZamijeniRed(int i,int j){
         double temp;
    for(int k=0;k<NCols();k++) {
        temp=matrix[i][k];
        matrix[i][k]=matrix[j][k];
        matrix[j][k]=temp;
    }
    }
    
    void Chop(double eps = -1) {
        if(eps<0) eps=GetEpsilon();
        for(int i=0; i<matrix.size(); i++)
            for(int j=0; j<matrix.at(i).size(); j++)
                if(fabs(matrix.at(i).at(j))<eps) matrix.at(i).at(j)=0;
    }
    
    bool  EqualTo(const Matrix &m, double eps=-1){
        if(this->NRows()!=m.NRows() || this->NCols()!=m.NCols())return false;
        if(eps<0) eps=GetEpsilon();
        for(int i=0; i<matrix.size(); i++)
            for(int j=0; j<matrix.at(i).size(); j++)
                if(fabs( matrix.at(i).at(j)-m.matrix.at(i).at(j))>=eps) return false;

        return true;
    }

    friend Matrix LeftDiv(Matrix m1, Matrix m2);
    friend Vector LeftDiv(Matrix m, Vector v);

    friend Matrix operator /(const Matrix &m, double s);
    Matrix &operator /=(double s){
        if(s<this->GetEpsilon()) throw std::domain_error("Division by zero");
        for(vector<double> &v:this->matrix)
            for(double &v : v)
                v/=s;
        return *this;
    }
    
    friend Matrix operator /(Matrix m1, Matrix m2);
    Matrix &operator /=(Matrix m){  
        if(m.NRows()!=m.NCols()) throw std::domain_error("Divisor matrix is not square");
        if(NCols()!=m.NRows()) throw std::domain_error("Incompatible formats");

        for(int k=0; k<m.NCols(); k++) {
            int p=k;
            for(int i=k+1; i<m.NCols(); i++) {
                if (fabs(m[k][i])>fabs(m[k][p])) p=i;}

                if(fabs(m[k][p])<m.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");
                if(p!=k) {
                    this->ZamijeniKolone(k,p);
                    m.ZamijeniKolone(k,p);
                }
                for(int i=k+1; i<m.NCols(); i++) {
                    double u = m[k][i]/ m[k][k];
                    for(int j=k+1; j<m.NRows(); j++) m[j][i]-=u* m[j][k];
                    for(int j=0; j<this->NRows(); j++) (*this)[j][i]-=u*(*this)[j][k];
                }
            }

        Matrix x(this->NRows(),m.NCols());
        for(int k=0; k<this->NRows(); k++) {
            for(int i=m.NCols()-1; i>=0; i--) {
                double s = (*this)[k][i];
                for(int j=i+1; j<m.NCols(); j++) s=s-m[j][i]*x[k][j];
                x[k][i]= s / m[i][i];
            }
        }

        *this=x;
        return *this;
    }

    friend double Det(Matrix m);
    double Det() const;

    void Invert(){
    if(this->NRows()!=this->NCols()) throw std::domain_error("Matrix is not square");  
    vector<int> w(this->NRows());
    for(int k=1; k<=this->NRows(); k++) {
        int p=k;
        for(int i=k+1; i<=this->NRows(); i++) {
            if (fabs((*this)(i,k))>fabs((*this)(p,k))) p=i;}

            if(fabs((*this)(p,k))<std::numeric_limits<double>::epsilon()*this->NRows()) throw std::domain_error("Matrix is singular");
            if(p!=k) {
                this->matrix.at(p-1).swap(this->matrix.at(k-1));
            }
            
            w.at(k-1)=p;
            double u=(*this)(k,k);
            (*this)(k,k)=1;
            for(int j=1;j<=this->NRows();j++) (*this)(k,j)/=u;
            for(int i=1;i<=this->NRows();i++){
                if(i!=k){
                    u=(*this)(i,k);
                    (*this)(i,k)=0;
                    for(int j=1;j<=this->NRows();j++) 
                    (*this)(i,j)-=u*(*this)(k,j);
                }
            }
    }
     
     for(int j=this->NRows();j>=1;j--){
         int p=w.at(j-1);
         if(p!=j-1)
         for(int i=1;i<=this->NRows();i++)
         std::swap((*this)(i,p),(*this)(i,j));
     }
    }
    friend Matrix Inverse(Matrix m);

    void ReduceToRREF(){
        int k=-1, l=-1, p;
    vector<bool> w(this->NCols(),false);
    double eps=GetEpsilon();
    double mi;
    
    while(k<NRows() && l<NCols()) {
        l++;
        k++;
        double v(0);
        while(v<eps && l<NCols()) {
            p=k;
            for(int i(k);i<NRows();i++) {
                if(std::fabs(matrix[i][l])>v){
                    v=std::fabs(matrix[i][l]);
                    p=i;}
                }
            if(v<eps) l++;
        }
        if(l<NCols()) {
            w[l]=true;
            if(p!=k)
                this->matrix.at(k).swap(this->matrix.at(p));
            mi=matrix[k][l];
            for(int j(l);j<NCols();j++)
                matrix[k][j]/=mi;
            for(int i=0;i<NRows();i++)
                if(i!=k) {
                    mi=matrix[i][l];
                    for(int j(l);j<NCols();j++) {
                        matrix[i][j]-=mi*matrix[k][j];
                    }
                }
        }
    }
    }
    friend Matrix RREF(Matrix m);

    int Rank() const {
        Matrix m(*this);
        m.ReduceToRREF();
        int rang{};
        double eps=this->GetEpsilon();
        for(int i=0;i<this->NRows();i++)
        for(int j=0;j<this->NCols();j++)
        if(fabs(m.matrix.at(i).at(j))>eps){rang++; break;}
        
        return rang;
    }
    friend int Rank(Matrix m);
};
void PrintMatrix(const Matrix &m, int width, double eps) {
    m.Print(width,eps);
}
Matrix operator +(const Matrix &m1, const Matrix &m2) {
    if(m1.NRows()!=m2.NRows() || m1.NCols()!=m2.NCols()) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(),m1.NCols());
    for(int i=0; i<m3.NRows(); i++)
        for(int j=0; j<m3.NCols(); j++)
            m3[i][j]=m1[i][j]+m2[i][j];
    return m3;
}
Matrix operator -(const Matrix &m1, const Matrix &m2) {
    if(m1.NRows()!=m2.NRows() || m1.NCols()!=m2.NCols()) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(),m1.NCols());
    for(int i=0; i<m3.NRows(); i++)
        for(int j=0; j<m3.NCols(); j++)
            m3[i][j]=m1[i][j]-m2[i][j];
    return m3;
}
Matrix operator *(const Matrix &m, double s) {
    Matrix m1(m.NRows(),m.NCols());
    for(int i=0; i<m1.NRows(); i++)
        for(int j=0; j<m1.NCols(); j++)
            m1[i][j]=m[i][j]*s;
    return m1;
}
Matrix operator *(double s, const Matrix &m) {
    return m*s;
}
Matrix operator *(const Matrix &m1, const Matrix &m2) {
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
Vector operator *(const Matrix &m, const Vector &v) {
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
int Rank(Matrix m){
return m.Rank();    
}
Matrix RREF(Matrix m){
 m.ReduceToRREF();
 return m;
}
Matrix Inverse(Matrix m){
    m.Invert();
    return m;
}
Matrix operator /(const Matrix &m, double s) {
    if(s<m.GetEpsilon()) throw std::domain_error("Division by zero");
    Matrix m2(m);
    return m2/=s;
}
Matrix operator /(Matrix m1, Matrix m2) {
 return m1/=m2;
}
double Matrix::Det() const {
    if(this->NRows()!=this->NCols()) throw std::domain_error("Matrix is not square");
    Matrix m2=*this;
    double v=1;
    for(int k=1; k<=m2.NRows(); k++) {
        int p=k;
        for(int i=k+1; i<=m2.NRows(); i++) {
            if (fabs(m2(i,k))>fabs(m2(p,k))) p=i;}

            if(fabs(m2(p,k))<m2.GetEpsilon()) return 0;
            if(p!=k) {
                m2.matrix.at(p-1).swap(m2.matrix.at(k-1));
                v=-v;
            }
            v*=m2(k,k);
            for(int i=k+1; i<=m2.NRows(); i++) {
                double u = m2(i,k)/m2(k,k);
                for(int j=k+1; j<=m2.NRows(); j++) m2(i,j)-=u*m2(k,j);
            }
        }
    return v;
}
double Det(Matrix m) {
    return m.Det();
}
Matrix LeftDiv(Matrix m1, Matrix m2) {
    if(m2.NRows()!=m1.NRows()) throw std::domain_error("Incompatible formats");
    if(m1.NRows()!=m1.NCols()) throw std::domain_error("Divisor matrix is not square");
    
    for(int k=1; k<=m1.NRows(); k++) {
        int p=k;
        for(int i=k+1; i<=m1.NRows(); i++){
            if (fabs(m1(i,k))>fabs(m1(p,k))) p=i;}

            if(fabs(m1(p,k))<epsilon) throw std::domain_error("Divisor matrix is singular");
            if(p!=k) {
                m1.matrix.at(p-1).swap(m1.matrix.at(k-1));
                m2.matrix.at(p-1).swap(m2.matrix.at(k-1));
            }
            for(int i=k+1; i<=m1.NRows(); i++) {
                double u = m1(i,k)/m1(k,k);
                for(int j=k+1; j<=m1.NRows(); j++) m1(i,j)-=u*m1(k,j);
                for(int j=1; j<=m2.NCols(); j++) m2(i,j)-=u*m2(k,j);
            }
    }

    Matrix x(m1.NRows(),m2.NCols());
    for(int k=1; k<=x.NCols(); k++) {
        for(int i=x.NRows(); i>=1; i--) {
            double s = m2(i,k);
            for(int j=i+1; j<=x.NRows(); j++) s=s-m1(i,j)*x(j,k);

            x(i,k)= s / m1(i,i);
        }
    }

    return x;
}
Vector LeftDiv(Matrix m1, Vector v){
    Matrix m2(v);
    Matrix x = LeftDiv(m1,m2);
    Vector v1(x.NRows());
    for(int i=0;i<x.NRows();i++)
    v1[i]=x.matrix[i][0];
    
    return v1;
}
Matrix Transpose(const Matrix &m) {
    Matrix m1=m;
    m1.Transpose();
    return m1;
}
double MatrixNorm(const Matrix &m) {
    return m.Norm();
}

//////////////////////
//// LUDECOMPOSER ///
////////////////////
class LUDecomposer {
        Matrix matrix;
        Vector zamijenjeniRedovi;
    public:
        LUDecomposer(Matrix m);
        void Solve(const Vector &b, Vector &x) const;
        Vector Solve(Vector b) const;
        void Solve(Matrix &b, Matrix &x) const;
        Matrix Solve(Matrix b) const;
        Matrix GetCompactLU() const { return matrix; }
        Matrix GetL() const;
        Matrix GetU() const;
        Vector GetPermuation() const { return zamijenjeniRedovi; }
};
LUDecomposer::LUDecomposer(Matrix m) : matrix(m), zamijenjeniRedovi(m.NRows()) {
    if(m.NRows()!=m.NCols()) throw std::domain_error("Matrix is not square");
    if(std::fabs(Det(m))<epsilon) throw std::domain_error("Matrix is singular");
    
    for(int j=0;j<m.NRows();j++) {
        double s;
        for(int i=0;i<=j;i++) {
            s=matrix[i][j];
            for(int k=0;k<i;k++)
                s-=matrix[i][k]*matrix[k][j];
            matrix[i][j]=s;
        }
        int p=j;
        for(int i=j+1;i<m.NRows();i++) {
            s=matrix[i][j];
            for(int k=0;k<j;k++) s= s - matrix[i][k]*matrix[k][j];
            matrix[i][j]=s;
            if(fabs(s)>fabs(matrix[p][j])) p=i;
        }
        if(p!=j) matrix.ZamijeniRed(p,j);
        zamijenjeniRedovi[j]=p;
        double mi=matrix[j][j];
        for(int i=j+1;i<m.NRows();i++)
            matrix[i][j]/=mi;
    }
}
Matrix LUDecomposer::GetL() const {
    Matrix L(matrix.NRows(),matrix.NCols());
    for(int i=0;i<matrix.NRows();i++) {
        for(int j=0;j<matrix.NCols();j++) {
            if(i==j) L[i][j]=1;
            else if(i>j) L[i][j]=matrix[i][j];
            else L[i][j]=0;
        }
    }
    return L;
}
Matrix LUDecomposer::GetU() const {
    Matrix U(matrix.NRows(),matrix.NCols());
    for(int i=0;i<matrix.NRows();i++) {
        for(int j=0;j<matrix.NCols();j++) {
            if(j>=i) U[i][j]=matrix[i][j];
        }
    }
    return U;
}
void LUDecomposer::Solve(Matrix &B, Matrix &X) const {
    if(B.NCols()!=matrix.NRows() || B.NRows()!=matrix.NCols() 
    || X.NRows()!=matrix.NCols() || X.NCols()!=matrix.NRows()) throw std::domain_error("Incompatible formats");
    X=B;
    // Supstitucija unaprijed 
    for(int k=0;k<X.NRows();k++) {
        for(int i=0;i<matrix.NRows();i++) {
            int p=zamijenjeniRedovi[i];
            double s=X[p][k];
            X[p][k]=X[i][k];
            for(int j=0;j<i;j++) s-=matrix[i][j]*X[j][k];
            X[i][k]=s;
        }
        
        //Supstitucija unazad
         for(int i=matrix.NCols()-1;i>=0;i--) {
            double s=X[i][k];
            for(int j=i+1;j<matrix.NCols();j++) {
                s-=matrix[i][j]*X[j][k];
            }
            X[i][k]=s/matrix[i][i];
        }
    }
} 
Matrix LUDecomposer::Solve(Matrix b) const {
    Matrix x(b.NRows(),b.NCols());
    this->Solve(b,x);
    return x;
}
void LUDecomposer::Solve(const Vector &b, Vector &x) const {
    if(x.NElems()!=matrix.NRows() || b.NElems()!=matrix.NRows()) throw std::domain_error("Incompatible formats");
    x=b;
    for(int i=0;i<x.NElems();i++) {
        int p=zamijenjeniRedovi[i];
        double s=x[p];
        x[p]=x[i];
        for(int j=0;j<i;j++)
            s-=matrix[i][j]*x[j];
        x[i]=s;
    }
    
    for(int i=b.NElems()-1;i>=0;i--) {
        double s=x[i];
        for(int j=i+1;j<b.NElems();j++) s-=matrix[i][j]*x[j];
            x[i]=s/matrix[i][i];
    }
}
Vector LUDecomposer::Solve(Vector b) const {
    Vector x(b.NElems());
    this->Solve(b,x);
    return x;
}


//////////////////////
//// QRDECOMPOSER ///
////////////////////
class QRDecomposer {
    private:
        Matrix matrix;
        Vector v;
    protected: 
        Vector MulQWithVector(Vector b, bool T) const;
        Matrix MulQWithMatrix(Matrix b, bool T) const;
    public:
        QRDecomposer(Matrix m);
        
        Matrix GetQ() const;
        Matrix GetR() const;
        
        Vector MulQWith(Vector v) const {return MulQWithVector(v,false);}
        Vector MulQTWith(Vector v) const {return MulQWithVector(v,true);}
        Matrix MulQWith(Matrix m) const {return MulQWithMatrix(m,false);}
        Matrix MulQTWith(Matrix m) const {return MulQWithMatrix(m,true);}
        
        void Solve(const Vector &b, Vector &x) const;
        Vector Solve(Vector b) const;
        void Solve(Matrix &b, Matrix &x) const;
        Matrix Solve(Matrix b) const;
};
QRDecomposer::QRDecomposer(Matrix m) : matrix(m), v(m.NCols()) {
    if(m.NRows()<m.NCols()) throw std::domain_error("Invalid matrix format");
    for(int k=0;k<matrix.NCols();k++) {
        double s=0;
        for(int i=k;i<matrix.NRows();i++) s+=matrix[i][k]*matrix[i][k];
        s=std::sqrt(s);
        double u=std::sqrt(s*(s+std::fabs(matrix[k][k])));
        if(std::fabs(u)<epsilon) throw std::domain_error("Matrix is singular");
        if(matrix[k][k]<0) s=-s;
        matrix[k][k]=(matrix[k][k]+s)/u;
        for(int i=k+1;i<matrix.NRows();i++) matrix[i][k]=matrix[i][k]/u;
        v[k]=-s;
        for(int j=k+1;j<matrix.NCols();j++) {
            s=0;
            for(int i=k;i<matrix.NRows();i++) s+=matrix[i][k]*matrix[i][j];
            for(int i=k;i<matrix.NRows();i++) matrix[i][j]-=s*matrix[i][k];
        }
    }
}
Matrix QRDecomposer::GetQ() const {
    Matrix Q(matrix.NRows(), matrix.NRows());
    for(int j=0;j<matrix.NRows();j++) {
        for(int i=0;i<matrix.NRows();i++) Q[i][j]=0;
        Q[j][j]=1;
        for(int k=matrix.NCols()-1;k>=0;k--) {
            double s=0;
            for(int i=k;i<matrix.NRows();i++) s+=matrix[i][k]*Q[i][j];
            for(int i(k);i<matrix.NRows();i++) Q[i][j]-=s*matrix[i][k];
        }
    }
    return Q;
}
Matrix QRDecomposer::GetR() const {
    Matrix R(matrix.NRows(),matrix.NCols());
    for(int i=0;i<R.NRows();i++)
        for(int j=0;j<R.NCols();j++){
            if(j<i) R[i][j]=0;
            else if (j>i) R[i][j]=matrix[i][j];
            else R[i][j] = v[i];}
    return R;
}
Matrix QRDecomposer::MulQWithMatrix(Matrix B, bool T) const {
    if(matrix.NRows()!=B.NRows()) throw std::domain_error("Incompatible formats");
    Vector v(matrix.NRows());
    for(int i=0;i<B.NCols();i++) {
        for(int j=0;j<v.NElems();j++) v[j]=B[j][i];
        v=MulQWithVector(v,T);
        for(int j=0;j<v.NElems();j++) B[j][i]=v[j];
    }
    return B;
}
Vector QRDecomposer::MulQWithVector(Vector b, bool T) const {
    if(matrix.NRows()!=b.NElems()) throw std::domain_error("Incompatible formats");
    
    double s{};
    if(!T) {
        for(int k=matrix.NCols()-1;k>=0;k--) {
            s=0;
            for(int i=k;i<matrix.NRows();i++) s+=matrix[i][k]*b[i];
            for(int i=k;i<matrix.NRows();i++) b[i]-=s*matrix[i][k];
        }
        for(int i=0;i<matrix.NRows();i++)
            if(std::fabs(b[i])<epsilon) b[i]=0;
    }
    
    else {
        for(int k=0;k<matrix.NCols();k++) {
            s=0;
            for(int i(k);i<matrix.NRows();i++) s+=matrix[i][k]*b[i];
            for(int i(k);i<matrix.NRows();i++) b[i]-=s*matrix[i][k];
        }
        for(int i=0;i<matrix.NRows();i++)
            if(std::fabs(b[i]) < epsilon) b[i]=0;
    }
    return b;
}
void QRDecomposer::Solve(const Vector &b, Vector &x) const {
    if(matrix.NRows()!=matrix.NCols()) throw std::domain_error("Matrix is not square");
    if(b.NElems()!=matrix.NRows() || b.NElems()!=x.NElems()) throw std::domain_error("Incompatible formats");
    Vector b1(MulQTWith(b));
    for(int i=matrix.NCols()-1;i>=0;i--) {
        double s=b1[i];
        for(int j=i+1;j<matrix.NCols();j++){
            if(j<i) continue;
            else if (j>i) s-=matrix[i][j]*x[j];
            else s-=v[i]*x[j];
        }
        
        x[i]=s/v[i];
    }
}
Vector QRDecomposer::Solve(Vector b) const {
    Vector v(b.NElems());
    Solve(b,v);
    return v;
}
void QRDecomposer::Solve(Matrix &B, Matrix &X) const {
    if(matrix.NRows()!=matrix.NCols()) throw std::domain_error("Matrix is not square");
    if(B.NRows()!=matrix.NCols() || B.NCols()!=matrix.NRows() 
    || X.NRows()!=matrix.NCols() || X.NCols()!=matrix.NRows()) throw std::domain_error("Incompatible formats");
    for(int i=0;i<B.NCols();i++) {
        Vector v1(B.NRows()), v2(B.NRows());
        for(int j=0;j<B.NRows();j++) v1[j]=B[j][i];
        Solve(v1,v2);
        for(int j=0;j<B.NRows();j++) X[j][i]=v2[j];
    }
}
Matrix QRDecomposer::Solve(Matrix B) const {
    Matrix M(B.NRows(), B.NCols());
    Solve(B,M);
    return M;
}

////////////////////////
//// TEST FUNCTIONS ///
//////////////////////

// TEST VECTOR
void classVectorChopEqualTest() {
    Vector v1{1,2,3,4,5,6};
    Vector v2(5);
    for(int i=0; i<5; i++) v2[i]=i*log(i+1);

    Vector v3{0.132,-1.2345,3.14159265,991.113,0.0003};
    Vector v4{11./8, 33./241,33.11,42.42,0.01, 874.423};
    
    Vector v5{1,2,3,4,5,6};
    Vector v6{epsilon,10e-21,112e-44};
    if(!v3.EqualTo(v4) && v1.EqualTo(v5)) cout << "Vector EqualTo test: Good" << endl;
    else cout << "Vector EqualTo test: Error" << endl;
    
    
    v6.Chop();
    v5.Chop(5);
    if (std::fabs(v5[1])>0) std::cout << "Vector Chop test: Error" << std::endl;
    else std::cout << "Vector Chop test: Good" << std::endl;
}

// TEST MATRIX
void classMatrixChopEqualTest(){
    {
    Matrix m1({{33, 0.4,23}, {11, 34, 2.33}, {0.34123, -133, 1337}});
    Matrix m2({{33, -4.5,24}, {11.9, 35, 3.33}, {1.34123, -134, 1340}});

    if (m1.EqualTo(m2, 5)) std::cout << "Matrix EqualTo test: Good" << std::endl;
    else std::cout << "Matrix EqualTo test: Error" << std::endl;
    }
    
    {
    Matrix m1({{33, 0.4,23}, {11, 34, 2.33}, {0.34123, -133, 1337}});
    m1.Chop(1.1);
    if (std::fabs(m1[2][0])>0) std::cout << "Matrix Chop test: Error" << std::endl;
    else std::cout << "Matrix Chop test: Good" << std::endl;
    }
        
}
void classMatrixLeftDivTest() {
    Matrix m1({{33, 0.4,23}, {11, 34, 2.33}, {0.34123, -133, 1337}});
    Matrix m2({{1}, {2}, {3}});
    Matrix m3(LeftDiv(m1, m2));
    
    Matrix m4({{0.024642},{0.0503545},{0.00724663}});
    if(!m3.EqualTo(m4,0.001)) cout << "Matrix LeftDiv test: Error"<< endl;
    else cout << "Matrix LeftDiv test: Good"<< endl;
    
    try {
        Matrix m4({{1231,3213, 12}, {34234, 987298, 0.32423}});
        LeftDiv(m1,m4);
        cout << "Matrix LeftDiv test (Incompatible) : Error"<< endl;
    }
    catch(std::domain_error) { cout << "Matrix LeftDiv test (Incompatible): Good"<< endl;}
    
    try {
        Matrix m4({{1231,3213, 12}, {34234, 987298, 0.32423}});
        Matrix m5({{0,1}, {3, 4}});
        LeftDiv(m4,m5);
    cout << "Matrix LeftDiv test (Divisor matrix not square): Error"<< endl;
    }
    catch(std::domain_error) { cout << "Matrix LeftDiv test (Divisor matrix not square): Good"<< endl;}
    
    try {
        Matrix m4({{1, 0, 0}, {-2, 0, 0},{4, 6, 1}});
        Matrix m5({{1, 1}, {2, 2}, {3, 3}});
        LeftDiv(m4,m5);
        cout << "Matrix LeftDiv test (Divisor matrix singular): Error"<< endl;
    }
    catch(std::domain_error) { cout << "Matrix LeftDiv test (Divisor matrix singular): Good" << endl;}
    
    Matrix m11({{33, 0.4,23}, {11, 34, 2.33}, {0.34123, -133, 1337}});
    Vector v({51, 53,35});
    Vector rez(LeftDiv(m11, v));
    Vector v2({{1.4391657},{1.0840523},{0.1336484}});
    if(!rez.EqualTo(v2,0.0001)) cout << "Matrix LeftDiv with Vector test: Error"<< endl;
    else cout <<"Matrix LeftDiv with Vector test: Good"<< endl;
    
}
void classMatrixRightDivTest() {
    Matrix m{{133, 51, 44}, {982, 112, 132}};
    Matrix m1(m/11);
    Matrix rez({{12.090909,4.6363636,4},{89.272727,10.181818,12}});
    if (m1.EqualTo(rez,0.0001)) cout << "Matrix division by double: Good" <<endl;
    else cout << "Matrix division by double: Error"<<endl;
    

    //Division by zero:
    try {
        m1/0;
        cout << "Matrix division by zero: Error"<<endl;
    }
    catch(std::domain_error) {cout << "Matrix division by zero: Good" <<endl;}
    
    
    Matrix m2({{1, 2, 1},{1, 1, 3},{1, 1, 1}});
    Matrix m3(m/m2);
    Matrix m4({{-82,-44.5,259.5},{-870, -425, 2277}});
    if (m4.EqualTo(m3,0.0001)) cout << "Matrix RightDiv test: Good" << endl;
    else cout << "Matrix RightDiv test: Error" << endl;
    
    
    //Incompatible format
    Matrix m5({{4, 5},{5, 4},{4, 5}});
    try {
        m5/m2;
    cout << "Matrix RightDiv test (Incompatible) : Error"<< endl;
    }
    catch(std::domain_error) { cout << "Matrix RightDiv test (Incompatible): Good"<< endl;}
    
    //Divisor matrix not square
    try {
        m2/m4;
    cout << "Matrix RightDiv test (Divisor matrix not square): Error"<< endl;
    }
    catch(std::domain_error) { cout << "Matrix RightDiv test (Divisor matrix not square): Good"<< endl;}
    
    //Divisor matrix is singular:
    try {
        Matrix m4({{1, 0, 0}, {-2, 0, 0},{4, 6, 1}});
        Matrix m5({{33, 33,0}, {22, 22,22}, {31, 31, 123}});
        m4/m5;
    cout << "Matrix RightDiv test (Divisor matrix singular): Error"<< endl;
    }
    catch(std::domain_error& e) { cout << "Matrix RightDiv test (Divisor matrix singular): Good" << endl;}
    
}
void classMatrixDetTest() {
    Matrix M({{123,0.3213,213},{321,   98,      0.3214},{3.14,   2.71,     3241}});
    if(fabs(M.Det()-38852385.672692678868)<epsilon) cout << "Matrix Determinant test: Good" <<endl;
    else cout << "Matrix Determinant test: Error" << endl;
}
void classMatrixInvertTest() {
    Matrix M({{123,0.3213,213},{321,   98,      0.3214},{3.14,   2.71,     3241}});
    Matrix rez({{0.008175,-0.0000119, -0.0005373}, {-0.0267772,   0.0102432,   0.0017588 }, {0.0000145,  -0.0000086,   0.0003076}});
    if (Inverse(M).EqualTo(rez, 0.001)) cout << "Matrix Invert test: Good" <<endl;
    else cout << "Matrix Inverttest: Error" << endl;
    
    try {
        Matrix m4({{99, 0, 0}, {-3, 0, 0},{78, 0, 1}});
        Inverse(m4);
        cout << "Matrix Invert singular test: Error" << endl;
    }
    catch(std::domain_error) {cout << "Matrix Invert singular test: Good" <<endl;}
    
    try {
        Matrix m4({{13,90}, {-3.14, 0.73},{61.11, 94.33}});
        Inverse(m4);
    cout << "Matrix Invert not square test: Error" << endl;
    }
    catch(std::domain_error) {cout << "Matrix Invert not square test: Good" <<endl;}
}
void classMatrixRREFTest() {
    Matrix M({{123,0.3213},{321,   98},{3.14,   2.71}});
    Matrix rez({{1,0}, {0,1}, {0,0}});
    if (RREF(M).EqualTo(rez)) std::cout << "Matrix RREF test: Good" << std::endl;
    else std::cout << "Matrix RREF test: Error" << std::endl;
}

// TEST LU FACTORIZATION
void classLUDecomposerTest() { 
    Matrix m({{123,0.3213,213},{321,   98,      0.3214},{3.14,   2.71,     3241}}); 
    LUDecomposer lu(m);
    Matrix rezU({{321, 98, 0.3214}, { 0,    -37.230102,  212.87685}, {0,    0,         3251.011}});
    Matrix rezL({{1, 0, 0}, {0.383178 , 1, 0}, {0.00978193,-0.0470418, 1}});
    if (lu.GetL().EqualTo(rezL, 0.001) && lu.GetU().EqualTo(rezU, 0.001)) cout << "LUDecomposer test: Good" << endl;
    else cout << "LUDecomposer test : Error" << endl;

    Matrix LiU({{321, 98, 0.3214}, { 0.383178,    -37.230102,  212.87685}, {0.00978193,-0.0470418,        3251.011}});
    if (lu.GetCompactLU().EqualTo(LiU, 0.1)) cout << "LUDecomposer CompactLU test: Good" << endl;
    else cout << "LUDecomposer CompactLU test : Error" << endl;

}
void clasLUDecomposerSolveLUTest() {
    Matrix A({{2, 4}, {3, -3}});
    LUDecomposer lu(A);
    Vector b({-2, 5});
    Vector x(2);
    lu.Solve(b, x);
    if (x.EqualTo({0.777778, -0.888889}, 0.0001)) std::cout << "LUDecomposer Solve test: Good" << std::endl;
    else std::cout << "LUDecomposer Solve test: Error" << std::endl;
    
    Matrix m1({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Matrix b1({{4, 4, 4}, {8, 7, 1}, {2, 2, 3}});
    Matrix x1(3,3);
    LUDecomposer lu1(m1);
    lu1.Solve(b1, x1);
    
    Matrix rez({{5.08108, 3.59459, -5.43243}, {0.27027, 0.648649, 2.89189}, {-0.540541, -0.297297, 1.21622}});
    if (lu1.Solve(b1).EqualTo(rez, 0.0001))  std::cout << "LUDecomposer Solve test (Matrix): Good" << std::endl;
    else std::cout << "LUDecomposer Solve test (Matrix): Error" << std::endl;
}

// TEST QR FACTORIZATION
void classQRDecomposerTest() {
    Matrix m({{123,0.3213,213},{321,   98,      0.3214},{3.14,   2.71,     3241}}); 
    QRDecomposer qr(m);
    Matrix r({{-343.773,  -91.6478,  -106.11}, {0,-34.8157,  24.1848}, {0,          0, -3246.17}});
    Matrix q({{-0.357794,0.932618,-0.0469717}, {-0.9337558,  -0.3568322, -0.0277659}, {-0.00913393,-0.0537946,  -0.99851}});
    if (qr.GetR().EqualTo(r, 0.1) && qr.GetQ().EqualTo(q, 0.1)) cout << "QRDecomposer test: Good" << endl;
    else cout << "QRDecomposer test: Error" << endl;
}
void classQRDecomposerSolveTest() {
    Matrix A({{2, 4}, {3, -3}});
    QRDecomposer qr(A);
    Vector b({-2, 5});
    Vector x(2);
    qr.Solve(b, x); 
    if (x.EqualTo({0.777778, -0.888889}, 0.0001)) std::cout << "QRDecomposer Solve test: Good" << std::endl;
    else std::cout << "QRDecomposer Solve test: Error" << std::endl;
    

    Matrix m1({{1, 2, 3}, {2, 2, 5}, {3, -3, 23}});
    Matrix b1({{4, 4, 4}, {8, 7, 1}, {2, 2, 3}});
    Matrix x1(3,3);
    QRDecomposer qr1(m1);
    qr1.Solve(b1, x1);
    
    Matrix xrez({{5.08108, 3.59459, -5.43243}, {0.27027, 0.648649, 2.89189}, {-0.540541, -0.297297, 1.21622}});
    if (qr1.Solve(b1).EqualTo(xrez, 0.0001)) std::cout << "QRDecomposer Solve test (matrix): Good" << std::endl;
    else std::cout << "QRDecomposer Solve test (matrix): Error" << std::endl;
}


int main (){
    classVectorChopEqualTest();
    
    classMatrixChopEqualTest();
    classMatrixLeftDivTest(); 
    classMatrixRightDivTest(); 
    classMatrixDetTest();
    classMatrixInvertTest(); 
    classMatrixRREFTest();
    
    classLUDecomposerTest();
    clasLUDecomposerSolveLUTest();
    
    classQRDecomposerTest();
    classQRDecomposerSolveTest();

    return 0;
}
