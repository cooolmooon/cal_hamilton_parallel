#include<fstream>
#include"spline.h"
#include"timer.h"
#include<cmath>
using namespace std;
//-----------------------------------------------------------
//EXPLAIN:the construction function
//-----------------------------------------------------------
spline::spline(const string&file){
    timer::tick("spline","interplation");
    ifstream Distribution(file);
    if(!Distribution){
        //failed to read the input
        cout<<"failed to read distribution"<<endl;
        return;
    }
    char c;
    string name;
    Distribution>>name>>cutoff
        >>name>>dr
        >>name>>mesh
        >>name>>l>>name;
    for(int i=0;i<10000;i++){
        m[i]=0;
        fvalues[i]=0;
    }
    for(int i=0;i<mesh;i++){
        Distribution>>fvalues[i]>>c;
    }
    Distribution.close();
    interplation();
    timer::tick("spline","interplation");
}
spline::~spline(){}
//-----------------------------------------------------------
//EXPLAIN:get the derivatives by chase method
//-----------------------------------------------------------
void spline::chase(double* alpha,double* beta,double* d){
    double* l, * u, * y,*c;
    l = new double[mesh];
    u = new double[mesh];
    y = new double[mesh];
    c = new double[mesh];
    c[0] = 1;
    for (int i = 1; i < mesh-1; i++) {
        c[i] = 2;
    }
    c[mesh-1] = 1;
    u[0] = 1;
    for (int i = 1; i <= mesh; i++) {
        l[i] = alpha[i] / u[i - 1];
        u[i] = c[i] - beta[i - 1] * l[i];
    }
    y[0] = d[0];
    for (int i = 1; i <= mesh-1; i++) {
        y[i] = d[i] - l[i] * y[i - 1];
    }
    m[mesh-1] = y[mesh-1] / u[mesh-1];
    for (int i = mesh - 1; i >= 0; i--) {
        m[i] = (y[i] - beta[i] * m[i + 1]) / u[i];
    }
    delete[] l,u,y,c;
}
//-----------------------------------------------------------
//EXPLAIN:called in the construction function
//-----------------------------------------------------------
void spline::interplation(){
    int i;
    double* alpha;
    double a = (fvalues[1]-fvalues[0])/dr; double b = 0.0;
    alpha = new double[mesh];
    for (i = 1; i <mesh-1; i++) {
        alpha[i] = 0.5;
    }
    alpha[mesh-1] = 0;
    double* d;
    d = new double[mesh];
    d[0] = a;
    for (i = 1; i < mesh-1; i++) {
        d[i] = 3/(2*dr) * (fvalues[i + 1] - fvalues[i - 1]);
    }
    d[mesh-1] = b;
    double* beta;
    beta = new double[mesh];
    beta[0] = 0;
    for (i = 1; i <mesh-1; i++) {
        beta[i] = 0.5;
    }
    chase(alpha, beta, d);
    delete[] alpha,beta,d;
}
//-----------------------------------------------------------
//EXPLAIN:calculate f(r) of given r
//-----------------------------------------------------------
double spline::cal(double r){
    int n=floor(r/dr);
    double p1=r/dr-n;
    double p2=n+1-r/dr;
    double ans=(1+2*p1)*p2*p2*fvalues[n]+(1+2*p2)*p1*p1*fvalues[n+1]
                +(r-n*dr)*p2*p2*m[n]+(r-(n+1)*dr)*p1*p1*m[n+1];
    return ans;
}