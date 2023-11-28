#ifndef SPLINE_H
#define SPLINE_H

#include<string>

using namespace std;

class spline{
    public:
//-----------------------------------------------------------
//EXPLAIN:parameters of distribution
//-----------------------------------------------------------
    int cutoff;
    int mesh;
    double dr;
    int l;
    double fvalues[10000];
//-----------------------------------------------------------
//EXPLAIN:the derivatives of each node
//-----------------------------------------------------------
    double m[10000];
    spline(const string& file);
    ~spline();
//-----------------------------------------------------------
//EXPLAIN:get the derivatives by chase method
//-----------------------------------------------------------
    void chase(double* alpha,double* beta,double* d);
//-----------------------------------------------------------
//EXPLAIN:get the derivatives
//-----------------------------------------------------------
    void interplation();
//-----------------------------------------------------------
//EXPLAIN:calculate the value of given r
//-----------------------------------------------------------
    double cal(double r);
    
};

#endif