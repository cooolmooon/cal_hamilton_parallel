#ifndef INPUT_H
#define INPUT_H

#include<string>
#include<iostream>

using namespace std;

class input{
    public:
    //global parameters
    int isHexahedral;
    int lx;
    int ly;
    int lz;
    int thetaxy;
    int thetayz;
    int thetaxz;
    int support_SH;
    string diago_lib;
    int support_Periodic_Boundary;
    int multi_parallel_strategies;
    string points_path;
    string v_path;
    string distribution_path;
    //points
    int npoints;
    double pos_x[50];
    double pos_y[50];
    double pos_z[50];
    //V
    int nx;
    int ny;
    int nz;
    double* v_value;

    input(const string& filename);
    ~input();

    friend ostream& operator<<(ostream&os,const input& in);

    private:
    bool read_file(const string& filename);
};
#endif