#ifndef INPUT_H
#define INPUT_H


#include<iostream>
#include<fstream>
#include<string>
#include<algorithm>
#include<sstream>

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
    int pos_x[50];
    int pos_y[50];
    int pos_z[50];
    //V
    int nx;
    int ny;
    int nz;
    double* v_value;

    input(const string& filename){
    isHexahedral=0;
    lx=1000;
    ly=1000;
    lz=1000;
    thetaxy=0;
    thetayz=0;
    thetaxz=0;
    support_SH=0;
    diago_lib="lapack";
    support_Periodic_Boundary=0;
    multi_parallel_strategies=0;
    points_path="";
    v_path="";
    distribution_path="";
    v_value=NULL;
    //read the input
    if(!input::read_file(filename)){
        cout << "Failed to read file " << filename << endl;
    }
    }
    ~input(){
        delete[] v_value;
    }

    friend ostream& operator<<(ostream&os,const input& in){
        os<<"distribution:"<<in.distribution_path<<endl
            <<"npoints:"<<in.npoints<<endl
            <<"point 2:"<<in.pos_x[1]<<in.pos_y[1]<<in.pos_z[1]<<endl
            <<"V 0 1 4:"<<in.v_value[53]<<endl;
        return os;
    }

    private:
    bool read_file(const string& filename){
        cout<<"reading input..."<<endl;
        ifstream infile(filename.c_str());
        if (!infile) {
            //failed to read the input
            cout<<"failed to read input"<<endl;
            return false;
        }
        string name;
        infile>>name>>isHexahedral
            >>name>>lx
            >>name>>ly
            >>name>>lz
            >>name>>thetaxy
            >>name>>thetayz
            >>name>>thetaxz
            >>name>>support_SH
            >>name>>diago_lib
            >>name>>support_Periodic_Boundary
            >>name>>multi_parallel_strategies
            >>name>>points_path
            >>name>>v_path
            >>name>>distribution_path;
        infile.close();
//-----------------------------------------------------------
//EXPLAIN:read_points
//-----------------------------------------------------------
        ifstream Point(points_path.c_str());
        npoints=0;
            if (!Point) {
            //failed to read the input
            cout<<"failed to read point"<<endl;
            return false;
        }
        string line;
        char c;
        while(getline(Point,line)){
            if (line.empty()) {
                continue;
            }
            istringstream in_line(line);
            in_line>>c>>pos_x[npoints]
                >>c>>pos_y[npoints]
                >>c>>pos_z[npoints]>>c;
            npoints+=1;
        }
        Point.close();
//-----------------------------------------------------------
//EXPLAIN:read_V
//-----------------------------------------------------------
        ifstream V(v_path.c_str());
        printf("reading V...\n");
        if (!V) {
            //failed to read the input
            cout<<"failed to read V"<<endl;
            return false;
        }
        V>>name>>nx>>name>>ny>>name>>nz>>name;
        v_value=new double[nx*ny*nz];
        for(int i=0;i<nx*ny*nz;i++){
            V>>v_value[i];
        }
        printf("finish reading V!\n");
        V.close();
        return true;
    }
};
#endif