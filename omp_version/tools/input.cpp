#include"input.h"
#include<iostream>
#include<fstream>
#include<string>
#include<algorithm>
#include<sstream>
#include"timer.h"

using namespace std;

input::input(const string&filename){
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

input::~input(){
    delete[] v_value;
};

bool input::read_file(const string&filename){
    timer::tick("input","read_file");
    ifstream infile(filename);
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
    char c;
//-----------------------------------------------------------
//EXPLAIN:read_points
//-----------------------------------------------------------
    ifstream Point(points_path);
    npoints=0;
        if (!Point) {
        //failed to read the input
        cout<<"failed to read point"<<endl;
        return false;
    }
    string line;
    while(getline(Point,line)){
        if (line.empty() || std::all_of(line.begin(), line.end(), ::isspace)) {
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
    ifstream V(v_path);
    cout<<"reading v..."<<endl;
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
    cout<<"finish reading V!"<<endl;
    V.close();
    timer::tick("input","read_file");
    return true;
}

ostream& operator<<(ostream&os,const input& in){
    os<<"distribution:"<<in.distribution_path<<endl
        <<"npoints:"<<in.npoints<<endl
        <<"point 2:"<<in.pos_x[1]<<in.pos_y[1]<<in.pos_z[1]<<endl
        <<"V 0 1 4:"<<in.v_value[53]<<endl;
    return os;
}