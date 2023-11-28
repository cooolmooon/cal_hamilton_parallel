#include"diago.h"
#include<lapacke.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<string>

using namespace std;

diago::diago(){}
diago::~diago(){}

void diago::lapack_diago(int n,double* mat){
    //save the eigenvalue 
    double ev[n];
    //initialize
    int info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',n,mat,n,ev);
    //chech whether the diagonazation succeed
    if(info!=0){
        cout<<"Failed to diagonalize matrix!"<<endl;
    }
//-----------------------------------------------------------
//EXPLAIN:creater folder to save the output files
//-----------------------------------------------------------
    string fn="output";
    std::stringstream ss;
    ss<<"test -d "<<fn<<" || mkdir "<<fn;
    system(ss.str().c_str());
    //output the result to file
    ofstream valuelog("./output/eigenvalues.log");
    valuelog<<"Eigenvalues are: "<<endl;
    for(int i=0;i<n;i++){
        valuelog<<ev[i]<<" ";
    }
    valuelog.close();
    ofstream vectorlog("./output/eigenvectors.log");
    vectorlog<<"Eigenvectors are: "<<endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            vectorlog<<setw(15)<<mat[n*j+i];
        }
        vectorlog<<endl;
    }
    vectorlog.close();
    return;
}
