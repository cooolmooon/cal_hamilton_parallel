#include"input.h"
#include"timer.h"
#include"spline.h"
#include"diago.h"
#include<omp.h>
using namespace std;
#include<sstream>
#include<cmath>
#include<iomanip>
int main(int argc,char **argv){
//-----------------------------------------------------------
//EXPLAIN:initialize
//-----------------------------------------------------------

    timer clock;
    clock.start();
    const int nthreads=8;
    omp_set_num_threads(nthreads);
//-----------------------------------------------------------
//EXPLAIN:read the input
//-----------------------------------------------------------
    input in("./input/INPUT.txt");
    spline Spline(in.distribution_path);
//-----------------------------------------------------------
//EXPLAIN:preprocessing
//-----------------------------------------------------------
    double rr=Spline.cutoff*Spline.cutoff;
    bool* flag=new bool[in.npoints*in.npoints];
    for(int i=0;i<in.npoints*in.npoints;i++){
        flag[i]=0;
    }
    for(int i=0;i<in.npoints;i++){
        for(int j=i+1;j<in.npoints;j++){
            if(((in.pos_x[i]-in.pos_x[j])*(in.pos_x[i]-in.pos_x[j])
            +(in.pos_y[i]-in.pos_y[j])*(in.pos_y[i]-in.pos_y[j])
            +(in.pos_z[i]-in.pos_z[j])*(in.pos_z[i]-in.pos_z[j]))
            >4*rr){
                flag[in.npoints*i+j]=1;
            }
        }
    }
    // ofstream fl("flag.txt");
    // fl<<"dimension: "<<in.npoints<<endl;
    // for(int i=0;i<in.npoints;i++){
    //     for(int j=0;j<in.npoints;j++){
    //         fl<<setprecision(10)<<setw(20)<<flag[in.npoints*i+j];
    //     }
    //     fl<<endl;
    // }
    // fl.close();
//-----------------------------------------------------------
//EXPLAIN:initialize the H matrix
//-----------------------------------------------------------
    double* matH;
    matH=new double[in.npoints*in.npoints];
    for(int i=0;i<in.npoints*in.npoints;i++){
        matH[i]=0;
    }
    double tx=(double)in.lx/(in.nx);
    double ty=(double)in.ly/(in.ny);
    double tz=(double)in.lz/(in.nz);
    double dv=tx*ty*tz;

    timer::tick("main","integral");
//-----------------------------------------------------------
//EXPLAIN:multiple points
//-----------------------------------------------------------
    if(in.npoints>=10){
        #pragma omp parallel
        {
            #pragma omp for reduction(+:matH[:in.npoints*in.npoints]) schedule(dynamic)
            for(int pointA=0;pointA<in.npoints;pointA++){
                for(int pointB=pointA;pointB<in.npoints;pointB++){
                    if(flag[in.npoints*pointA+pointB])continue;
                    for(int i=0;i<in.nx;i++){
                        double xA=(in.pos_x[pointA]-i*tx)*(in.pos_x[pointA]-i*tx);
                        double xB=(in.pos_x[pointB]-i*tx)*(in.pos_x[pointB]-i*tx);
                        if(xA>rr||xB>rr)continue;
                        for(int j=0;j<in.ny;j++){
                            double yA=(in.pos_y[pointA]-j*ty)*(in.pos_y[pointA]-j*ty);
                            double yB=(in.pos_y[pointB]-j*ty)*(in.pos_y[pointB]-j*ty);
                            if((xA+yA)>rr||(xB+yB)>rr)continue;
                            for(int k=0;k<in.nz;k++){
                                double zA=(in.pos_z[pointA]-k*tz)*(in.pos_z[pointA]-k*tz);
                                double zB=(in.pos_z[pointB]-k*tz)*(in.pos_z[pointB]-k*tz);
                                if((xA+yA+zA)>rr||(xB+yB+zB)>rr)continue;
                                else{
                                    double fA,fB;
                                    fA=Spline.cal(sqrt(xA+yA+zA));
                                    if(pointA==pointB){
                                        fB=fA;
                                        matH[in.npoints*pointA+pointB]+=fA*fB*in.v_value[i*in.ny*in.nz+j*in.nz+k]*dv;
                                    }
                                    else{
                                        fB=Spline.cal(sqrt(xB+yB+zB));
                                        double tmp=fA*fB*in.v_value[i*in.ny*in.nz+j*in.nz+k]*dv;
                                        matH[in.npoints*pointA+pointB]+=tmp;
                                        matH[in.npoints*pointB+pointA]+=tmp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
// -----------------------------------------------------------
// EXPLAIN:npoints<10
// -----------------------------------------------------------
    else{
        for(int pointA=0;pointA<in.npoints;pointA++){
            for(int pointB=pointA;pointB<in.npoints;pointB++){
                if(flag[in.npoints*pointA+pointB])continue;
                #pragma omp parallel
                {
                    #pragma omp for reduction(+:matH[:in.npoints*in.npoints]) schedule(dynamic)
                    for(int i=0;i<in.nx;i++){
                        double xA=(in.pos_x[pointA]-i*tx)*(in.pos_x[pointA]-i*tx);
                        double xB=(in.pos_x[pointB]-i*tx)*(in.pos_x[pointB]-i*tx);
                        if(xA>rr||xB>rr)continue;
                        for(int j=0;j<in.ny;j++){
                            double yA=(in.pos_y[pointA]-j*ty)*(in.pos_y[pointA]-j*ty);
                            double yB=(in.pos_y[pointB]-j*ty)*(in.pos_y[pointB]-j*ty);
                            if((xA+yA)>rr||(xB+yB)>rr)continue;
                            for(int k=0;k<in.nz;k++){
//-----------------------------------------------------------
//EXPLAIN:rA and rB represent square r from the node to pointA(B)
//-----------------------------------------------------------
                                double zA=(in.pos_z[pointA]-k*tz)*(in.pos_z[pointA]-k*tz);
                                double zB=(in.pos_z[pointB]-k*tz)*(in.pos_z[pointB]-k*tz);
                                if((xA+yA+zA)>rr||(xB+yB+zB)>rr)continue;
                                double fA,fB;
                                fA=Spline.cal(sqrt(xA+yA+zA));
                                if(pointA==pointB){
                                    fB=fA;
                                    matH[in.npoints*pointA+pointB]+=fA*fB*in.v_value[i*in.ny*in.nz+j*in.nz+k]*dv;
                                }
                                else{
                                    fB=Spline.cal(sqrt(xB+yB+zB));
                                    double tmp=fA*fB*in.v_value[i*in.ny*in.nz+j*in.nz+k]*dv;
                                    matH[in.npoints*pointA+pointB]+=tmp;
                                    matH[in.npoints*pointB+pointA]+=tmp;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    timer::tick("main","integral");
//-----------------------------------------------------------
//EXPLAIN:output the matrix
//-----------------------------------------------------------
    ofstream of("./output/matrix.txt");
    of<<"dimension: "<<in.npoints<<endl;
    for(int i=0;i<in.npoints;i++){
        for(int j=0;j<in.npoints;j++){
            of<<setprecision(10)<<setw(20)<<matH[in.npoints*i+j];
        }
        of<<endl;
    }
    of.close();
    clock.end();
//-----------------------------------------------------------
//EXPLAIN:diagonization
//-----------------------------------------------------------
    diago::lapack_diago(in.npoints,matH);
    delete[] matH;
    delete[] flag;
    return 0;
}