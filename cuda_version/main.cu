#include<iostream>
#include"./tools/input.h"
#include"./tools/spline.h"
#include<sstream>
#include<cmath>
#include<stdio.h>
#include<iomanip>
#include"cuda_runtime.h"

using namespace std;
//-----------------------------------------------------------
//EXPLAIN:constant used by interplation
//-----------------------------------------------------------
__device__ double dr;
__device__ int cutoff; 
__device__ int width;
__device__ double fvalues[2000];
__device__ double m[2000];
__device__ int dev_point[2500]; 
__device__ float sum=0.0;
//-----------------------------------------------------------
//EXPLAIN:constant used by integral
//-----------------------------------------------------------
__device__ double tx;
__device__ double ty;
__device__ double tz;
__device__ int npoint;
__device__ int nx;
__device__ int ny;
__device__ int nz;
//-----------------------------------------------------------
//EXPLAIN:interplation
//-----------------------------------------------------------
__device__ double interplate(double r){
    int n=floor(r/dr);
    double p1=r/dr-n;
    double p2=n+1-r/dr;
    double ans=(1+2*p1)*p2*p2*fvalues[n]+(1+2*p2)*p1*p1*fvalues[n+1]
                +(r-n*dr)*p2*p2*m[n]+(r-(n+1)*dr)*p1*p1*m[n+1];
    return ans;
}
//-----------------------------------------------------------
//EXPLAIN:test
//-----------------------------------------------------------
__global__ void print_data(int* dev_area){
    printf("(%d,%d),(%d,%d),(%d,%d),\n",dev_point[0],dev_point[1],dev_point[2],dev_point[3],dev_point[4],dev_point[5]);
    printf("sum=%f\n",sum);
}
//-----------------------------------------------------------
//EXPLAIN:integral
//-----------------------------------------------------------
__global__ void integral(int* dev_area,int*dev_x,int* dev_y,int* dev_z,double* dev_V,float* dev_matH){
    //allocate dynamic SMEM
    // extern __shared__ int arr[];
    int k=blockIdx.x;
    int thdid=threadIdx.x;
    int pA=dev_point[2*k];
    int pB=dev_point[2*k+1];
    int sx=dev_area[npoint*pA+pB*6];
    int sy=dev_area[npoint*pA+pB*6+2];
    int sz=dev_area[npoint*pA+pB*6+4];
    int task_num=width*width/4;
    int start=thdid*task_num;
    // arr[0]=sx;
    // arr[1]=sy;
    // arr[2]=sz;
    // arr[3]=width;
    // arr[4]=dev_x[pA];
    // arr[5]=dev_y[pA];
    // arr[6]=dev_z[pA];
    // arr[7]=dev_x[pB];
    // arr[8]=dev_y[pB];
    // arr[9]=dev_z[pB];
    // arr[10]=ny;
    // arr[11]=nz;
    //printf("%d,%d,%d,%d,%d,%d\n",k,thdid,pA,pB,task_num,start);
    for(int i=0;i<task_num;i++){
        // int x=((start+i)/(arr[3]*arr[3])+arr[0]);
        // int y=(((start+i)%(arr[3]*arr[3]))/arr[3]+arr[1]);
        // int z=((start+i)%arr[3]+arr[2]);
        // //if(thdid==width)printf("block:%d, id:%d,%d, %d, %d\n",k,thdid,x,y,z);
        // double rA=sqrt((x*tx-arr[4])*(x*tx-arr[4])+
        // (y*ty-arr[5])*(y*ty-arr[5])+
        // (z*tz-arr[6])*(z*tz-arr[6]));
        // double rB=sqrt((x*tx-arr[7])*(x*tx-arr[7])+
        // (y*ty-arr[8])*(y*ty-arr[8])+
        // (z*tz-arr[9])*(z*tz-arr[9]));
        // double fA=interplate(rA);
        // double fB=interplate(rB);
        // float tmp=fA*fB*dev_V[x*arr[10]*arr[11]+y*arr[11]+z]*tx*ty*tz;
        int x=((start+i)/(width*width)+sx);
        int y=(((start+i)%(width*width))/width+sy);
        int z=((start+i)%width+sz);
        //if(thdid==width)printf("block:%d, id:%d,%d, %d, %d\n",k,thdid,x,y,z);
        double rA=sqrt((x*tx-dev_x[pA])*(x*tx-dev_x[pA])+
        (y*ty-dev_y[pA])*(y*ty-dev_y[pA])+
        (z*tz-dev_z[pA])*(z*tz-dev_z[pA]));
        double rB=sqrt((x*tx-dev_x[pB])*(x*tx-dev_x[pB])+
        (y*ty-dev_y[pB])*(y*ty-dev_y[pB])+
        (z*tz-dev_z[pB])*(z*tz-dev_z[pB]));
        double fA=interplate(rA);
        double fB=interplate(rB);
        float tmp=fA*fB*dev_V[x*ny*nz+y*nz+z]*tx*ty*tz;
        atomicAdd(dev_matH+pA*npoint+pB,tmp);
        __syncthreads();
    }
}
//-----------------------------------------------------------
//EXPLAIN:locate the start and end pos
//-----------------------------------------------------------
void compare(const int&x1,const int&x2,const int&cutoff,const int&lx,int&sx,int&ex,
            const int& width,const double& tx,const int&nx){
    sx=x1-cutoff;
    if(sx<0){
        sx=0;
        ex=width;
    }else{
        sx=(int)sx/tx;
        ex=sx+width;
        if(ex>nx-1){
            ex=nx-1;
            sx=ex-width;
        }
    }
}
int main()
{
//-----------------------------------------------------------
//EXPLAIN:read input
//-----------------------------------------------------------
    input in("./input/INPUT.txt");
    spline Spline(in.distribution_path);
//-----------------------------------------------------------
//EXPLAIN:copy static data to device
//-----------------------------------------------------------
    size_t bytes1=Spline.mesh*sizeof(double);
    cudaMemcpyToSymbol(dr,&Spline.dr,sizeof(double));
    cudaMemcpyToSymbol(cutoff,&Spline.cutoff,sizeof(int));
    cudaMemcpyToSymbol(fvalues,Spline.fvalues,bytes1);
    cudaMemcpyToSymbol(m,Spline.m,bytes1);
    double host_tx=(double)in.lx/in.nx;
    double host_ty=(double)in.ly/in.ny;
    double host_tz=(double)in.lz/in.nz;
    cudaMemcpyToSymbol(tx,&host_tx,sizeof(double));
    cudaMemcpyToSymbol(ty,&host_ty,sizeof(double));
    cudaMemcpyToSymbol(tz,&host_tz,sizeof(double));
    int host_width=2*Spline.cutoff/host_tx+2;
    if(host_width%2==1)host_width+=1;
    cudaMemcpyToSymbol(width,&host_width,sizeof(int));
//-----------------------------------------------------------
//EXPLAIN:preprocess
//-----------------------------------------------------------
    int* cal_area;
    int cal_point[2500];//2k stand for i, 2k+1 stand for j, k is blockIdx
    cal_area=new int[6*in.npoints*in.npoints];//the area that might intersect:160x160x160
    int cal_cnt=0;
    for(int i=0;i<6*in.npoints*in.npoints;i++){
        cal_area[i]=0;
    }
    for(int i=0;i<in.npoints;i++){
        int x1=in.pos_x[i];
        int y1=in.pos_y[i];
        int z1=in.pos_z[i];
        for(int j=i;j<in.npoints;j++){
            int x2=in.pos_x[j];
            int y2=in.pos_y[j];
            int z2=in.pos_z[j];
            if(abs(x1-x2)>Spline.cutoff||
            abs(y1-y2)>Spline.cutoff||
            abs(z1-z2)>Spline.cutoff){
                continue;
            }
            int sx,sy,sz,ex,ey,ez;
//-----------------------------------------------------------
//EXPLAIN:locate the start and end pos
//-----------------------------------------------------------
            compare(x1,x2,Spline.cutoff,in.lx,sx,ex,host_width,host_tx,in.nx);
            compare(y1,y2,Spline.cutoff,in.ly,sy,ey,host_width,host_tx,in.nx);
            compare(z1,z2,Spline.cutoff,in.lz,sz,ez,host_width,host_tx,in.nx);
            cal_area[in.npoints*i+6*j]=sx;
            cal_area[in.npoints*i+6*j+1]=ex;
            cal_area[in.npoints*i+6*j+2]=sy;
            cal_area[in.npoints*i+6*j+3]=ey;
            cal_area[in.npoints*i+6*j+4]=sz;
            cal_area[in.npoints*i+6*j+5]=ez;
            cal_point[cal_cnt*2]=i;
            cal_point[cal_cnt*2+1]=j;
            cal_cnt+=1;
        }
    }
    // for(int i=0;i<in.npoints;i++){
    //     for(int j=i;j<in.npoints;j++){
    //         for(int k=0;k<6;k++){
    //             std::cout<<cal_area[in.npoints*i+6*j+k]<<" ";
    //         }
    //         std::cout<<endl;
    //     }
    // }
    int* dev_area;
    size_t bytes2=6*in.npoints*in.npoints*sizeof(int);
    cudaMalloc((int**)&dev_area,bytes2);
    cudaMemcpy(dev_area,cal_area,bytes2,cudaMemcpyHostToDevice);
    size_t bytes3=2000*sizeof(int);
    cudaMemcpyToSymbol(dev_point,cal_point,bytes3);
    cudaMemcpyToSymbol(npoint,&in.npoints,sizeof(int));
    cudaMemcpyToSymbol(nx,&in.nx,sizeof(int));
    cudaMemcpyToSymbol(ny,&in.ny,sizeof(int));
    cudaMemcpyToSymbol(nz,&in.nz,sizeof(int));
    int *dev_x;
    int *dev_y;
    int *dev_z;
    size_t bytes6=in.npoints*sizeof(int);
    cudaMalloc((int**)&dev_x,bytes6);
    cudaMalloc((int**)&dev_y,bytes6);
    cudaMalloc((int**)&dev_z,bytes6);
    cudaMemcpy(dev_x,in.pos_x,bytes6,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_y,in.pos_y,bytes6,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_z,in.pos_z,bytes6,cudaMemcpyHostToDevice);
//-----------------------------------------------------------
//EXPLAIN:matrix H
//-----------------------------------------------------------
    float* host_matH;
    float* dev_matH;
    host_matH=new float[in.npoints*in.npoints];
    size_t bytes4=in.npoints*in.npoints*sizeof(float);
    cudaMalloc((float**)&dev_matH,bytes4);
    double* dev_V;
    size_t bytes5=in.nx*in.ny*in.nz*sizeof(double);
    cudaMalloc((double**)&dev_V,bytes5);
    cudaMemcpy(dev_V,in.v_value,bytes5,cudaMemcpyHostToDevice);
	//number of blocks
	int num_blk=cal_cnt;
	//number of threads
	int num_thd=host_width*4;
    // int num_thd=512;
    // print_data<<<1,1>>>(dev_area);
//-----------------------------------------------------------
//EXPLAIN:integral
//-----------------------------------------------------------
    cudaEvent_t start,stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);
    integral<<<num_blk,num_thd>>>(dev_area,dev_x,dev_y,dev_z,dev_V,dev_matH);
    cudaDeviceSynchronize();
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime,start,stop);
    std::cout<<"time used for integral is: "<<elapsedTime<<"ms"<<endl;
    // print_data<<<1,1>>>(dev_area);
//-----------------------------------------------------------
//EXPLAIN:copy result to host
//-----------------------------------------------------------
    cudaMemcpy(host_matH,dev_matH,bytes4,cudaMemcpyDeviceToHost);
    // float host_sum=2;
    // cudaMemcpyFromSymbol(&host_sum,sum,sizeof(float));
    // std::cout<<"sum= "<<host_sum<<endl;
    for(int i=0;i<in.npoints;i++){
        for(int j=0;j<i;j++){
            host_matH[i*in.npoints+j]=host_matH[j*in.npoints+i];
        }
    }
//-----------------------------------------------------------
//EXPLAIN:creater folder to save the output files
//-----------------------------------------------------------
    string fn="output";
    std::stringstream ss;
    ss<<"test -d "<<fn<<" || mkdir "<<fn;
    system(ss.str().c_str());
//-----------------------------------------------------------
//EXPLAIN:output result
//-----------------------------------------------------------
    ofstream ofs("./output/matrix.txt");
    ofs<<"dimension: "<<in.npoints<<endl;
    for(int i=0;i<in.npoints;i++){
        for(int j=0;j<in.npoints;j++){
            ofs<<setw(15)<<host_matH[i*in.npoints+j]<<" ";
        }
        ofs<<endl;
    }
    cudaFree(dev_area);
    cudaFree(dev_matH);
    cudaFree(dev_V);
	cudaDeviceReset();
    delete[] cal_area;
    delete[] host_matH;

    return 0;
}