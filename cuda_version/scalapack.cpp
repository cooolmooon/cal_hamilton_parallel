#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<mpi.h>
#include<iomanip>
#include<stdio.h>

using namespace std;

extern "C" void blacs_get_(int*, int*, int*);
extern "C" void blacs_pinfo_(int*, int*);
extern "C" void blacs_gridinit_(int*, char*, int*, int*);
extern "C" void blacs_gridinfo_(int*, int*, int*, int*, int*);
extern "C" void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern "C" void blacs_gridexit_(int*);
extern "C" int numroc_(int*, int*, int*, int*, int*);
extern "C" void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja,
                 int *desca, double *w, double *z, int *iz, int *jz, int *descz,
                 double *work, int *lwork, int *info);
extern "C" void pdelset_(double* mat, const int* i, const int* j, const int* desc, const double* a);

int main(int argc,char** argv){
    int n;
    double* mat;
    int rank,size;
    int zero=0;
    int one=1;
//-----------------------------------------------------------
//EXPLAIN:set the block parameters
//-----------------------------------------------------------
    int nprow=2;
    int npcol=2;
    int nb=1;
    char jobz='V';
    char uplo='U';
    char layout='R';
//-----------------------------------------------------------
//EXPLAIN:read input matrix
//-----------------------------------------------------------
    fstream in("./output/matrix.txt");
    string name;
    in>>name>>n;
    mat=new double[n*n];
    for(int i=0;i<n*n;i++){
        in>>mat[i];
    }
    in.close();
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    //initialize BLACS
    int ictxt,myrow,mycol;
    int _zero=0;
    blacs_pinfo_(&rank,&size);
    blacs_get_(&_zero,&_zero,&ictxt);
    blacs_gridinit_(&ictxt,&layout,&nprow,&npcol);
    blacs_gridinfo_(&ictxt,&nprow,&npcol,&myrow,&mycol);

    //compute the size of local matrix
    int np=numroc_(&n,&nb,&myrow,&zero,&nprow);//local row
    int nq=numroc_(&n,&nb,&mycol,&zero,&npcol);//local col

    printf("Proc %d/%d for MPI, proc %d/%d for BLACS in position (%d,%d)/(%d,%d) with local matrix %dx%d, global matrix %d, block size %d\n",rank,size,rank,size,myrow,mycol,nprow,npcol,np,nq,n,nb);

    double*A=(double*)calloc(np*nq,sizeof(double));//save the local matrix
    double*Z=(double*)calloc(np*nq,sizeof(double));//save the eigenvectors
    double*W=(double*)calloc(n,sizeof(double));//save the eigenvalues

    //create descriptor
    int descA[9];
    int descZ[9];
    int info;
    int lddA=np>1 ? np : 1;
    descinit_(descA,&n,&n,&nb,&nb,&zero,&zero,&ictxt,&lddA,&info);
    descinit_(descZ,&n,&n,&nb,&nb,&zero,&zero,&ictxt,&lddA,&info);
    //call pdelset routine to distribute block matrix
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            //read matrix
            pdelset_(A,&i,&j,descA,&(mat[j-1+(i-1)*n]));
        }
    }

    double* work=(double*)calloc(2,sizeof(double));
    int lwork=-1;
    pdsyev_(&jobz,&uplo,&n,A,&one,&one,descA,W,Z,&one,&one,descZ,work,&lwork,&info);
    lwork=(int)work[0];
    delete[]work;
    work=(double*)calloc(lwork,sizeof(double));
    if(work==NULL)cout<<"error occurs in pdsyev in proc "<<rank<<endl;
    pdsyev_(&jobz,&uplo,&n,A,&one,&one,descA,W,Z,&one,&one,descZ,work,&lwork,&info);
//-----------------------------------------------------------
//EXPLAIN:creater folder to save the output files
//-----------------------------------------------------------
    string fn="output";
    std::stringstream ss;
    ss<<"test -d "<<fn<<" || mkdir "<<fn;
    system(ss.str().c_str());
    if (rank == 0) {
        //output the result to file
        ofstream valuelog("./output/scal_eigenvalues.log");
        valuelog<<"Eigenvalues are: "<<endl;
        for(int i=0;i<n;i++){
            valuelog<<W[i]<<" ";
        }
    valuelog.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //stringstream log;
    //log<<rank<<"vector.log";
    //ofstream of(log.str().c_str());
    //for(int i=0;i<np;i++){
    //    for(int j=0;j<nq;j++){
    //        of<<Z[np*i+j]<<' ';
    //    }
    //    of<<endl;
    //}
    //of.close();
    //gather the eigenvectors
    int send[4];
    send[0]=np;
    send[1]=nq;
    send[2]=myrow;
    send[3]=mycol;
    int* recv;
    recv=new int[4*size];
    //gather block information
    MPI_Allgather(send,4,MPI_INT,recv,4,MPI_INT,MPI_COMM_WORLD);
    //send the eigenvector to mat
    int start_row=0;
    int start_col=0;
    //suppose the block is divided according to rank
    //get the start position
    int last=-1;
    for(int i=0;i<rank;i++){
        if(recv[4*i+2]<recv[rank*4+2]){
            if(recv[4*i+2]!=last){
                last=recv[4*i+2];
                start_row+=recv[4*i];
            }
        }
    }
    for(int i=0;i<rank;i++){
        start_col+=recv[4*i+1];
    }
    start_col%=n;
    for(int i=0;i<np;i++){
        for(int j=0;j<nq;j++){
            mat[(i+start_row)*n+start_col+j]=Z[j*np+i];
        }
    }
    //send the position information
    int send_pos[2];
    send_pos[0]=start_row;
    send_pos[1]=start_col;
    int* recv_pos;
    recv_pos=new int[2*size];
    MPI_Allgather(send_pos,2,MPI_INT,recv_pos,2,MPI_INT,MPI_COMM_WORLD);
    //send the blocks
    if(rank!=0){
        for(int i=0;i<np;i++){
            MPI_Send(mat+(i+start_row)*n+start_col,nq,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
        }
    }else{
        MPI_Status status;
        for(int i=1;i<size;i++){
            for(int j=0;j<recv[i*4];j++){
                MPI_Recv(mat+(j+recv_pos[2*i])*n+recv_pos[2*i+1],nq,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //output the eigenvectors
    if(rank==0){
        ofstream vectorlog("./output/scal_eigenvectors.log");
        vectorlog<<"eigenvectors are:"<<endl;
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                vectorlog<<setw(15)<<mat[n*j+i];
            }
            vectorlog<<endl;
        }
    }
    //finalize
    delete[]A,W,Z,work,recv,recv_pos;
    blacs_gridexit_(&ictxt);
    MPI_Finalize();
}