#ifndef DIAGO_H
#define DIAGO_H

using namespace std;

class diago{
    public:
    diago();
    ~diago();

    static void lapack_diago(int n,double* mat);
};

#endif