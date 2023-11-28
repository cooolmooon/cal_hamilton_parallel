#include "timer.h"
#include <omp.h>
using namespace std;

int timer::max_timer=100;
int timer::num_timer=0;
string*timer::class_name;
string*timer::func_name;
unsigned long long* timer::calls;
double* timer::run_time;
double* timer::start_time;
int timer::start_flag=-1;
bool timer::print_flag=1;
bool timer::delete_flag=false;


//HINT:check whether it needs a ";" here
timer::timer(){}
timer::~timer(){}
//-----------------------------------------------------------
//EXPLAIN:start/end the recording of a function
//-----------------------------------------------------------
void timer::tick(const string&c_name,const string&f_name){
    int timer_location=0;
//-----------------------------------------------------------
//EXPLAIN:find the input name
//-----------------------------------------------------------
    for(;timer_location<num_timer;timer_location++){
        if(c_name==class_name[timer_location]
        &&f_name==func_name[timer_location]){
            break;
        }
    }
//-----------------------------------------------------------
//EXPLAIN:does not find the name, add it to the list
//-----------------------------------------------------------
    if(timer_location==num_timer){
        num_timer++;
        class_name[timer_location]=c_name;
        func_name[timer_location]=f_name;
    }

    if(num_timer>=max_timer){
        cout<<endl<<"Error! Numbers of timers exceed maximum!"<<endl;
        return;
    }

//-----------------------------------------------------------
//EXPLAIN:the start_time is initialized as -1,
//start_time==start_flag means a new call for a function
//NOTE:when a function ends,we should reset start_time as -1
//-----------------------------------------------------------
    if(start_time[timer_location]==start_flag){
        start_time[timer_location]=omp_get_wtime();
        calls[timer_location]++;
    }
    else{
        run_time[timer_location]+=omp_get_wtime()-start_time[timer_location];
        start_time[timer_location]=(double)start_flag;
    }
    return;
}
//-----------------------------------------------------------
//EXPLAIN:start recording time
//-----------------------------------------------------------
void timer::start(void){
    start_time=new double[max_timer];
    run_time=new double[max_timer];
    class_name=new string[max_timer];
    func_name=new string[max_timer];
    calls=new unsigned long long[max_timer];
    delete_flag=true;

    for(int i=0;i<max_timer;i++){
        start_time[i]=(double)start_flag;
        run_time[i]=0;
        class_name[i]="\0";
        func_name[i]="\0";
        calls[i]=0;
    }

    timer::tick("","total");
    return;
}

//-----------------------------------------------------------
//EXPLAIN:stop recording time
//-----------------------------------------------------------
void timer::end(void){
    timer::tick("","total");
    print_all();
    if(delete_flag){
        delete[] start_time;
        delete[] run_time;
        delete[] class_name;
        delete[] func_name;
        delete[] calls;
    }
}


void timer::print_all(void){
    if(!print_flag)return;
    const double smallest = 0.00000000000001;

    cout<<"\n--------- CLASS NAME----------"
        <<" NAME ----- TIME(sec) ---- CALLS ----- AVG ----- PER%"
        <<endl;
    bool* print_flag=new bool[max_timer];
    for(int i=0;i<max_timer;i++) print_flag[i]=false;

    for(int i=0;i<max_timer;i++){
        int k=0;
        double max_tmp=-1.0;
        for(int j=0;j<max_timer;j++){
            //ignore the function that is output already
            if(print_flag[j])continue;

            if(max_tmp<run_time[j]){
                k=j;
                max_tmp=run_time[j];
            }
        }
        print_flag[k]=true;
//-----------------------------------------------------------
//EXPLAIN:ignore the 0 parts
//-----------------------------------------------------------
		if ((run_time[k] >= 0 && run_time[k] < smallest) ||
		        (run_time[k] < 0 && run_time[k] > -smallest))
		{
			continue;
		}


        const long double total_run_time=run_time[k];
        const double average_run_time=total_run_time/calls[k];
        cout.setf(ios::showpoint);
//TODO:the output format may go wrong
        cout<<right
            <<setw(20)<<class_name[k]
            <<setw(15)<<func_name[k]
            <<setw(15)<<setprecision(6)<<total_run_time
            <<setw(10)<<calls[k]
            <<setw(13)<<setprecision(5)<<average_run_time
            <<setw(9)<<fixed<<setprecision(1)<<total_run_time/run_time[0]*100<<"%"
            <<endl;
        cout.unsetf(ios::fixed);    
    }
    cout<<"-----------------------------------------------------------------------------------"<<endl;
    delete[] print_flag;
    return;
}

//-----------------------------------------------------------
//EXPLAIN:set print_flag
//-----------------------------------------------------------
void timer::set(int n){
    if(n)print_flag=true;
    else{
        print_flag=false;
    }
    return;
}

