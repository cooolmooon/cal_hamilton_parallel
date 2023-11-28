#ifndef TIMER_H
#define TIMER_H

#include<ctime>
#include<string>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<iostream>
using namespace std;

//-----------------------------------------------------------
//EXPLAIN:definition of class "timer.h"
//used for recording the running time of different functions
//-----------------------------------------------------------

class timer{
    public:
    timer();
    ~timer();

//-----------------------------------------------------------
//EXPLAIN:used twice for a function
//-----------------------------------------------------------
    static void tick(const string& _class_name,const string& _func_name);

//-----------------------------------------------------------
//EXPLAIN:start and end a timer
//record the total running time
//-----------------------------------------------------------
    static void start(void);
    static void end(void);

//-----------------------------------------------------------
//EXPLAIN:output all the timer info
//-----------------------------------------------------------
    static void print_all(void);
    static void set(int n);
    private:

//-----------------------------------------------------------
//EXPLAIN:
//max_clock:maximum numbers of timers
//num_clock:numbers of timers now
//calls:the number of calls for a given function
//-----------------------------------------------------------

    static int max_timer;
    static int num_timer;
    static int start_flag;
    static bool print_flag;
    static bool delete_flag;
    static string*class_name;
    static string*func_name;
    static unsigned long long*calls;
    static double*start_time;
    static double*run_time;
};

#endif