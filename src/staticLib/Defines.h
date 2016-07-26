#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <set>
#include <string.h>
#include <math.h>
#include <iomanip>
#if defined(_WIN32) || defined(_WIN64)
#include <io.h>
#include <direct.h>
#include <Windows.h>
#define strcasecmp _stricmp         // use _stricmp() in Windows
#pragma warning(disable : 4996)     // ignore MSVC's concerns about safety of IO library
#endif
using namespace std;

//***********************************************************//
//**********          macro switch define          **********//
//***********************************************************//
#define KAKOU
//#define ADD_ATTRIBUTE
#define USE_SIFT_HIST

#define GCC_EXPANSION

//#define CDVS_SW_BITANSWER

//***********************************************************//
//**********          macro value define           **********//
//***********************************************************//
#define CDVS_PI     3.14159265358979323846f
#define CDVS_PI2    (CDVS_PI * 2)
//#define CDVS_INDEX_LIMIT    10000
#ifdef CDVS_SW_BITANSWER
    #define FEATURE_VERIFY_NUM   100
    #define RETRIEVE_VERIFY_NUM  100
#endif
#define MAX_MATCHES         500
#define MAX_TRUTHES         500
#define MAX_FILENAME_LENGTH 256
#define MAX_STRING_LENGTH   (MAX_FILENAME_LENGTH * MAX_TRUTHES)

#define MAX_DSC_LENGTH      (32 * 1024)			// Maximum size of descriptor buffer.
#define HM_SIFT_QUANT_SIZE  128				    // The size of quantized data (128 = 16*8, H mode)
#define MAX_NUM_FEATURES    1600

//macro typedef
typedef char FILENAME[MAX_FILENAME_LENGTH];

//***********************************************************//
//**********      macro operation define           **********//
//***********************************************************//
#define CDVS_MIN(a,b)  ((a) > (b) ? (b) : (a))
#define CDVS_MAX(a,b)  ((a) < (b) ? (b) : (a))
#define CDVS_ABS(a)    ((a) > 0 ? (a) : -(a))
#define log2(x)        (log(x) / VL_LOG_OF_2)

#if defined(_OPENMP)
    #define OMP_SET_LOCK 		                    \
        if(omp_in_parallel())   omp_set_lock(&lock);
    #define OMP_UNSET_LOCK 		                    \
        if(omp_in_parallel())   omp_unset_lock(&lock);
    #define OMP_INIT_LOCK                           \
        omp_init_lock(&lock);
    #define OMP_DESTROY_LOCK                        \
        omp_destroy_lock(&lock);
    #define OMP_SET_THREAD_NUMS(nthreads)           \
        omp_set_num_threads(nthreads)
#else
    #define OMP_SET_LOCK
    #define OMP_INIT_LOCK
    #define OMP_UNSET_LOCK
    #define OMP_DESTROY_LOCK
    #define OMP_SET_THREAD_NUMS
#endif
