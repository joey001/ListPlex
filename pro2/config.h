#pragma once
#define CAP_LIMIT1 (1<<12)
#define CAP_LIMIT2 (1<<12)
#define PSIZE_LIMIT  512 //Maximum possible size of K-plex
#define K_LIMIT 4

// cilkarts cilk++
#if defined(CILK)
#include <cilk.h>
#include <cassert>
#define parallel_main cilk_main
#define parallel_for cilk_for
#define parallel_for_1 _Pragma("cilk_grainsize = 1") cilk_for
#define parallel_for_256 _Pragma("cilk_grainsize = 256") cilk_for
static int getWorkers() { return -1; }
static void setWorkers(int n) { }

// intel cilkplus
#elif defined(CILKP)
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <sstream>
#include <iostream>
#include <cstdlib>
#define parallel_for cilk_for
#define parallel_main main
// #pragma cilk grainsize = value


static int getWorkers() {
  return __cilkrts_get_nworkers();
}
static void setWorkers(int n) {
  __cilkrts_end_cilk();
  //__cilkrts_init();
  std::stringstream ss; ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}
static int getNumber(){
  return __cilkrts_get_worker_number();
}

// openmp
#elif defined(OPENMP)
#include <omp.h>
#define parallel_main main
#define cilk_spawn
#define cilk_sync 
#define parallel_for _Pragma("omp parallel for") for
#define parallel_for_1 _Pragma("omp parallel for schedule (static,1)") for
#define parallel_for_256 _Pragma("omp parallel for schedule (static,256)") for

static int getWorkers() { return omp_get_max_threads(); }
static void setWorkers(int n) { omp_set_num_threads(n); }
static int getNumber(){ return omp_get_thread_num(); }

// serial c++
#else
#define parallel_main main
#define parallel_for for
#define parallel_for_1 for
#define parallel_for_256 for
#define cilk_sync  
#define cilk_spawn  

static int getWorkers() { return 1; }
static void setWorkers(int n) { }
static int getNumber(){ return 0; }

#endif


#include <limits.h>
#if defined(LONG)
typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#define UINT_T_MAX ULONG_MAX
#else
typedef int intT;
typedef unsigned int uintT;
#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX
#endif


#ifdef PROMPT
#include <cstdarg>
#define prompt(format, ...) \
	printf("[INFO]: %s:%d: " format , __FILE__, __LINE__, ##__VA_ARGS__)
#define promptln(format, ...) \
	printf("[INFO]: %s:%d: \n" format , __FILE__, __LINE__, ##__VA_ARGS__)
#define shortprompt(format, ...) \
    printf(format, ##__VA_ARGS__)
#define showContext()  printContextInfo()
#else
#define prompt(format, ...) 
#define promptln(format, ...)
#define shortprompt(format, ...)
#define showContext()
#endif