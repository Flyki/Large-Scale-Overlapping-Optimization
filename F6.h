#ifndef _F6_H
#define _F6_H

#include "Benchmarks.h"

class F6:public Benchmarks{
 public:
  F6();
  double compute(double* x) ;
  ~F6();
};

#endif

