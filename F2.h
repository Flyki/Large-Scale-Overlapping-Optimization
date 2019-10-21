#ifndef _F2_H
#define _F2_H

#include "Benchmarks.h"

class F2:public Benchmarks{
 public:
  F2();
  double compute(double* x) ;
  ~F2();
};

#endif

