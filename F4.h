#ifndef _F4_H
#define _F4_H

#include "Benchmarks.h"

class F4:public Benchmarks{
 public:
  F4();
  double compute(double* x) ;
  ~F4();
};

#endif

