# Large-Scale-Overlapping-Optimization

This is the large-scale overlapping benchmark functions created based on f13 and f14 of the CEC2013 LSGO Benchmark.

There are totally 12 functions here.

The first four functions f1 to f4 are built based on the schwefel function; f5 to f8 are built based on the elliptic function; f9 to f12 are built based on the rastrigin function. Schwefel function and elliptic function are unimodal; Rastrigin is multimodal which is harder to be optimized than the other two functions.

Functions with odd indeces are conforming functions; Functions with even indeces are conflicting functions.

f1, f2, f5, f6, f9, f10 have non-uniform subcomponent sizes. f3, f4, f7, f8, f11, f12 have uniform subcomponent size. Each function has 905 decision variables that can be divided into 20 subcomponents.

f1 and f2 of this benchmark are exactly f13 and f14 of CEC2013 LSGO Benchmark.

6 extra weight files are offered that be used to replace the original weight files. The weight file of each function is "FX-w.txt".
