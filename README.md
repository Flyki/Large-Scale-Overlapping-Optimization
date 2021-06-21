# Large-Scale-Overlapping-Optimization

1) This is the large-scale overlapping benchmark functions created based on f13 and f14 of the CEC2013 LSGO Benchmark with the CBCCO algorithm.
Here are totally 12 functions.

2) The first four functions f1 to f4 are built based on the schwefel function; f5 to f8 are built based on the elliptic function; f9 to f12 are built based on the rastrigin function. Schwefel function and elliptic function are unimodal; Rastrigin is multimodal which is harder to be optimized than the other two functions.

3) Functions with odd indeces are conforming functions; Functions with even indeces are conflicting functions.

4) f1, f2, f5, f6, f9, f10 have non-uniform subcomponent sizes. f3, f4, f7, f8, f11, f12 have uniform subcomponent size. Each function has 905 decision variables that can be divided into 20 subcomponents.

5) f1 and f2 of this benchmark are exactly f13 and f14 of CEC2013 LSGO Benchmark.

6) 6 extra weight files are offered that can be used to replace the original weight files. The weight file of each function is "FX-w.txt".

7) CBOG_CBD represents the CBCCO algorithm and CBOCC is the main file.

8) Two CMAESO files are just used to packing the CMA-ES algorithm. The original code of CMA-ES can be found http://cma.gforge.inria.fr/cmaes_sourcecode_page.html#C. Please use the C++ version.

9) An example of the partition file is shown in 1po.txt and 1oo.txt which represent the grouping file and overlapping variable file, respectively.
