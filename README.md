## Nbody simulation using Barnes-hut algorithm and a pointer based quadtree

Dependencies : [*fmt*](https://github.com/fmtlib/fmt), [*openmp*](https://www.openmp.org/)

fmt : used to generate name for files containing generated initial conditions

openmp : used to parallelize the force calculations once the quadtree is built


To compile : 

`cd build`

`cmake ..`

`make`
