## Nbody simulation using Barnes-hut algorithm and a pointer based quadtree

Dependencies : [*fmt*](https://github.com/fmtlib/fmt),  [*openmp*](https://www.openmp.org/)

 - fmt : used to generate name for files containing generated initial conditions

- openmp : used to parallelize the force calculations once the quadtree is built, each iteration


The hardest part of any Nbody simulation is setting the inital conditions.

Here the user can chose among three different way of settings :
- 1) a slighlty modified [Plummer model](https://en.wikipedia.org/wiki/Plummer_model) in order to achieve some 
sort of stability and not having a majority of particles flying off to infinity

- 2) Particles having only orbital velocities and orbiting a heavy mass at the center

- 3) Particles having no initial velocties, commonly called "cold start"

## To compile : 

`git clone https://github.com/Julien-Ar/nbody-Barnes-Hut.git`

`cd nbody-Barnes-Hut && mkdir build`

`cd build`

`cmake ..`

`make`

To execute :

`./main`

## To add next :
- Clean up code and add user inputs about the simulation's parameters
- [Morton code](https://en.wikipedia.org/wiki/Z-order_curve) encoding and sorting of particles' positions
in order to paralellize construction of the quadtree

- OpenGL realtime rendering or after calculations are completed
