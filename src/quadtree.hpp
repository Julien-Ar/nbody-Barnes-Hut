#ifndef QUADTREE_H
#define QUADTREE_H

#include <memory>
#include <vector>
#include "particle.hpp"
#include "vect2.hpp"


class QuadTree {
    public:

        Vect2 quadrant_center;
        Vect2 quadrant_top_left;
        Vect2 quadrant_bottom_right;

        float size_side;
        // total number of particles inside this quadrant
        int nb_particles;
        float total_mass;

        Vect2 mass_barycenter;
        

        // store adresses of particles that are inside that node 
        std::vector< std::shared_ptr<Particle> > vec_particle_ptr;

        bool isLeaf;
        bool update_bary;
        std::shared_ptr<QuadTree> TopLeft_subtree;
        std::shared_ptr<QuadTree> TopRight_subtree;
        std::shared_ptr<QuadTree> BottomLeft_subtree;
        std::shared_ptr<QuadTree> BottomRight_subtree;



        QuadTree();
        QuadTree(const QuadTree & );
        QuadTree( Vect2 &, 
                  Vect2 &, 
                  Vect2 &, 
                  int, 
                  float, 
                  Vect2 &, 
                  std::vector<std::shared_ptr<Particle>> &, 
                  bool, 
                  bool, 
                  std::shared_ptr<QuadTree>, 
                  std::shared_ptr<QuadTree> , 
                  std::shared_ptr<QuadTree>  , 
                  std::shared_ptr<QuadTree> );

        void insert( std::shared_ptr<Particle> &  );
        void computeAccelerationOn(std::shared_ptr<Particle> &, float theta);
        
};



#endif