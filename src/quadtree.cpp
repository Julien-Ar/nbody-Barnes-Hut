#include "quadtree.hpp"
#include "particle.hpp"
#include <iostream>
#include <math.h>

QuadTree::QuadTree(){

    quadrant_center = Vect2(0, 0);
    quadrant_top_left = Vect2(0, 0);
    quadrant_bottom_right = Vect2(0, 0);
    size_side = quadrant_bottom_right.axis_x - quadrant_top_left.axis_x;
    //std::cout << "Default construc | size side : " << size_side << std::endl;
    nb_particles = 0;
    total_mass = 0;
    mass_barycenter = Vect2(0, 0);

    vec_particle_ptr = {};

    isLeaf = true;
    update_bary = true;

    TopLeft_subtree = nullptr;
    TopRight_subtree = nullptr;
    BottomLeft_subtree = nullptr;
    BottomRight_subtree = nullptr;

}

QuadTree::QuadTree(const QuadTree & Q){

    quadrant_center = Q.quadrant_center;
    quadrant_top_left = Q.quadrant_top_left;
    quadrant_bottom_right = Q.quadrant_bottom_right;
    size_side = Q.size_side;
    //std::cout << "Const ref construc | size side : " << size_side << std::endl;

    nb_particles = Q.nb_particles;
    total_mass = Q.total_mass;
    mass_barycenter = Q.mass_barycenter;

    vec_particle_ptr = Q.vec_particle_ptr;

    isLeaf = Q.isLeaf;
    update_bary = Q.update_bary;
    TopLeft_subtree = Q.TopLeft_subtree;
    TopRight_subtree = Q.TopRight_subtree;
    BottomLeft_subtree = Q.BottomLeft_subtree;
    BottomRight_subtree = Q.BottomRight_subtree;

}

QuadTree::QuadTree( Vect2  &q_center, 
                    Vect2 &q_TL, 
                    Vect2 &q_BR, 
                    int nb_part, 
                    float tot_mass, 
                    Vect2& mass_bary, 
                    std::vector< std::shared_ptr<Particle> > & vec_part, 
                    bool is_Leaf,
                    bool update_bary1, 
                    std::shared_ptr<QuadTree> TL, 
                    std::shared_ptr<QuadTree> TR , 
                    std::shared_ptr<QuadTree> BL,
                    std::shared_ptr<QuadTree> BR){
                    
    quadrant_center = q_center;
    quadrant_top_left = q_TL;
    quadrant_bottom_right = q_BR;
    size_side = quadrant_bottom_right.axis_x - quadrant_top_left.axis_x;
    //std::cout << "Parameter construc | size side : " << size_side << std::endl;

    nb_particles = nb_part;
    total_mass = tot_mass;
    mass_barycenter = mass_bary;

    vec_particle_ptr = vec_part;

    isLeaf = is_Leaf;
    update_bary = update_bary1;
    TopLeft_subtree = TL;
    TopRight_subtree = TR;
    BottomLeft_subtree = BL;
    BottomRight_subtree = BR;
}




void QuadTree::insert(  std::shared_ptr<Particle>  & particle ){
    
    //std::cout<<"Quad center : ("<<quadrant_center.axis_x <<","<<quadrant_center.axis_x<< ") | Size of size / 2 : " << size_side/2 <<std::endl;

    //std::cout << "Inserting particle pos (" << particle->get_pos_x() <<","<<particle->get_pos_y() <<")"<<std::endl;
    float center_x, center_y, top_left_x, top_left_y, bot_right_x, bot_right_y;
    bool new_isLeaf = true;
    int new_nb_particle = 1;
    bool new_update_bary = true;
    Vect2 new_quadrant_center ;
    Vect2 new_quadrant_top_left;
    Vect2 new_quadrant_bottom_right;
    ++nb_particles;

    if(update_bary){
        mass_barycenter = (mass_barycenter *total_mass  + 
                        particle->position*particle->mass)*(1./(total_mass+particle->mass));  
    }
    total_mass += particle->mass;
    //std::cout << "Pos bary : ("<< mass_barycenter.axis_x<<","<<mass_barycenter.axis_y<<") : nb part :" << nb_particles<<std::endl;

    if(vec_particle_ptr.empty() && (isLeaf==true)){
        // only for root node
        //std::cout << "Root node" << std::endl;
        vec_particle_ptr.push_back(particle);
    }

    else {
        // TopLeft location
        if ((particle->get_pos_x() < quadrant_center.axis_x) && (particle->get_pos_y() > quadrant_center.axis_y)){
            //std::cout << "  Location  : TopLeft" <<std::endl;
            if (TopLeft_subtree == nullptr){    
                //std::cout << "  Creating TopLeft subtree.." << std::endl;
                center_x = quadrant_center.axis_x - (quadrant_center.axis_x-quadrant_top_left.axis_x)/2;
                center_y = quadrant_center.axis_y + (quadrant_top_left.axis_y-quadrant_center.axis_y)/2;
                top_left_x = quadrant_top_left.axis_x;
                top_left_y = quadrant_top_left.axis_y;
                bot_right_x = quadrant_center.axis_x;
                bot_right_y = quadrant_center.axis_y;

                new_quadrant_center.axis_x = center_x, 
                new_quadrant_center.axis_y = center_y, 

                new_quadrant_top_left.axis_x = top_left_x;
                new_quadrant_top_left.axis_y = top_left_y;
                
                new_quadrant_bottom_right.axis_x = bot_right_x;
                new_quadrant_bottom_right.axis_y = bot_right_y;
                

                std::vector < std::shared_ptr<Particle> > vec_particle_ptr_TL {particle}; 
   
                TopLeft_subtree = std::make_shared<QuadTree> (QuadTree(new_quadrant_center, new_quadrant_top_left, new_quadrant_bottom_right,
                                                                        new_nb_particle, particle->mass, particle->position,
                                                                        vec_particle_ptr_TL, new_isLeaf, new_update_bary ,nullptr, nullptr, nullptr, nullptr)) ;

            }
            else{
                //std::cout << "  Existing particle in TopLeft, inserting in TopLeft..." << std::endl;
                TopLeft_subtree->insert(particle);
            }       
        }
        // TopRight
        else if((particle->get_pos_x() > quadrant_center.axis_x) && (particle->get_pos_y() > quadrant_center.axis_y)) {
            //std::cout << "  Location  : TopRight" <<std::endl;
            if (TopRight_subtree == nullptr){  

                //std::cout << "  Creating TopRight subtree.." << std::endl;
                center_x = quadrant_center.axis_x + (quadrant_bottom_right.axis_x-quadrant_center.axis_x)/2;
                center_y = quadrant_center.axis_y + (quadrant_top_left.axis_y-quadrant_center.axis_y)/2;
                top_left_x = quadrant_center.axis_x;
                top_left_y = quadrant_top_left.axis_y;
                bot_right_x = quadrant_bottom_right.axis_x;
                bot_right_y = quadrant_center.axis_y;

                new_quadrant_center.axis_x = center_x, 
                new_quadrant_center.axis_y = center_y, 

                new_quadrant_top_left.axis_x = top_left_x;
                new_quadrant_top_left.axis_y = top_left_y;
                
                new_quadrant_bottom_right.axis_x = bot_right_x;
                new_quadrant_bottom_right.axis_y = bot_right_y;
                
                std::vector < std::shared_ptr<Particle> > vec_particle_ptr_TR {particle}; 

                TopRight_subtree = std::make_shared<QuadTree> (QuadTree(new_quadrant_center, new_quadrant_top_left, new_quadrant_bottom_right,
                                                                         new_nb_particle, particle->mass, particle->position,
                                                                        vec_particle_ptr_TR, new_isLeaf, new_update_bary,  nullptr, nullptr, nullptr, nullptr)) ;

            }
            else{
                //std::cout << "  Existing particle in TopRight, inserting in TopRight..." << std::endl;
                TopRight_subtree->insert(particle);
            } 

        }
        // BottomLeft
        else if((particle->get_pos_x() < quadrant_center.axis_x) && (particle->get_pos_y() < quadrant_center.axis_y)) {
            //std::cout << "  Location  : BottomLeft" <<std::endl;
            if (BottomLeft_subtree == nullptr){    
                //std::cout << "  Creating BottomLeft subtree.." << std::endl;

                center_x = quadrant_center.axis_x - (quadrant_center.axis_x-quadrant_top_left.axis_x)/2;
                center_y = quadrant_center.axis_y - (quadrant_center.axis_y-quadrant_bottom_right.axis_y)/2;
                top_left_x = quadrant_top_left.axis_x;
                top_left_y = quadrant_center.axis_y;
                bot_right_x = quadrant_center.axis_x;
                bot_right_y = quadrant_bottom_right.axis_y;
                
                new_quadrant_center.axis_x = center_x, 
                new_quadrant_center.axis_y = center_y, 

                new_quadrant_top_left.axis_x = top_left_x;
                new_quadrant_top_left.axis_y = top_left_y;
                
                new_quadrant_bottom_right.axis_x = bot_right_x;
                new_quadrant_bottom_right.axis_y = bot_right_y;
                
                std::vector < std::shared_ptr<Particle> > vec_particle_ptr_BL {particle}; 

                BottomLeft_subtree = std::make_shared<QuadTree> (QuadTree(new_quadrant_center, new_quadrant_top_left, new_quadrant_bottom_right,
                                                                        new_nb_particle, particle->mass, particle->position,
                                                                        vec_particle_ptr_BL, new_isLeaf, new_update_bary, nullptr, nullptr, nullptr, nullptr)) ;

            }
            else{
                //std::cout << "  Existing particle in BottomLeft, inserting in BottomLeft..." << std::endl;
                BottomLeft_subtree->insert(particle);
            } 
         }
        // BottomRight
        else if((particle->get_pos_x() > quadrant_center.axis_x) && (particle->get_pos_y() < quadrant_center.axis_y)) {
            //std::cout << "  Location  : BottomRight" <<std::endl;
            if (BottomRight_subtree == nullptr){    
                //std::cout << "  Creating BottomRight subtree.." << std::endl;

                center_x = quadrant_center.axis_x + (quadrant_bottom_right.axis_x - quadrant_center.axis_x)/2;
                center_y = quadrant_center.axis_y - (quadrant_center.axis_y - quadrant_bottom_right.axis_y)/2;
                top_left_x = quadrant_center.axis_x;
                top_left_y = quadrant_center.axis_y;
                bot_right_x = quadrant_bottom_right.axis_x;
                bot_right_y = quadrant_bottom_right.axis_y;

                new_quadrant_center.axis_x = center_x, 
                new_quadrant_center.axis_y = center_y, 

                new_quadrant_top_left.axis_x = top_left_x;
                new_quadrant_top_left.axis_y = top_left_y;
                
                new_quadrant_bottom_right.axis_x = bot_right_x;
                new_quadrant_bottom_right.axis_y = bot_right_y;
                
                std::vector < std::shared_ptr<Particle> > vec_particle_ptr_BR {particle}; 

                BottomRight_subtree = std::make_shared<QuadTree> (QuadTree(new_quadrant_center, new_quadrant_top_left, new_quadrant_bottom_right,
                                                                         new_nb_particle, particle->mass, particle->position,
                                                                        vec_particle_ptr_BR, new_isLeaf, new_update_bary, nullptr, nullptr, nullptr, nullptr)) ;

            }
            else{
                //std::cout << "  Existing particle in BottomRight, inserting in BottomRight..." << std::endl;
                BottomRight_subtree->insert(particle);
            } 
        }
        // we have to reaffect the previous particle occupying this quadrant a new quadrant 
        if (vec_particle_ptr.size() != 0 ) {
            auto & previous_particle = vec_particle_ptr[0];
            vec_particle_ptr.pop_back(); // the if-condition prevents accessing vec_particle[0],which doesnt exist at the next call of insert
            --nb_particles; // incremented a next call 
            total_mass -= previous_particle->mass; // incremented with same masse at next call
            update_bary = false; // the current barycenter is correct, no need to update it
            isLeaf = false; // not a leaf anymore 
            insert(previous_particle);
            update_bary = true;
            
        }
         
    }   
}

void QuadTree::computeAccelerationOn(std::shared_ptr<Particle> & particle, float theta){


    
    float dist_from_bary = pow(pow(mass_barycenter.axis_x - particle->position.axis_x, 2) + 
                                pow(mass_barycenter.axis_y - particle->position.axis_y, 2), 0.5);
    // std::cout << "Dist from bary : " << dist_from_bary<< std::endl;
    // std::cout << "size side : " << size_side << std::endl;

    if(isLeaf){
        auto & p1 = vec_particle_ptr[0];
        if(p1 != particle){
            //std::cout <<"Different particle" << std::endl;
            particle->computeAccelerationDueTo(p1);
        }
    }
    else{
        float MAC = size_side/dist_from_bary; // Multipole-Acceptance-Criterion
        if(MAC < theta){
            //std::cout << "MAC < theta : " << MAC << "| : "<<theta <<std::endl;

            std::shared_ptr<Particle> mass_bary_particle = std::make_shared<Particle>(Particle(mass_barycenter, 
                                                                                                Vect2(0,0), 
                                                                                                Vect2(0,0), 
                                                                                                total_mass));
            particle->computeAccelerationDueTo(mass_bary_particle);
        }
        else{
            //std::cout << "MAC > theta : " << MAC << "| : "<<theta <<std::endl;
            if(TopLeft_subtree != nullptr){
                TopLeft_subtree->computeAccelerationOn(particle, theta);
            }
            if(BottomRight_subtree != nullptr){
                BottomRight_subtree->computeAccelerationOn(particle, theta);
            }
            if(BottomLeft_subtree != nullptr){
                BottomLeft_subtree->computeAccelerationOn(particle, theta);
            }
            if(BottomRight_subtree != nullptr){
                BottomRight_subtree->computeAccelerationOn(particle, theta);
            }        
        }
    }
}


















