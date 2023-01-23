#include "utils.hpp"
#include "quadtree.hpp"
#include <omp.h>
#include <chrono>
#include <ctime>


// gravitational constant with mass = 1 Mo and distance = 1 pc and velocities in km/s
float G = 0.0043;
//float G = 1; // in henon units 
float epsilon = 100; // softening parameer in force calculations 
                    // (usualy < 1, but set at a higher value to keep
                    // particles from having a too great accel from close encounters)

auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());


int main(){
    std::cout<< ctime(&timenow) << std::endl;

    // ------------- Initial conditions model parameters ---------------

    // ------------- Mandatory parameters -------------------
    constexpr float M = 1e8;        // total mass of particles with masses (black hole not counted)
    constexpr float Nbody = 10000;  // Number of particle with masses that interact with each other
    constexpr float Ntracer = 5000; // Nb of particle without masses that interacts only with particles that have mass
    constexpr float theta = 1;      // threshold in Barnes-hut algo below which controls wether or not to approximate
                                    // a cluster of body as a single body in the force calculation

    // ------------- Plummer model parameters ---------------
    constexpr float a = 50; // radius kpc in plummer model
    constexpr float a_multiplier = 1; // for particle without masses

    // ------------- Cold start and kepler model parameters ---------------
    constexpr float R = 1000; //radius of the galaxy with a black hole at the center


    float max_r = a*a_multiplier*pow(pow(0.99, -2./3)-1, -0.5); // max radius allowed por particle's position
    //float max_r = 10;
    std::vector< std::shared_ptr<Particle> > vec_ptr_particles = create_initial_conditions_plummer( a, 
                                                                                           M,
                                                                                            Nbody,
                                                                                            Ntracer,
                                                                                              a_multiplier,
                                                                                            6,
                                                                                            true);

    //std::vector< std::shared_ptr<Particle> > vec_ptr_particles = read_IC_from_file("../src/behalf_initCond_50a_100000000000000.0Mt_10000N.txt");    
    // std::vector< std::shared_ptr<Particle> > vec_ptr_particles = create_initial_conditions_cold_start(M, 
    //                                                                                                   Nbody,
    //                                                                                                   Ntracer, 
    //                                                                                                   R);
    

    constexpr float Tmax = 30;
    float t= 0;
    constexpr float dt = 0.01;
    //constexpr int save_frequence = int(1/dt);
    constexpr int save_frequence = 2;
    constexpr float iter_max = Tmax/dt;

    int iteration = 0;
    float iteration_percent = 0.0;
    bool launch_simu = true;
    bool bruteforce = false;
        // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();
 
    if(launch_simu){
        if(bruteforce){
            // Regular way of computing pairwise interactions in O(N²)
            while ( t < Tmax ){
                std::cout << "t : " << t <<std::endl;

                for ( auto &p1 : vec_ptr_particles){
                    p1->computePositionAtHalfTimeStep(dt); // drift
                }
                // i and j are used to skip the iteration in which p1 = p2
                int i = 0;
                int j;
                for ( auto & p1 : vec_ptr_particles){
                    p1->resetAcceleration();
                    j = 0;
                    for(auto & p2 : vec_ptr_particles){
                        if(i != j){
                            if(p2->mass > 0){
                                p1->computeAccelerationDueTo(p2);
                            }
                        }
                        ++j;
                    }
                    ++i;
                    p1->computeVelocityAndPosition(dt); // Kick and Drift
                }            
                if( iteration % save_frequence == 0){
                    std::cout << " == Saving iteration " << iteration << " | " << iteration_percent << std::endl;
                    SaveSolution(t, iteration, vec_ptr_particles); 
                }
                ++iteration;    
                iteration_percent =  (iteration/iter_max)*100;
                t += dt;
            }
        }
        else{
            // Barnes and Hut algorithm O(NlogN)
            while ( t < Tmax ){
                std::cout << "t : " << t <<std::endl;
                QuadTree rootQt = QuadTree();
                rootQt.quadrant_center = Vect2(10.,10.);
                rootQt.quadrant_top_left = Vect2(-max_r, max_r);
                rootQt.quadrant_bottom_right = Vect2( max_r, -max_r);
                rootQt.size_side = 2*max_r;
                //std::cout << "Number of particles left : " << vec_ptr_particles.size() << std::endl;
                for ( auto it = vec_ptr_particles.begin(); it != vec_ptr_particles.end(); ++it){  
                    
                    std::shared_ptr<Particle> & particle= *it;
                    particle->computePositionAtHalfTimeStep(dt); // Drift every particle once
                    particle->resetAcceleration();

                    if((std::abs(particle->get_pos_x()) < max_r) && (std::abs(particle->get_pos_y()) < max_r) ){
                        if(particle->mass > 0.0){
                            rootQt.insert(particle);
                        }        
                    }
                    else{
                        std::cout << "Particle "<< *it << " removed"<<std::endl;
                        vec_ptr_particles.erase(it);
                    }
                }
                //prinf_info_quadtree(rootQt);
                #pragma omp parallel for
                for(auto & particle : vec_ptr_particles){
                    rootQt.computeAccelerationOn(particle, theta);
                    particle->computeVelocityAndPosition(dt);  // kick and drift
                }
                
                if( iteration % save_frequence == 0){
                    std::cout << " == Saving iteration " << iteration << " | " << iteration_percent << " %"<< std::endl;                 
                    SaveSolution(t, iteration, vec_ptr_particles); 
                }
                ++iteration;    
                iteration_percent =  (iteration/iter_max)*100;
                t += dt;
            }
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Duration of while loop : " <<duration.count() << " seconds"<<std::endl; 
    }
}




