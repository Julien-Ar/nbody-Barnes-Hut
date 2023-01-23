#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <string>
#include <sstream>
#include <ostream>
#include <iomanip>
#include <tgmath.h>
#include <random>
#include <math.h>
#include <fmt/core.h>
#include <filesystem>

#include "quadtree.hpp"
#include "particle.hpp"
#include "vect2.hpp"

// Defines several initial conditions models along some
// other functions : - one that read a file containing initial conditions
//                   - another that save snapchots of each particle positions in a file 
//                     at a given timestep
//                   - another that prints out info about the quadtree


std::vector< std::shared_ptr<Particle> > create_initial_conditions_cold_start(float M, 
                                                                             float N_mass,
                                                                             float N_tracer,
                                                                             float R,
                                                                             int seed = 4321){
    // Create cold start initial conditions :   position randomly sampled inside a 
    //                                          circle of radii R and no velocity
    std::vector< std::shared_ptr<Particle> > vec_ptr_particles;

    std::string current_path = std::filesystem::current_path();
    std::string file_name;
    std::string core_file_name = fmt::format("coldStart_{:g}M_{:g}Nm_{:g}Nt_{:g}R_{:g}seed", 
                                                        M,  N_mass, N_tracer, R, seed);

    int count = 0; // count the number of time a file with the same current initial conditions has been generated
    // for loop to check for same file
    for(const auto & file_tmp : std::filesystem::directory_iterator(current_path)){
        std::string str_file_tmp(file_tmp.path());
        if (str_file_tmp.find(core_file_name) != std::string::npos) {
            ++count;
        }
    }
    if(count == 0) {
        file_name = core_file_name + ".txt";
    }
    else{
        file_name = core_file_name + "_" + std::to_string(count) + ".txt";
    }
    std::ofstream file;
    file.open(file_name);
    file << "NbMass\n";
    file << N_mass+1 << "\n";
    float SMBH_mass = 10*M;
    file << 0. << " " << 0. << " " << 0. << " " << 0. << " "<< SMBH_mass << "\n";
    std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(0., 0.) ) ;
    std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
    std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
    std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle( *vect2_pos_ptr, 
                                                                                    *vect2_vel_ptr,
                                                                                    *vect2_acc_ptr, 
                                                                                    SMBH_mass));
    vec_ptr_particles.push_back(particle_ptr);

    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    generator.seed(seed);
    std::uniform_real_distribution<float>  unif(0., 1.);
    std::uniform_real_distribution<float>  unif_angle(0., 2*M_PI); // sample a float from 0 to pi
    
    float mi = M/N_mass; // mass of individual star
    for(int i = 0 ; i < N_mass; ++i ){
        float theta = unif_angle(generator);
        float U = unif(generator);
        float r = U*R;
        float x = r*cos(theta);
        float y = r*sin(theta);
        std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(x, y) ) ;
        std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(0.0, 0.0) );
        std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0.0, .0) );
        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                            *vect2_acc_ptr, mi));
        file << x << " " << y << " " << 0. << " " << 0. << " "<< mi << "\n";

        vec_ptr_particles.push_back(particle_ptr);
    }
    for(int i = 0 ; i < N_tracer; ++i ){
        float theta = unif_angle(generator);
        float U = unif(generator);
        float r = U*R;
        float x = r*cos(theta);
        float y = r*sin(theta);
        std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(x, y) ) ;
        std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(0.0, 0.0) );
        std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0.0, 0.0) );
        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                            *vect2_acc_ptr, 0.0));
        file << x << " " << y << " " << 0. << " " << 0. << " "<< 0.0<< "\n";

        vec_ptr_particles.push_back(particle_ptr);
    } 
    file.close(); 
    return(vec_ptr_particles);
}


std::vector< std::shared_ptr<Particle> > create_initial_conditions_kepler(float M, 
                                                                             float N_mass,
                                                                             float N_tracer,
                                                                             float R,
                                                                             int seed = 4321){
    // Create cold start initial conditions :   position randomly sampled inside a 
    //                                          circle of radii R and orbital velocity
        // TO DO : COMPUTE ORBITAL VELOCITY
    
    std::vector< std::shared_ptr<Particle> > vec_ptr_particles;

    std::string current_path = std::filesystem::current_path();
    std::string file_name;
    std::string core_file_name = fmt::format("kepler_{:g}G_{:g}M_{:g}Nm_{:g}Nt_{:g}R_{:g}seed", 
                                                        G, M,  N_mass, N_tracer, R, seed);

    int count = 0; // count the number of time a file with the same current initial conditions has been generated
    // for loop to check for same file
    for(const auto & file_tmp : std::filesystem::directory_iterator(current_path)){
        std::string str_file_tmp(file_tmp.path());
        if (str_file_tmp.find(core_file_name) != std::string::npos) {
            ++count;
        }
    }
    if(count == 0) {
        file_name = core_file_name + ".txt";
    }
    else{
        file_name = core_file_name + "_" + std::to_string(count) + ".txt";
    }
    std::ofstream file;
    file.open(file_name);
    file << "NbMass\n";
    file << N_mass+1 << "\n";
    float SMBH_mass = 10*M;
    file << 0. << " " << 0. << " " << 0. << " " << 0. << " "<< SMBH_mass << "\n";
    std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(0., 0.) ) ;
    std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
    std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
    std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle( *vect2_pos_ptr, 
                                                                                    *vect2_vel_ptr,
                                                                                    *vect2_acc_ptr, 
                                                                                    SMBH_mass));
    vec_ptr_particles.push_back(particle_ptr);

    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    generator.seed(seed);
    std::uniform_real_distribution<float>  unif(0., 1.);
    std::uniform_real_distribution<float>  unif_angle(0., 2*M_PI); // sample a float from 0 to pi
    
    float mi = M/N_mass; // mass of individual star
    


    for(int i = 0 ; i < N_mass; ++i ){
        float theta = unif_angle(generator);
        float U = unif(generator);
        float r = U*R;
        float x = r*cos(theta);
        float y = r*sin(theta);
        std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(x, y) ) ;
        std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(0.0, 0.0) );
        std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0.0, .0) );
        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                            *vect2_acc_ptr, mi));
        file << x << " " << y << " " << 0 << " " << 0 <<" "<<mi << "\n";  

        vec_ptr_particles.push_back(particle_ptr);
    }
    for(int i = 0 ; i < N_tracer; ++i ){
        float theta = unif_angle(generator);
        float U = unif(generator);
        float r = U*R;
        float x = r*cos(theta);
        float y = r*sin(theta);
        std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(x, y) ) ;
        std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(0.0, 0.0) );
        std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0.0, 0.0) );
        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                            *vect2_acc_ptr, 0.0));
        file << x << " " << y << " " << 0 << " " << 0 <<" "<< 0.0 << "\n";  
        vec_ptr_particles.push_back(particle_ptr);
    }  
    file.close();
    return(vec_ptr_particles);
}


float g_q(float q){
    return( pow(q, 2)*pow(1 - pow(q, 2), 3.5) );
}

std::vector< std::shared_ptr<Particle> > create_initial_conditions_plummer( float a, 
                                                                    float M, 
                                                                    float N_mass, 
                                                                    float N_tracer =0, 
                                                                    float a_multiplier_tracer = 1,
                                                                    int seed = 1234,
                                                                    bool SMBH=false){
    // create initial conditions using the Plummer model and store them in a file
    // a : scale factor for the core of the galaxy, related to the particles with masses
    // M : total mass of the galaxy
    // N_mass : number of body with masses
    // N_tracer : number of body without masses 
    // a_multiplier_tracer : by how much to multiply "a", in order to scatter away from the core the particle
    //                       without masses or group them near it
    std::cout << "Creating initial conditions..." << std::endl;
    float bool_SMBH_to_float = static_cast<float>(SMBH);
    float mi = M/N_mass;
    if((M <=1) && (G <= 1) ){
        // henon units
        a = 3*M_PI/16.;
    } 
    float X1, r, theta, x, y, Ve, X4, X5, V, U, vx, vy, SMBH_mass;
    std::vector< std::shared_ptr<Particle> > vec_ptr_particles;

    std::string current_path = std::filesystem::current_path();
    std::string file_name;
    std::string core_file_name = fmt::format("plummer_{:g}G_{:g}a_{:g}M_{:g}Nm_{:g}Nt_{:g}amult_{:g}seed_{:g}SBMH", 
                                                        G, a, M, N_mass, N_tracer, a_multiplier_tracer, float(seed), bool_SMBH_to_float);

    int count = 0; // count the number of time a file with the same current initial conditions has been generated
    // for loop to check for same file
    for(const auto & file_tmp : std::filesystem::directory_iterator(current_path)){
        std::string str_file_tmp(file_tmp.path());
        if (str_file_tmp.find(core_file_name) != std::string::npos) {
            ++count;
        }
    }
    if(count == 0) {
        file_name = core_file_name + ".txt";
    }
    else{
        file_name = core_file_name + "_" + std::to_string(count) + ".txt";
    }

    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    generator.seed(seed);
    std::uniform_real_distribution<float>  unif(0., 1.);
    std::uniform_real_distribution<float>  unif_angle(0., 2*M_PI); // sample a float from 0 to pi

    std::ofstream file;
    file.open(file_name);
    file << "NbMass\n";
    
    if(SMBH){
        std::cout << "SMBH created" << std::endl;
        SMBH_mass = 0.001*M;
        file << std::scientific << N_mass+1 << "\n"; 
        file << 0. << " " << 0. << " " << 0. << " " << 0. << " "<< SMBH_mass << "\n";
        std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(0., 0.) ) ;
        std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
        std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle( *vect2_pos_ptr, 
                                                                                        *vect2_vel_ptr,
                                                                                        *vect2_acc_ptr, 
                                                                                        SMBH_mass));
        vec_ptr_particles.push_back(particle_ptr);

    }
    else{
        file << std::scientific << N_mass << "\n"; 
    }
    for(int i = 0; i<N_mass; ++i){
        X1 = unif(generator);
        while (X1 > 0.99) { X1 = unif(generator); }
        r = a*pow(pow(X1, -2./3)-1, -0.5);
        theta = unif_angle(generator);
        x = r*cos(theta); // position
        y = r*sin(theta); // position
        if(SMBH){
            Ve = sqrt(2*G*(M+SMBH_mass)/a)*pow(1 + pow(r/a, 2), -0.25); // escape velocity
        }
        else{
            Ve = sqrt(2*G*M/a)*pow(1 + pow(r/a, 2), -0.25); // escape velocity
        }
        X4 = unif(generator);
        X5 = unif(generator);
        while( (0.1*X5) > g_q(X4) ) {
            X4 = unif(generator);
            X5 = unif(generator);      
        }
        //bool X4_X5 = ((0.1*X5) > g_q(X4));
        //std::cout <<  (0.1*X5) << " " << g_q(X4)<<" " <<X4_X5 << std::endl;
        V = X4*Ve;
        U = unif(generator);
        if( U < 0.5 ) { V = -V ;}
        
        theta = unif_angle(generator);
        vx = V*cos(theta); // velocity
        vy = V*sin(theta); // velocity
        file << x << " " << y << " " << vx << " " << vy << " " << mi<< "\n";

        std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(x, y) ) ;
        std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(vx, vy) );
        std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                            *vect2_acc_ptr, mi));
        vec_ptr_particles.push_back(particle_ptr);
    }
    for(int i = 0; i<N_tracer; ++i){
        X1 = unif(generator);
        while (X1 > 0.99) { X1 = unif(generator); }
        r = (a_multiplier_tracer*a)*pow(pow(X1, -2./3)-1, -0.5); // extend the radius for tracer particles
        theta = unif_angle(generator);
        x = r*cos(theta); // position
        y = r*sin(theta); // position
        if(SMBH){
            Ve = sqrt(2*G*(M+SMBH_mass)/(a*a_multiplier_tracer))*pow(1 + pow(r/(a*a_multiplier_tracer), 2), -0.25); // escape velocity
        }
        else{
            Ve = sqrt(2*G*M/a)*pow(1 + pow(r/a, 2), -0.25); // escape velocity
        }
        X4 = unif(generator);
        X5 = unif(generator);
        //std::cout <<  (0.1*X5) << " " << g_q(X4) << std::endl;
        while( (0.1*X5) > g_q(X4) ) {
            //std::cout <<  (0.1*X5) << " " << g_q(X4) << std::endl;
            X4 = unif(generator);
            X5 = unif(generator);      
        }
        V = X4*Ve;
        U = unif(generator);
        if( U < 0.5 ) { V *= -1 ;}
        theta = unif_angle(generator);
        vx = V*cos(theta); // velocity
        vy = V*sin(theta); // velocity
        file << x << " " << y << " " << vx << " " << vy <<" " <<0.0<<"\n";  
        std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(x, y) ) ;
        std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(vx, vy) );
        std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );
        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                            *vect2_acc_ptr, 0.));
        vec_ptr_particles.push_back(particle_ptr);

    }
    file.close();
    std::cout << "Initial conditions created" << std::endl;

    return(vec_ptr_particles);
}





std::vector< std::shared_ptr<Particle> > read_IC_from_file(std::string file_name){  

    // Read initial conidtions file written as : 
    // G a M NbMass NbTracer
    // 0.043 100 10**7 1000 10
    // x y vx vy    the first NbMass lines with positions and velocities corresponds to body with masses
    // . . .. .. 
    std::cout << "Readind initial conditions from file : "<<file_name << std::endl;

    std::vector< std::shared_ptr<Particle> > vec_ptr_particles;

    float x, y, vx, vy, trash;
    float mi; 
    int NbMass; // Nbmass : number of particles with masses = mi
    std::ifstream file;
    file.open(file_name);
    std::string line;
    bool first_line = true;
    bool second_line = true;
    while ( file.is_open()){
        if (file.eof()){
            file.close();
        }
        else{
            int i = 0; // line position
            while( std::getline(file, line)){
                if(first_line) { 
                    // skip first line
                    first_line = false; 
                } 
                else if (second_line) {
                    second_line = false;
                    std::istringstream sString( line );
                    sString >> NbMass;
                    ;
                } 
                else {
                    std::istringstream sString(line);
                    sString >> x >> y >> vx >> vy >>mi;

                    std::shared_ptr<Vect2> vect2_pos_ptr = std::make_shared<Vect2>( Vect2(x, y) ) ;
                    std::shared_ptr<Vect2> vect2_vel_ptr = std::make_shared<Vect2>( Vect2(vx, vy) );
                    std::shared_ptr<Vect2> vect2_acc_ptr = std::make_shared<Vect2>( Vect2(0., 0.) );

                    if (i+1 <= NbMass ){
                        //std::cout << "Read line with particle with masses | i = " << i << std::endl;
                        // first block of data corresponds to body with the same mass
                        std::shared_ptr<Particle> particle_ptr = std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                                         *vect2_acc_ptr, mi));
                        vec_ptr_particles.push_back(particle_ptr);
                    }
                    else{
                        //std::cout << "Read line with particle without masses | i = " << i << std::endl;

                        // second block of data corresponds to tracer particles, bodies without masses
                        std::shared_ptr<Particle> particle_ptr= std::make_shared<Particle> ( Particle(*vect2_pos_ptr, *vect2_vel_ptr,
                                                                         *vect2_acc_ptr, 0.));
                        vec_ptr_particles.push_back(particle_ptr);
                    }  
                    ++i;
                }
            }

        }
        
    }
    std::cout << "Initial conditions read." << std::endl;

    return(vec_ptr_particles);
}   


void SaveSolution(float t, int iteration, const std::vector<std::shared_ptr<Particle>>& vec_particles_ptr ){
    // Each iteration, we save a snapchot of every particle's position
    // in a file.
    
    std::ostringstream oString;
    oString << "gravity_" << std::setfill('0') << std::setw(7) << iteration << ".particles"; 
       
    std::ofstream file_out;
    file_out.open(oString.str());
    file_out << "# X Y Z" << std::endl;
    for (auto & p_mass : vec_particles_ptr){
        file_out << p_mass->get_pos_x() << " " << p_mass->get_pos_y() << " " << 0  << std::endl;
    };
 
     file_out.close();
}




void prinf_info_quadtree(const QuadTree & root){
    std::cout<<" ----- \n";
    std::cout << "Mass barycenter : ("<<root.mass_barycenter.axis_x
              <<","<<root.mass_barycenter.axis_y
              <<" | Nb part : "<<root.nb_particles<<" | Total mass : "<< root.total_mass 
              << " | size vec part : "<<root.vec_particle_ptr.size()<<"\n";
    if(root.isLeaf){
          std::cout << "isLeaf : part pos : ("<<root.vec_particle_ptr[0]->get_pos_x()
              <<","<<root.vec_particle_ptr[0]->get_pos_y()
              << " size vec part : " << root.vec_particle_ptr.size()<<"\n";  
    }
    auto children = std::vector<std::shared_ptr<QuadTree>>{root.TopLeft_subtree,
                                                           root.TopRight_subtree,
                                                           root.BottomLeft_subtree,
                                                           root.BottomRight_subtree};
    int i = 0;
    for(auto & child : children){
        if(i==0){
            std::cout << "\n SubTopLeft :\n";
            if(child != nullptr){
                prinf_info_quadtree(*child);
            }
            else{
                std::cout <<" Subtopleft nullptr\n";
            }
        }
        if(i==1){
            std::cout << "\n SubTopRight :\n";
            if(child != nullptr){
                prinf_info_quadtree(*child);
            }
            else{
                std::cout <<" SubTopRight nullptr\n";
            }
        }
        if(i==2){
            std::cout << "\n SubBottomLeft :\n";
            if(child != nullptr){
                prinf_info_quadtree(*child);
            }
            else{
                std::cout <<" SubBottomLeft nullptr\n";
            }
        }
        if(i==3){
            std::cout << "\n SubBottomRight :\n";
            if(child != nullptr){
                prinf_info_quadtree(*child);
            }
            else{
                std::cout <<" SubBottomRight nullptr\n";
            }
        }

        ++i;
    }


}






#endif