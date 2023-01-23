#include <math.h>
#include <iostream>
#include "particle.hpp"



Particle::Particle(){
    position = Vect2(0,0);
    velocity = Vect2(0,0);
    acceleration = Vect2(0,0);
    mass = 0;
    //Ek = 0;
    //U = 0;
};

Particle::Particle(const Vect2 & pos, const Vect2 & vel, 
                    const Vect2  & accel, float mass1){
  position = pos;
  velocity = vel;
  acceleration = accel;
  mass = mass1;
};

Particle::Particle(const Particle & p1){
    position = p1.position;
    velocity = p1.velocity;
    acceleration = p1.acceleration;
    mass = p1.mass;
    //Ek = p1.Ek;
    //U = p1.U;
};


void Particle::computeAccelerationDueTo(const std::shared_ptr<Particle>  & p1){
    // We compute the acceleration on the current particle object
    // due to Particle p1
    
    float dx = get_pos_x() - p1->get_pos_x();
    float dy = get_pos_y() - p1->get_pos_y();
    float distance = pow(dx*dx + dy*dy, 0.5);
    float Fx;
    float Fy;
    if(distance < epsilon){
        Fx = -(G*(p1->mass)*dx)/pow((distance*distance + epsilon*epsilon), 1.5); 
        Fy = -(G*(p1->mass)*dy)/pow((distance*distance + epsilon*epsilon), 1.5); 
    }
    else{        
        Fx = -(G*(p1->mass)*dx)/pow(distance , 3.); // Force intensity
        Fy = -(G*(p1->mass)*dy)/pow(distance , 3.); // Force intensity
    }
    

    Vect2 acceleration_change( Fx, Fy );

    acceleration = acceleration + acceleration_change;
    //U += G*(p1.mass)*mass/distance;
    
};


void Particle::computePositionAtHalfTimeStep(float dt){
    
    position = position + velocity * (dt/2);

};


void Particle::computeVelocityAndPosition(float dt){
    // Leap-Frog integration scheme
    
    //DKD Leap Frog 
    velocity = velocity + acceleration*(dt); // Kick
    position = position + velocity*(dt/2); // Drift 


    
};


void Particle::resetAcceleration(){
    //set_acceleration(0, 0);
    acceleration.axis_x = 0.;
    acceleration.axis_y = 0.;
    //U = 0;
};



float Particle::get_velocity_value(const Particle &P)  {
    float vx = P.get_vel_x();
    float vy = P.get_vel_y();
    float v = sqrt(vx*vx + vy*vy);
    return v;  
}


void Particle::set_position(float x, float y){
    position.axis_x = x;
    position.axis_y = y;
}

void Particle::set_velocity(float vel_x, float vel_y){
    velocity.axis_x = vel_x;
    velocity.axis_y = vel_y;
}

void Particle::set_acceleration(float acc_x, float acc_y){
    acceleration.axis_x = acc_x;
    acceleration.axis_y = acc_y;
}
