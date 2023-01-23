#ifndef PARTICLE_H
#define PARTICLE_H

#include "vect2.hpp"
#include <memory>

extern float G;
extern float epsilon;
class Particle {
    public :
        Vect2 position;
        Vect2 velocity;
        Vect2 acceleration;
        float mass;
        //float Ek = 0.5*mass*( (velocity.norm2()) * (velocity.norm2()) ); // kinetic energy
        //float U; // potential energy
        


        Particle();
        Particle(const Vect2 &, const Vect2 &, const Vect2 &, float);
        Particle(const Particle &);

        void computePositionAtHalfTimeStep(float dt);
        void computeVelocityAndPosition(float dt);
        void computeAccelerationDueTo(const std::shared_ptr<Particle> & );
        void resetAcceleration();
        float get_pos_x() const { return position.axis_x;};
        float get_pos_y() const { return position.axis_y;};
        float get_vel_x() const { return velocity.axis_x;};
        float get_vel_y() const { return velocity.axis_y;};
        float get_acc_x() const { return acceleration.axis_x;};
        float get_acc_y() const { return acceleration.axis_y;};
        float get_velocity_value(const Particle &);
        void set_position(float , float );
        void set_velocity(float , float );
        void set_acceleration(float , float );

};

#endif


