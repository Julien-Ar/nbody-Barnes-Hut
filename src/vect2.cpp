#include <iostream>
#include <math.h>
#include "vect2.hpp"

Vect2::Vect2(){
    axis_x = 0;
    axis_y = 0;
};

Vect2::Vect2(const Vect2 & other_vect2){
    axis_x = other_vect2.axis_x;
    axis_y = other_vect2.axis_y;
};

Vect2::Vect2(double x, double y){
    axis_x = x;
    axis_y = y;
};


void Vect2::print(){
    std::cout << "axis_x : "<< axis_x <<"\n"<< "axis_y : " << axis_y << std::endl;
}



Vect2 Vect2::operator+(const Vect2 & other_vect2 ){
    Vect2 result;
    result.axis_x = axis_x + other_vect2.axis_x;
    result.axis_y = axis_y + other_vect2.axis_y;
    return result; 
}   

Vect2 Vect2::operator-(const Vect2 & other_vect2 ){
    Vect2 result;
    result.axis_x = axis_x - other_vect2.axis_x;
    result.axis_y = axis_y - other_vect2.axis_y;
    return result; 
}   


Vect2 Vect2::operator*(double t){
    Vect2 result;
    result.axis_x = axis_x*t;
    result.axis_y = axis_y*t;
    return result; 
}   


Vect2 Vect2::operator*(const Vect2 & other_vect2){
    Vect2 result;
    result.axis_x = axis_x*other_vect2.axis_x;
    result.axis_y = axis_y*other_vect2.axis_y;
    return result;
}

double Vect2::norm2(){
    return sqrt(axis_x*axis_x +axis_y*axis_y);
}