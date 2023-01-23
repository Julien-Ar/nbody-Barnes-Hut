#ifndef VECT2_H
#define VECT2_H



class Vect2{
    public :
        // struct to hold 2 elements along a x-axis and y-axis,
        // and overloading some operators for ease of use
        double axis_x;
        double axis_y;

        Vect2();
        Vect2(const Vect2&);
        Vect2(double, double);
        
        void print();
        // operators that create new object
        Vect2 operator+(const Vect2&);
        Vect2 operator-(const Vect2&);
        Vect2 operator*(double);
        Vect2 operator*(const Vect2 &);
        double norm2();

};


#endif
