add_executable(main main.cpp particle.cpp vect2.cpp quadtree.cpp quadtree.hpp utils.hpp particle.hpp) 
find_package(OpenMP)
#target_link_libraries(main fmt) 
if(OpenMP_CXX_FOUND)
    target_link_libraries(main fmt OpenMP::OpenMP_CXX)
endif()
