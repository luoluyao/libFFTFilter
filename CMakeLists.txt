cmake_minimum_required (VERSION 2.6) 
project (FFTFilter) 

include_directories (${FFTFilter_SOURCE_DIR}/include  
                     ${FFTFilter_SOURCE_DIR}/kissfft) 

ADD_SUBDIRECTORY (kissfft)

add_library (${PROJECT_NAME} src/fft_filter.cpp)
target_link_libraries (${PROJECT_NAME} kissfft)  

ADD_SUBDIRECTORY (test)