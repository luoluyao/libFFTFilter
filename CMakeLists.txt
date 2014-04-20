cmake_minimum_required (VERSION 2.6) 
project (FFTFilter) 

include_directories (${FFTFilter_SOURCE_DIR}/include) 

add_library( kissfft src/kiss_fft.c)

add_library (fft_filter src/fft_filter.cpp)
target_link_libraries (fft_filter kissfft)  

