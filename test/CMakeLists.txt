
enable_testing()
find_package(GTest REQUIRED)
find_package (Threads)

include_directories (${GTEST_INCLUDE_DIR}) 

add_executable(test_fft_filter test_fft_filter.cpp)
target_link_libraries(test_fft_filter ${PROJECT_NAME} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

add_test(
    NAME test_fft_filter
    COMMAND test_fft_filter
)

add_executable(test_fft test_fft.cpp)
target_link_libraries(test_fft kissfft ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

add_test(
    NAME test_fft
    COMMAND test_fft
)