cmake_minimum_required(VERSION 3.10)

project(LinearAlgebra)
set(CMAKE_CXX_STANDARD 11)
add_executable(LinearAlebgra Testing/LinearAlgebraTests.cpp LinearAlgebra/Vector.h LinearAlgebra/Point.h LinearAlgebra/Matrix.h TestDrivers/TestDriver.h TestDrivers/LinearAlgebraTestDriver.h)