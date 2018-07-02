#include <iostream>
#include "../TestDrivers/LinearAlgebraTestDriver.h"

using namespace linear_algebra;

int main() {

    LinearAlgebraTestDriver* driver = new LinearAlgebraTestDriver();
    driver->run();

    std::cout << driver->getOutput() << std::endl;

    delete driver;

    return 0;
}