//
// Created by molberding on 6/28/2018.
//

#ifndef LINEARALGEBRA_TESTDRIVER_H
#define LINEARALGEBRA_TESTDRIVER_H

#include <string>
#include <sstream>

class TestDriver {
protected:
    int passed;
    int failed;
    int total;
    std::stringstream output;
public:
    virtual void run() = 0;

    virtual void init(std::string driverName) {
        passed = 0;
        failed = 0;
        total = 0;

        output = std::stringstream();
        for(int i = 0; i < 80; i++) {
            output << "=";
        }
        output << std::endl << "Start of " << driverName << " TestDrivers\n";
        for(int i = 0; i < 80; i++) {
            output << "=";
        }
        output << std::endl;
    }

    virtual std::string getOutput() {
        output << "Test Run Results were:";
        output << "\n\tpassed: " << passed;
        output << "\n\tfailed: " << failed;
        output << "\n\ttotal: " << total << std::endl;
        return output.str();
    }

    virtual int getFailures() {
        return failed;
    }

    virtual int getPasses() {
        return passed;
    }

    template<typename T> bool assert(T expected, T actual) {
        bool result = expected == actual;
        if(result) {
            passed++;
        } else {
            output << "TestFailed\n__________\n";
            output << "\texpected: " << expected << "\n\tactual  : " << actual << std::endl;
            failed++;
        }
        total++;
        return result;
    }

    template<typename T> bool assert(T expected, T actual, const char* message) {
        bool result = expected == actual;
        if(result) {
            passed++;
        } else {
            output << "TestFailed\n__________\n";
            output << "\texpected: " << expected << "\n\tactual  : " << actual << std::endl;
            output << "\tmessage : " << message << std::endl;
            failed++;
        }
        total++;
        return result;
    }
};

#endif //LINEARALGEBRA_TESTDRIVER_H
