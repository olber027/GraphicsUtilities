//
// Created by molberding on 6/28/2018.
//

#ifndef LINEARALGEBRA_LINEARALGEBRATESTDRIVER_H
#define LINEARALGEBRA_LINEARALGEBRATESTDRIVER_H

#include "TestDriver.h"
#include "../LinearAlgebra/Point.h"
#include "../LinearAlgebra/Vector.h"
#include "../LinearAlgebra/Matrix.h"
#include <vector>

using namespace linear_algebra;

class LinearAlgebraTestDriver : public TestDriver {
private:

public:

    LinearAlgebraTestDriver() {
        init("Linear Algebra");
    }

    void VectorConstructorTests() {
        Vector vec;
        assert(0, vec.getSize());
        vec = Vector(5);
        assert(5, vec.getSize());
        vec = Vector(5, 1);
        assert(5, vec.getSize());
        assert(1.0, vec[0]);
        vec = Vector(1.0, 2.0);
        assert(2, vec.getSize());
        assert(true, vec.is2d());
        assert(2.0, vec[1]);
        vec = Vector(1.0, 2.0, 3.0);
        assert(3, vec.getSize());
        assert(true, vec.is3d());
        assert(3.0, vec[2]);
        double temp[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        vec = Vector(temp, 6);
        assert(6, vec.getSize());
        std::vector<double> stdVec;
        for(int i = 0; i < 10; i++) {
            stdVec.push_back(i);
        }
        vec = Vector(stdVec);
        assert(10, vec.getSize());
        assert(9.0, vec[9]);
    }

    void PointConstructorTests() {
        Point point;
        assert(0, point.getSize());
        point = Point(5);
        assert(5, point.getSize());
        point = Point(5, 1);
        assert(5, point.getSize());
        assert(1.0, point[0]);
        point = Point(1.0, 2.0);
        assert(2, point.getSize());
        assert(true, point.is2d());
        assert(2.0, point[1]);
        point = Point(1.0, 2.0, 3.0);
        assert(3, point.getSize());
        assert(true, point.is3d());
        assert(3.0, point[2]);
        double temp[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        point = Point(temp, 6);
        assert(6, point.getSize());
        std::vector<double> stdVec;
        for(int i = 0; i < 10; i++) {
            stdVec.push_back(i);
        }
        point = Point(stdVec);
        assert(10, point.getSize());
        assert(9.0, point[9]);
    }

    void MatrixConstructorTests() {

    }

    void VectorOperatorTests() {
        Vector vec(1.0, 1.0);
        Vector result = vec * 2;
        assert(2.0, result[0]);
        assert(2.0, result[1]);
        result /= 2;
        assert(1.0, result[0]);
        assert(1.0, result[1]);
        result *= 2;
        assert(2.0, result[0]);
        assert(2.0, result[1]);
        result = result / 2;
        assert(1.0, result[0]);
        assert(1.0, result[1]);
        result = 2 * result;
        assert(2.0, result[0]);
        assert(2.0, result[1]);
        result = vec + vec;
        assert(2.0, result[0]);
        assert(2.0, result[1]);
        result -= vec;
        assert(1.0, result[0]);
        assert(1.0, result[1]);
        result += vec;
        assert(2.0, result[0]);
        assert(2.0, result[1]);
        result.subtract(vec);
        assert(1.0, result[0]);
        assert(1.0, result[1]);
        result.add(vec);
        assert(2.0, result[0]);
        assert(2.0, result[1]);

        vec = Vector(1.0, 0.0, 0.0);
        result = Vector(0.0, 1.0, 0.0);
        double temp = vec*result;
        assert(0.0, temp);
        temp = vec.dot(result);
        assert(0.0, temp);
        Vector crossVec = vec.cross(result);
        assert(Vector(0.0, 0.0, 1.0), crossVec);
    }

    void PointOperatorTests() {
        Point point(1.0, 1.0);
        Vector vec(2.0, 2.0);
        point *= 2;
        assert(2.0, point[0]);
        assert(2.0, point[1]);
        point /= 2;
        assert(1.0, point[0]);
        assert(1.0, point[1]);
        point = 2 * point;
        assert(2.0, point[0]);
        assert(2.0, point[1]);
        point = point / 2;
        assert(1.0, point[0]);
        assert(1.0, point[1]);
        point = point * 2;
        assert(2.0, point[0]);
        assert(2.0, point[1]);
        point -= vec;
        assert(0.0, point[0]);
        assert(0.0, point[1]);
        point += vec;
        assert(2.0, point[0]);
        assert(2.0, point[1]);
        point = point - vec;
        assert(0.0, point[0]);
        assert(0.0, point[1]);
        point = point + vec;
        assert(2.0, point[0]);
        assert(2.0, point[1]);
        point.subtract(vec);
        assert(0.0, point[0]);
        assert(0.0, point[1]);
        point.add(vec);
        assert(2.0, point[0]);
        assert(2.0, point[1]);
    }

    void MatrixOperatorTests() {

    }



    void run() {
        VectorConstructorTests();
        VectorOperatorTests();
        PointConstructorTests();
        PointOperatorTests();
    }
};

#endif //LINEARALGEBRA_LINEARALGEBRATESTDRIVER_H
