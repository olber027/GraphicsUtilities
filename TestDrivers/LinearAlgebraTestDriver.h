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
#include <random>
#include <ctime>

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
        vec = Vector(5, 1.0);
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
        point = Point(5, 1.0);
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
        Matrix mat(2, 2);
        assert(2, mat.numCols());
        assert(2, mat.numRows());
        assert(0.0, mat[0][0]);
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
        vec[0] = 75.0;
        assert(75.0, vec[0]);
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
        Matrix ident = Matrix::identity(2, 2);
        Matrix other(2,2);
        other[0][0] = 1.0;
        other[0][1] = 2.0;
        other[1][0] = 3.0;
        other[1][1] = 4.0;
        Matrix result = ident * other;
        assert(1.0, result[0][0]);
        assert(2.0, result[0][1]);
        assert(3.0, result[1][0]);
        assert(4.0, result[1][1]);
        Vector vec(2, 5.0);
        Vector resultVec = ident * vec;
        assert(5.0, vec[0]);
        assert(5.0, vec[1]);
        Point point(2, 5.0);
        Point resultPoint = ident * point;
        assert(5.0, point[0]);
        assert(5.0, point[1]);
        Vector row0, row1, row2, row3;
        row0.pushBack(5.0, -2.0, 2.0, 7.0);
        row1.pushBack(1.0, 0.0, 0.0, 3.0);
        row2.pushBack(-3.0, 1.0, 5.0, 0.0);
        row3.pushBack(3.0, -1.0, -9.0, 4.0);
        Matrix mat4(4,4);
        mat4.setRow(0, row0);
        mat4.setRow(1, row1);
        mat4.setRow(2, row2);
        mat4.setRow(3, row3);
        Matrix inv = mat4.getInverse();
        Matrix mat = (mat4*inv);
        assert(Matrix::identity(4,4), mat);
    }

    void run() {
        VectorConstructorTests();
        VectorOperatorTests();
        PointConstructorTests();
        PointOperatorTests();
        MatrixConstructorTests();
        MatrixOperatorTests();
    }
};

#endif //LINEARALGEBRA_LINEARALGEBRATESTDRIVER_H
