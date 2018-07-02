//
// Created by molberding on 6/25/2018.
//

#ifndef GRAPHICSUTILITIES_MATRIX_H
#define GRAPHICSUTILITIES_MATRIX_H

#include "Vector.h"
#include "Point.h"

namespace linear_algebra {
    class Matrix {
    private:
        int rows, cols;
        Vector* vals;

        void destroy() {
            if(vals != nullptr) {
                delete [] vals;
                rows = 0;
                cols = 0;
            }
        }

        void initialize(const int rowSize, const int colSize) {
            destroy();
            rows = rowSize;
            cols = colSize;
            vals = new Vector[rows];
            for(int i = 0; i < rows; i++) {
                vals[i] = Vector(cols, 0);
            }
        }

    public:
        Matrix() : rows(0), cols(0), vals(nullptr) {}
        Matrix(const int rowSize, const int colSize) : Matrix() {
            initialize(rowSize, colSize);
        }
        Matrix(const int rowSize, const int colSize, const double fill) : Matrix(rowSize, colSize) {
            for(int i = 0; i < rows; i++) {
                vals[i] = Vector(cols, fill);
            }
        }
        Matrix(const int rowSize, const int colSize, double** source) : Matrix(rowSize, colSize) {
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    vals[i][j] = source[i][j];
                }
            }
        }
        Matrix(const Matrix& other) {
            if(this != &other) {
                vals = new Vector[rows];
                for(int i = 0; i < rows; i++) {
                    vals[i] = other[i];
                }
            }
        }
        Matrix(const Matrix&& other) noexcept {
            rows = other.rows;
            cols = other.cols;
            vals = other.vals;
        }

        Vector& operator[](const int rowIndex) const {
            return vals[rowIndex];
        }

    };
}

#endif //GRAPHICSUTILITIES_MATRIX_H
