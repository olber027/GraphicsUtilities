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
                vals[i] = Vector(cols, 0.0);
            }
        }

        void checkBounds(const int index) const {
            if(index < 0 || index >= rows) {
                std::stringstream str;
                str << "index out of range: " << index;
                throw std::out_of_range(str.str());
            }
        }

        template<typename... Args>
        void checkBounds(const int index, Args... args) const {
            if(index < 0 || index >= rows) {
                std::stringstream str;
                str << "index out of range: " << index;
                throw std::out_of_range(str.str());
            }
            checkBounds(args...);
        }

        static void checkDimensionCompatibility(const Matrix& mat, const Vector& vec) {
            if(mat.numCols() != vec.numRows()) {
                throwDimensionError(mat, vec);
            }
        }
        static void checkDimensionCompatibility(const Matrix& mat, const Point& point) {
            if(mat.numCols() != point.numRows()) {
                throwDimensionError(mat, point);
            }
        }
        static void checkDimensionCompatibility(const Matrix& a, const Matrix& b) {
            if(a.numCols() != b.numRows()) {
                throwDimensionError(a, b);
            }
        }
        static void checkDimensionCompatibility(const Vector& vec, const Matrix& mat) {
            if(vec.numCols() != mat.numRows()) {
                throwDimensionError(mat, vec);
            }
        }
        static void checkDimensionCompatibility(const Point& point, const Matrix& mat) {
            if(point.numCols() != mat.numRows()) {
                throwDimensionError(mat, point);
            }
        }
        static void checkDimensionEquality(const Matrix& a, const Matrix& b) {
            if(a.rows != b.rows || a.cols != b.cols) {
                throwDimensionError(a, b);
            }
        }
        static void throwDimensionError(const Matrix& mat, const Vector& vec) {
            std::stringstream str;
            str << "The given vector was not compatible with the given matrix." << std::endl;
            str << "The dimensions were: " << std::endl;
            str << "Matrix - " << mat.numRows() << "x" << mat.numCols() << std::endl;
            str << "Vector - " << vec.numRows() << "x" << vec.numCols() << std::endl;
            throw std::runtime_error(str.str());
        }
        static void throwDimensionError(const Matrix& mat, const Point& point) {
            std::stringstream str;
            str << "The given point was not compatible with the given matrix." << std::endl;
            str << "The dimensions were: " << std::endl;
            str << "Matrix - " << mat.numRows() << "x" << mat.numCols() << std::endl;
            str << "Point  - " << point.numRows() << "x" << point.numCols() << std::endl;
            throw std::runtime_error(str.str());
        }
        static void throwDimensionError(const Matrix& a, const Matrix& b) {
            std::stringstream str;
            str << "The two matrices were not compatibile." << std::endl;
            str << "The given dimensions were: " << std::endl;
            str << "Matrix A - " << a.numRows() << "x" << a.numCols() << std::endl;
            str << "Matrix B - " << b.numRows() << "x" << b.numCols() << std::endl;
            throw std::runtime_error(str.str());
        }

        Matrix getMinorMatrix(const int rowIndex, const int columnIndex) const {
            Matrix result(rows - 1, cols - 1);
            int resultI = 0;
            for(int i = 0; i < rows; i++) {
                if(i == rowIndex) {
                    continue;
                }
                int resultJ = 0;
                for(int j = 0; j < cols; j++) {
                    if(j == columnIndex) {
                        continue;
                    }
                    result[resultI][resultJ] = vals[i][j];
                    resultJ++;
                }
                resultI++;
            }
            return result;
        }

    public:
        Matrix() : rows(0), cols(0), vals(nullptr) {}
        Matrix(const int rowSize, const int colSize, const double fill) {
            rows = rowSize;
            cols = colSize;
            vals = new Vector[rows];
            for(int i = 0; i < rows; i++) {
                vals[i] = Vector(cols, fill);
            }
        }
        Matrix(const int rowSize, const int colSize) : Matrix(rowSize, colSize, 0.0) { }
        Matrix(const int rowSize, const int colSize, double** source) : Matrix(rowSize, colSize, 0.0) {
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    vals[i][j] = source[i][j];
                }
            }
        }
        Matrix(const Matrix& other) : Matrix() {
            rows = other.rows;
            cols = other.cols;
            vals = new Vector[rows];
            for(int i = 0; i < rows; i++) {
                vals[i] = other[i];
            }
        }
        Matrix(const Matrix&& other) noexcept {
            rows = other.rows;
            cols = other.cols;
            vals = other.vals;
        }
        Matrix& operator=(const Matrix& other) {
            if(this != &other) {
                destroy();
                rows = other.rows;
                cols = other.cols;
                vals = new Vector[rows];
                for(int i = 0; i < rows; i++) {
                    vals[i] = other[i];
                }
            }
            return *this;
        }

        friend std::ostream& operator<<(std::ostream& out, const Matrix& mat) {
            for(int i = 0; i < mat.rows; i++) {
                out << mat.vals[i] << std::endl;
            }
            return out;
        }

        // Matrix Accessors
        Vector& operator[](const int rowIndex) const {
            checkBounds(rowIndex);
            return vals[rowIndex];
        }
        Vector& getRow(const int rowIndex) const {
            checkBounds(rowIndex);
            return vals[rowIndex];
        }
        Vector  getCol(const int columnIndex) const {
            Vector result(rows);
            for(int i = 0; i < rows; i++) {
                result[i] = vals[i][columnIndex];
            }
            return result;
        }

        // Matrix Multiplcation
        Vector operator*(const Vector& vec) const {
            checkDimensionCompatibility(*this, vec);
            Vector result(cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    result[j] += vals[i][j] * vec[j];
                }
            }
            return result;
        }
        friend Vector operator*(const Vector& vec, const Matrix& mat) {
            checkDimensionCompatibility(vec, mat);
            Vector result(mat.cols);
            for(int i = 0; i < mat.rows; i++) {
                for(int j = 0; j < mat.cols; j++) {
                    result[j] += mat.vals[i][j] * vec[j];
                }
            }
            return result;
        }
        Point operator*(const Point& point) const {
            checkDimensionCompatibility(*this, point);
            Point result(cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    result[j] += vals[i][j] * point[j];
                }
            }
            return result;
        }
        friend Point operator*(const Point& point, const Matrix& mat) {
            checkDimensionCompatibility(point, mat);
            Point result(mat.cols);
            for(int i = 0; i < mat.rows; i++) {
                for(int j = 0; j < mat.cols; j++) {
                    result[j] += mat.vals[i][j] * point[j];
                }
            }
            return result;
        }
        Matrix operator*(const Matrix& mat) const {
            checkDimensionCompatibility(*this, mat);
            Matrix result(rows, mat.cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < mat.cols; j++) {
                    for(int k = 0; k < cols; k++) {
                        result[i][j] += vals[i][k] * mat.vals[k][j];
                    }
                }
            }
            return result;
        }
        template<typename T> Matrix operator*(const T scale) const {
            Matrix result(rows, cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    result[i][j] = vals[i][j] * scale;
                }
            }
            return result;
        }
        template<typename T> friend Matrix operator*(const T scale, const Matrix& mat) {
            Matrix result(mat.rows, mat.cols);
            for(int i = 0; i < mat.rows; i++) {
                for(int j = 0; j < mat.cols; j++) {
                    result[i][j] = mat[i][j] * scale;
                }
            }
            return result;
        }
        template<typename T> Matrix& operator*=(const T scale) {
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    vals[i][j] *= scale;
                }
            }
            return *this;
        }

        // Matrix Division
        template<typename T> Matrix  operator/(const T scale) const {
            Matrix result(rows, cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    result[i][j] = vals[i][j] / scale;
                }
            }
            return result;
        }
        template<typename T> Matrix& operator/=(const T scale) {
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    vals[i][j] /= scale;
                }
            }
            return *this;
        }

        // Matrix Addition
        Matrix  operator+(const Matrix& mat) const {
            checkDimensionEquality(*this, mat);
            Matrix result(rows, cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    result[i][j] = vals[i][j] + mat.vals[i][j];
                }
            }
            return result;
        }
        Matrix& operator+=(const Matrix& mat) {
            checkDimensionEquality(*this, mat);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    vals[i][j] += mat.vals[i][j];
                }
            }
            return *this;
        }

        // Matrix Subtraction
        Matrix  operator-(const Matrix& mat) const {
            checkDimensionEquality(*this, mat);
            Matrix result(rows, cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    result[i][j] = vals[i][j] - mat.vals[i][j];
                }
            }
            return result;
        }
        Matrix& operator-=(const Matrix& mat) {
            checkDimensionEquality(*this, mat);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    vals[i][j] -= mat.vals[i][j];
                }
            }
            return *this;
        }

        // Matrix Equality
        bool operator==(const Matrix& mat) const {
            return equals(mat);
        }
        bool operator!=(const Matrix& mat) const {
            return !equals(mat);
        }
        bool equals(const Matrix& mat) const {
            return equals(mat, 0.0000000001);
        }
        bool equals(const Matrix& mat, const double precision) const {
            if(rows != mat.rows || cols != mat.cols) {
                return false;
            }
            for(int i = 0; i < rows; i++) {
                if(!vals[i].equals(mat.vals[i], precision)) {
                    return false;
                }
            }
            return true;
        }

        // Matrix Status
        int    numRows() const {
            return rows;
        }
        int    numCols() const {
            return cols;
        }
        double min() const {
            checkBounds(0);
            double min = vals[0][0];
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    double val = vals[i][j];
                    if(val < min) {
                        min = val;
                    }
                }
            }
            return min;
        }
        double max() const {
            checkBounds(0);
            double max = vals[0][0];
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    double val = vals[i][j];
                    if(val > max) {
                        max = val;
                    }
                }
            }
            return max;
        }
        bool   isSquare() const {
            return rows == cols;
        }
        bool   isInvertible() const {
            return determinant() != 0;
        }
        double determinant() const {
            if(!isSquare()) {
                throw std::runtime_error("Cannot calculate the determinant of a non-square matrix.");
            }
            if(rows == 2 && cols == 2) {
                return (vals[0][0] * vals[1][1] - vals[1][0] * vals[0][1]);
            }
            double sum = 0;
            double sign = 1.0;
            for(int i = 0; i < cols; i++) {
                Matrix subMatrix = getMinorMatrix(0, i);
                sum += vals[0][i] * subMatrix.determinant() * sign;
                sign *= -1;
            }
            return sum;
        }
        Matrix getTransposedMatrix() const {
            Matrix newMatrix(cols, rows);
            for(int i = 0; i < rows; i++) {
                newMatrix.setColumn(i, vals[i]);
            }
            return newMatrix;
        }
        Matrix getAdjugateMatrix() const {
            if(!isSquare()) {
                throw std::runtime_error("cannot find a cofactor matrix for a non-square matrix");
            }
            Matrix result(rows, cols);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    Matrix minor = getMinorMatrix(i, j);
                    result[j][i] = minor.determinant() * ((i+j) % 2 == 0 ? 1 : -1);
                }
            }
            return result;
        }
        Matrix getInverse() const {
            if(!isSquare()) {
                throw std::runtime_error("cannot find the inverse of a non-square matrix");
            }
            double det = determinant();
            if(det == 0) {
                throw std::runtime_error("Cannot find the inverse of a non-singular matrix");
            }
            Matrix adjugateMatrix = getAdjugateMatrix();
            return adjugateMatrix / det;
        }

        // Matrix operations
        Matrix& fillMatrix(const double fillValue) {
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    vals[i][j] = fillValue;
                }
            }
            return *this;
        }
        Matrix& fillRow(const int rowIndex, const double fillValue) {
            checkBounds(rowIndex);
            for(int i = 0; i < cols; i++) {
                vals[rowIndex][i] = fillValue;
            }
            return *this;
        }
        Matrix& fillColumn(const int columnIndex, const double fillValue) {
            // don't need to do bounds checking since it gets handled in the vector
            for(int i = 0; i < rows; i++) {
                vals[i][columnIndex] = fillValue;
            }
            return *this;
        }
        Matrix& setRow(const int rowIndex, const Vector& newRow) {
            if(newRow.getSize() != cols) {
                throwDimensionError(*this, newRow);
            }
            checkBounds(rowIndex);
            for(int i = 0; i < cols; i++) {
                vals[rowIndex][i] = newRow[i];
            }
            return *this;
        }
        Matrix& setColumn(const int columnIndex, const Vector& newCol) {
            if(newCol.getSize() != rows) {
                throwDimensionError(*this, newCol);
            }
            for(int i = 0; i < rows; i++) {
                vals[i][columnIndex] = newCol[i];
            }
            return *this;
        }
        Matrix& changeSize(const int newRows, const int newCols) {
            if(newCols != cols) {
                for(int i = 0; i < rows && i < newRows; i++) {
                    vals[i].changeSize(newCols);
                }
                cols = newCols;
            }
            if(newRows != rows) {
                Vector* newVals = new Vector[newRows];
                if(newRows < rows) {
                    for(int i = 0; i < newRows; i++) {
                        newVals[i] = vals[i];
                    }
                } else {
                    for(int i = 0; i < rows; i++) {
                        newVals[i] = vals[i];
                    }
                    for(int i = rows; i < newRows; i++) {
                        newVals[i] = Vector(cols);
                    }
                }
                rows = newRows;
                vals = newVals;
            }
            return *this;
        }
        Matrix& addRow() {
            return changeSize(rows+1, cols);
        }
        Matrix& addRow(const Vector& newRow) {
            changeSize(rows+1, cols);
            vals[rows-1] = newRow;
            return *this;
        }
        Matrix& addColumn() {
            return changeSize(rows, cols+1);
        }
        Matrix& addColumn(const Vector& newColumn) {
            changeSize(rows, cols+1);
            for(int i = 0; i < rows; i++) {
                vals[i][cols-1] =  newColumn[i];
            }
            return *this;
        }
        Matrix& removeRow(const int rowIndex) {
            Matrix temp(rows-1, cols);
            for(int i = 0; i < rowIndex; i++) {
                for(int j = 0; j < cols; j++) {
                    temp[i][j] = vals[i][j];
                }
            }
            for(int i = rowIndex+1; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                    temp[i-1][j] = vals[i][j];
                }
            }
            *this = temp;
            return *this;
        }
        Matrix& removeColumn(const int columnIndex) {
            Matrix temp(rows, cols-1);
            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < columnIndex; j++) {
                    temp[i][j] = vals[i][j];
                }
            }
            for(int i = 0; i < rows; i++) {
                for(int j = columnIndex+1; j < cols; j++) {
                    temp[i][j-1] = vals[i][j];
                }
            }
            *this = temp;
            return *this;
        }
        Vector  getPartialRow(const int rowIndex, const int startIndex, const int endIndex) const {
            checkBounds(rowIndex);
            return vals[rowIndex].getValues(startIndex, endIndex);
        }
        Vector  getPartialColumn(const int colIndex, const int startIndex, const int endIndex) const {
            checkBounds(startIndex, endIndex);
            return getCol(colIndex).getValues(startIndex, endIndex);
        }
        Matrix& transpose() {
            *this = getTransposedMatrix();
            return *this;
        }
        Matrix& invert() {
            *this = getInverse();
            return *this;
        }

        // Static Matrix functions
        static Matrix identity(const int numRows, const int numCols) {
            if(numRows != numCols) {
                throw std::runtime_error("cannot create a non-square identity matrix.");
            }
            Matrix result(numRows, numCols);
            for(int i = 0; i < numRows; i++) {
                result[i][i] = 1.0;
            }
            return result;
        }

    };
}

#endif //GRAPHICSUTILITIES_MATRIX_H
