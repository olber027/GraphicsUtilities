//
// Created by molberding on 6/25/2018.
//

#ifndef GRAPHICSUTILITIES_VECTOR_H
#define GRAPHICSUTILITIES_VECTOR_H

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <vector>

namespace linear_algebra {

    // The vectors, by default, will be created as a column vector but they can be transposed to a row vector.
    // an implicit transposition is done for vector-vector operations if necessary.
    class Vector {
    private:
        int size;
        double* vals;
        int memorySize;
        int rows, cols;

        void destroy() {
            if(vals != nullptr) {
                delete [] vals;
                size = 0;
                vals = nullptr;
                memorySize = 0;
                rows = 0;
                cols = 0;
            }
        }

        void initialize(const int newSize) {
            destroy();
            size = newSize;
            vals = new double[size];
            memorySize = size;
            rows = size;
            cols = 1;
        }

        void extend(const int numValsToAdd, const bool addToFront) {
            if(numValsToAdd <= 0) {
                return;
            }

            int newSize = size + numValsToAdd;
            double* newVals = nullptr;

            if(memorySize <= newSize) {
                if (memorySize <= 0) {
                    memorySize = 1;
                }
                while (memorySize <= newSize) {
                    memorySize *= 2;
                }
                newVals = new double[memorySize];
            }

            if(newVals != nullptr) {
                int offset = 0;
                if(addToFront) {
                    offset += numValsToAdd;
                }
                int i = offset;
                for (i; i < size + offset; i++) {
                    newVals[i] = vals[i - offset];
                }
                delete[] vals;
                vals = newVals;
            } else if(addToFront) {
                for(int i = size; i >= 0; i--) {
                    vals[i + numValsToAdd] = vals[i];
                }
            }

            size += numValsToAdd;
            if(cols > rows) {
                cols = size;
            } else {
                rows = size;
            }
        }

        void extend(const int numValsToAdd) {
            extend(numValsToAdd, false);
        }

        void checkBounds(const int index) const {
            if(index < 0 || index >= size) {
                std::stringstream str;
                str << "index out of range: " << index;
                throw std::out_of_range(str.str());
            }
        }

        template<typename... Args>
        void checkBounds(const int index, Args... args) const {
            if(index < 0 || index >= size) {
                std::stringstream str;
                str << "index out of range: " << index;
                throw std::out_of_range(str.str());
            }
            checkBounds(args...);
        }

        void throwSizeMismatchError(const Vector& vec1, const Vector& vec2) const {
            std::stringstream str;
            str << "Vectors must be of the same size to perform calculations." << std::endl;
            str << "Given vectors were: " << std::endl;
            str << vec1 << std::endl;
            str << vec2 << std::endl;
            throw std::runtime_error(str.str());
        }

        double round(const double val, const double precision) const {
            if(fabs(val - std::round(val)) < precision) {
                return std::round(val);
            }
            return val;
        }

    public:

        Vector() : size(0), vals(nullptr), memorySize(0), rows(0), cols(0) {}
        Vector(const int initialSize) : Vector() {
            initialize(initialSize);
        }
        Vector(const int initialSize, const bool rowVector) : Vector() {
            initialize(initialSize);
            if(rowVector) {
                transpose();
            }
        }
        Vector(const int initialSize, const double fillValue) : Vector(initialSize) {
            fillVector(fillValue);
        }
        Vector(const int initialSize, const double fillValue, const bool rowVector) : Vector(initialSize) {
            fillVector(fillValue);
            if(rowVector) {
                transpose();
            }
        }
        Vector(double a, double b) : Vector(2) {
            vals[0] = a;
            vals[1] = b;
        }
        Vector(double a, double b, double c) : Vector(3) {
            vals[0] = a;
            vals[1] = b;
            vals[2] = c;
        }
        Vector(const double* source, const int sourceSize) : Vector(sourceSize) {
            for(int i = 0; i < sourceSize; i++) {
                vals[i] = source[i];
            }
        }
        Vector(const std::vector<double> source) : Vector(source.size()) {
            for(int i = 0; i < source.size(); i++) {
                vals[i] = source[i];
            }
        }
        Vector(const Vector& other) : Vector() {
            initialize(other.size);
            for(int i = 0; i < size; i++) {
                vals[i] = other.vals[i];
            }
            if(other.isRowVector()) {
                transpose();
            }
        }
        Vector(Vector&& other) noexcept {
            vals = other.vals;
            size = other.size;
            memorySize = other.memorySize;
        }
        ~Vector() {
            destroy();
        }

        double  getValue(const int index) const {
            checkBounds(index);
            return vals[index];
        }
        Vector& setValue(const int index, const double newVal) {
            checkBounds(index);
            vals[index] = newVal;
            return *this;
        }
        Vector  getValues(const int startIndex, const int endIndex) const {
            checkBounds(startIndex, endIndex);
            Vector result;
            for(int i = startIndex; i <= endIndex; i++) {
                result.pushBack(vals[i]);
            }
            return result;
        }

        Vector& pushBack(const double newVal) {
            extend(1);
            vals[size-1] = newVal;
            return *this;
        }
        template<typename... Args>
        Vector& pushBack(const double newVal, Args... args) {
            pushBack(newVal);
            pushBack(args...);
            return *this;
        }
        Vector& pushBack(const Vector& newVals) {
            int oldSize = size;
            extend(newVals.size);
            for(int i = 0; i < newVals.size; i++) {
                vals[oldSize + i] = newVals[i];
            }
            return *this;
        }
        Vector& pushBack(const double* newVals, const int newValsSize) {
            int oldSize = size;
            extend(newValsSize);
            for(int i = 0; i < newValsSize; i++) {
                vals[oldSize + i] = newVals[i];
            }
            return *this;
        }
        Vector& pushBack(const std::vector<double>& newVals) {
            int oldSize = size;
            extend(newVals.size());
            for(int i = 0; i < newVals.size(); i++) {
                vals[oldSize + i] = newVals[i];
            }
            return *this;
        }

        Vector& pushFront(const double newVal) {
            extend(1, true);
            vals[0] = newVal;
            return *this;
        }
        Vector& pushFront(const Vector& newVals) {
            extend(newVals.size, true);
            for(int i = 0; i < newVals.size; i++) {
                vals[i] = newVals[i];
            }
            return *this;
        }
        Vector& pushFront(const double* newVals, const int newValsSize) {
            extend(newValsSize, true);
            for(int i = 0; i < newValsSize; i++) {
                vals[i] = newVals[i];
            }
            return *this;
        }
        Vector& pushFront(const std::vector<double>& newVals) {
            extend(newVals.size());
            for(int i = 0; i < newVals.size(); i++) {
                vals[i] = newVals[i];
            }
            return *this;
        }

        double& operator[](const int index) const {
            checkBounds(index);
            return vals[index];
        }

        // vector status
        bool   is2d() const {
            return size == 2;
        }
        bool   is3d() const {
            return size == 3;
        }
        bool   isUnit() const {
            return magnitude() == 1;
        }
        bool   isZero() const {
            for(int i = 0; i < size; i++) {
                if(vals[i] != 0) {
                    return false;
                }
            }
            return true;
        }
        bool   isEmpty() const {
            return size == 0;
        }
        bool   isRowVector() const {
            if(rows == 1) {
                return true;
            }
            return false;
        }
        bool   isColumnVector() const {
            if(cols == 1) {
                return true;
            }
            return false;
        }
        int    numRows() const {
            return rows;
        }
        int    numCols() const {
            return cols;
        }
        int    getSize() const {
            return size;
        }
        int    getMemorySize() const {
            return memorySize;
        }
        double magnitude() const {
            double sum = 0;
            for(int i = 0; i < size; i++) {
                sum += (vals[i] * vals[i]);
            }
            return sqrt(sum);
        }
        double min() const {
            if(size == 0) {
                std::runtime_error("Cannot find the min of an empty vector");
            }
            double result = vals[0];
            for(int i = 1; i < size; i++) {
                if(vals[i] < result) {
                    result = vals[i];
                }
            }
            return result;
        }
        double max() const {
            if(size == 0) {
                std::runtime_error("Cannot find the max of an empty vector");
            }
            double result = vals[0];
            for(int i = 1; i < size; i++) {
                if(vals[i] > result) {
                    result = vals[i];
                }
            }
            return result;
        }

        friend std::ostream& operator<<(std::ostream& out, const Vector& vec) {
            if(vec.size > 0) {
                out << "[ " << vec.vals[0];
                for(int i = 1; i < vec.size; i++) {
                    out << ", " << vec.vals[i];
                }
                out << "]";
            } else {
                out << "[ ]";
            }
            return out;
        }

        Vector& operator=(const Vector& rhs) {
            if(this != &rhs) {
                initialize(rhs.size);
                for(int i = 0; i < size; i++) {
                    vals[i] = rhs.vals[i];
                }
            }
            return *this;
        }
        Vector& operator=(const std::vector<double>& rhs) {
            initialize(rhs.size());
            for(int i = 0; i < size; i++) {
                vals[i] = rhs[i];
            }
            return *this;
        }

        bool operator==(const Vector& rhs) const {
            return equals(rhs);
        }
        bool operator!=(const Vector& rhs) const {
            return !equals(rhs);
        }
        bool equals(const Vector& rhs) const {
            return equals(rhs, 0.0000000001);
        }
        bool equals(const Vector& rhs, const double precision) const {
            if(rhs.size != size) {
                return false;
            }
            for(int i = 0; i < size; i++) {
                if(round(vals[i], precision) != round(rhs.vals[i], precision)) {
                    return false;
                }
            }
            return true;
        }

        // scalar multiplication
        template<typename T> Vector  operator*(const T scale) const {
            Vector result(size);
            for(int i = 0; i < size; i++) {
                result.vals[i] = vals[i] * scale;
            }
            return result;
        }
        template<typename T> friend Vector operator*(const T scale, const Vector& vec) {
            Vector result(vec.size);
            for(int i = 0; i < vec.size; i++) {
                result.vals[i] = vec.vals[i] * scale;
            }
            return result;
        }
        template<typename T> Vector& operator*=(const T scale) {
            for(int i = 0; i < size; i++) {
                vals[i] *= scale;
            }
            return *this;
        }

        // scalar division
        template<typename T> Vector  operator/(const T scale) const {
            Vector result(size);
            for(int i = 0; i < size; i++) {
                result.vals[i] = vals[i] / scale;
            }
            return result;
        }
        template<typename T> Vector& operator/=(const T scale) {
            for(int i = 0; i < size; i++) {
                vals[i] /= scale;
            }
            return *this;
        }

        // dot product
        double operator*(const Vector& rhs) const {
            if(size != rhs.size) {
                throwSizeMismatchError(*this, rhs);
            }
            double sum = 0;
            for(int i = 0; i < size; i++) {
                sum += (vals[i] * rhs.vals[i]);
            }
            return sum;
        }
        double dot(const Vector& other) const {
            if(size != other.size) {
                throwSizeMismatchError(*this, other);
            }
            double sum = 0;
            for(int i = 0; i < size; i++) {
                sum += (vals[i] * other.vals[i]);
            }
            return sum;
        }

        // cross product
        Vector cross(const Vector& rhs) const {
            if(!is3d() || !rhs.is3d()) {
                throw std::runtime_error("cross product of anything not a 3d vector is meaningless");
            }

            Vector result(3);
            result.vals[0] = (vals[1] * rhs.vals[2]) - (vals[2] * rhs.vals[1]);
            result.vals[1] = (vals[2] * rhs.vals[0]) - (vals[0] * rhs.vals[2]);
            result.vals[2] = (vals[0] * rhs.vals[1]) - (vals[1] * rhs.vals[0]);

            return result;
        }

        // vector addition
        Vector  operator+(const Vector& rhs) const {
            if(size != rhs.size) {
                throwSizeMismatchError(*this, rhs);
            }
            Vector result(size);
            for(int i = 0; i < size; i++) {
                result.vals[i] = vals[i] + rhs.vals[i];
            }
            return result;
        }
        Vector& operator+=(const Vector& rhs) {
            if(size != rhs.size) {
                throwSizeMismatchError(*this, rhs);
            }
            for(int i = 0; i < size; i++) {
                vals[i] += rhs.vals[i];
            }
            return *this;
        }
        Vector& add(const Vector& other) {
            if(size != other.size) {
                throwSizeMismatchError(*this, other);
            }
            for(int i = 0; i < size; i++) {
                vals[i] += other.vals[i];
            }
            return *this;
        }

        // vector subtraction
        Vector  operator-(const Vector& rhs) const {
            if(size != rhs.size) {
                throwSizeMismatchError(*this, rhs);
            }
            Vector result(size);
            for(int i = 0; i < size; i++){
                result.vals[i] = vals[i] - rhs.vals[i];
            }
            return result;
        }
        Vector& operator-=(const Vector& rhs) {
            if(size != rhs.size) {
                throwSizeMismatchError(*this, rhs);
            }
            for(int i = 0; i < size; i++) {
                vals[i] -= rhs.vals[i];
            }
            return *this;
        }
        Vector& subtract(const Vector& other) {
            if(size != other.size) {
                throwSizeMismatchError(*this, other);
            }
            for(int i = 0; i < size; i++) {
                vals[i] -= other.vals[i];
            }
            return *this;
        }

        // vector operations
        Vector& transpose() {
            int temp = rows;
            rows = cols;
            cols = temp;
            return *this;
        }
        Vector& normalize() {
            double mag = magnitude();
            for(int i = 0; i < size; i++) {
                vals[i] /= mag;
            }
            return *this;
        }
        Vector  getNormalized() const {
            Vector result(size);
            double mag = magnitude();
            for(int i = 0; i < size; i++) {
                result.vals[i] = vals[i] / mag;
            }
            return result;
        }
        Vector& toUnitVector() {
            normalize();
            return *this;
        }
        Vector  getUnitVector() const {
            return getNormalized();
        }
        Vector& abs() {
            for(int i = 0; i < size; i++) {
                if(vals[i] < 0) {
                    vals[i] *= -1;
                }
            }
            return *this;
        }
        Vector& changeSize(const int newSize) {
            if(newSize <= 0) {
                destroy();
                return *this;
            }
            if(newSize > size) {
                extend(newSize);
            } else if(newSize < size) {
                double* temp = new double[newSize];
                for(int i = 0; i < newSize; i++) {
                    temp[i] = vals[i];
                }
                destroy();
                size = newSize;
                vals = temp;
                memorySize = size;
            }
            return *this;
        }
        Vector& fillVector(const double fillValue) {
            for(int i = 0; i < size; i++) {
                vals[i] = fillValue;
            }
            return *this;
        }
        Vector& zero() {
            return fillVector(0);
        }
        Vector& toVector4() {
            // This is meant to be used for graphics purposes. most transforms need the 3d vector to have a fourth
            // dimension with a value of 1. so... this does that.
            changeSize(4);
            vals[3] = 1;
            return *this;
        }
        Vector& toVector3() {
            // This is here to undo the toVector4 command
            changeSize(3);
            return *this;
        }
        void    clear() {
            destroy();
        }
    };
}

#endif //GRAPHICSUTILITIES_VECTOR_H
