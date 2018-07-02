//
// Created by molberding on 6/25/2018.
//

#ifndef GRAPHICSUTILITIES_POINT_H
#define GRAPHICSUTILITIES_POINT_H

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <vector>
#include "Vector.h"

namespace linear_algebra {

    class Point {
    private:
        int size;
        double* vals;
        int memorySize;

        void destroy() {
            if(vals != nullptr) {
                delete [] vals;
                size = 0;
                vals = nullptr;
            }
        }

        void initialize(const int newSize) {
            destroy();
            size = newSize;
            vals = new double[size];
            memorySize = size;
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

        void throwSizeMismatchError(const Point& point, const Vector& vec) const {
            std::stringstream str;
            str << "Vector and point must be of the same size to perform calculations." << std::endl;
            str << "Given vector and point were: " << std::endl;
            str << point << std::endl;
            str << vec << std::endl;
            throw std::runtime_error(str.str());
        }

    public:

        Point() : size(0), vals(nullptr), memorySize(0) {}
        Point(const int initialSize) : Point() {
            initialize(initialSize);
        }
        Point(const int initialSize, const double fillValue) : Point(initialSize) {
            fillPoint(fillValue);
        }
        Point(double a, double b) : Point(2) {
            vals[0] = a;
            vals[1] = b;
        }
        Point(double a, double b, double c) : Point(3) {
            vals[0] = a;
            vals[1] = b;
            vals[2] = c;
        }
        Point(const double* source, const int sourceSize) : Point(sourceSize) {
            for(int i = 0; i < size; i++) {
                vals[i] = source[i];
            }
        }
        Point(const std::vector<double>& source) : Point(source.size()) {
            for(int i = 0; i < size; i++) {
                vals[i] = source[i];
            }
        }
        Point(const Point& other) : Point() {
            if(&other != this) {
                initialize(other.size);
                for(int i = 0; i < size; i++) {
                    vals[i] = other.vals[i];
                }
            }
        }
        Point(const Point&& other) noexcept {
            size = other.size;
            vals = other.vals;
            memorySize = other.memorySize;
        }
        ~Point() {
            destroy();
        }

        double getValue(const int index) const {
            checkBounds(index);
            return vals[index];
        }
        Point& setValue(const int index, const double newVal) {
            checkBounds(index);
            vals[index] = newVal;
            return *this;
        }
        Point getValues(const int startIndex, const int endIndex) const {
            checkBounds(startIndex, endIndex);
            Point result;
            for(int i = startIndex; i <= endIndex; i++) {
                result.pushBack(vals[i]);
            }
            return result;
        }

        Point& pushBack(const double newVal) {
            extend(1);
            vals[size-1] = newVal;
            return *this;
        }
        Point& pushBack(const Point& newVals) {
            int oldSize = size;
            extend(newVals.size);
            for(int i = 0; i < newVals.size; i++) {
                vals[oldSize + i] = newVals[i];
            }
            return *this;
        }
        Point& pushBack(const double* newVals, const int newValsSize) {
            int oldSize = size;
            extend(newValsSize);
            for(int i = 0; i < newValsSize; i++) {
                vals[oldSize + i] = newVals[i];
            }
            return *this;
        }
        Point& pushBack(const std::vector<double>& newVals) {
            int oldSize = size;
            extend(newVals.size());
            for(int i = 0; i < newVals.size(); i++) {
                vals[oldSize + i] = newVals[i];
            }
            return *this;
        }

        Point& pushFront(const double newVal) {
            extend(1, true);
            vals[0] = newVal;
            return *this;
        }
        Point& pushFront(const Point& newVals) {
            extend(newVals.size, true);
            for(int i = 0; i < newVals.size; i++) {
                vals[i] = newVals[i];
            }
            return *this;
        }
        Point& pushFront(const double* newVals, const int newValsSize) {
            extend(newValsSize, true);
            for(int i = 0; i < newValsSize; i++) {
                vals[i] = newVals[i];
            }
            return *this;
        }
        Point& pushFront(const std::vector<double>& newVals) {
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

        // Point status
        bool is2d() const {
            return size == 2;
        }
        bool is3d() const {
            return size == 3;
        }
        bool isZero() const {
            for(int i = 0; i < size; i++) {
                if(vals[i] != 0) {
                    return false;
                }
            }
            return true;
        }
        bool isEmpty() const {
            return size == 0;
        }
        int getSize() const {
            return size;
        }
        int getMemorySize() const {
            return memorySize;
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

        friend std::ostream& operator<<(std::ostream& out, const Point& point) {
            if(point.size > 0) {
                out << "[ " << point.vals[0];
                for(int i = 1; i < point.size; i++) {
                    out << ", " << point.vals[i];
                }
                out << "]";
            } else {
                out << "[ ]";
            }
            return out;
        }

        Point& operator=(const Point& rhs) {
            if(this != &rhs) {
                initialize(rhs.size);
                for(int i = 0; i < size; i++) {
                    vals[i] = rhs.vals[i];
                }
            }
            return *this;
        }
        Point& operator=(const std::vector<double>& rhs) {
            initialize(rhs.size());
            for(int i = 0; i < size; i++) {
                vals[i] = rhs[i];
            }
            return *this;
        }

        bool operator==(const Point& rhs) {
            if(size != rhs.size) {
                return false;
            }
            for(int i = 0; i < size; i++) {
                if(vals[i] != rhs.vals[i]) {
                    return false;
                }
            }
            return true;
        }
        bool operator!=(const Point& rhs) {
            if(size != rhs.size) {
                return true;
            }
            for(int i = 0; i < size; i++) {
                if(vals[i] != rhs.vals[i]) {
                    return true;
                }
            }
            return false;
        }

        // scalar multiplication
        template<typename T> Point operator*(const T scale) const {
            Point result(size);
            for(int i = 0; i < size; i++) {
                result.vals[i] = vals[i] * scale;
            }
            return result;
        }
        template<typename T> friend Point operator*(const T scale, const Point& point) {
            Point result(point.size);
            for(int i = 0; i < point.size; i++) {
                result.vals[i] = point.vals[i] * scale;
            }
            return result;
        }
        template<typename T> Point& operator*=(const T scale) {
            for(int i = 0; i < size; i++) {
                vals[i] *= scale;
            }
            return *this;
        }

        // scalar division
        template<typename T> Point operator/(const T scale) const {
            Point result(size);
            for(int i = 0; i < size; i++) {
                result.vals[i] = vals[i] / scale;
            }
            return result;
        }
        template<typename T> Point& operator/=(const T scale) {
            for(int i = 0; i < size; i++) {
                vals[i] /= scale;
            }
            return *this;
        }

        // vector addition
        Point operator+(const Vector& vec) const {
            if(size != vec.getSize()) {
                throwSizeMismatchError(*this, vec);
            }
            Point result(size);
            for(int i = 0; i < size; i++) {
                result.vals[i] = vals[i] + vec[i];
            }
            return result;
        }
        Point& operator+=(const Vector& vec) {
            if(size != vec.getSize()) {
                throwSizeMismatchError(*this, vec);
            }
            for(int i = 0; i < size; i++) {
                vals[i] += vec[i];
            }
            return *this;
        }
        Point& add(const Vector& vec) {
            if(size != vec.getSize()) {
                throwSizeMismatchError(*this, vec);
            }
            for(int i = 0; i < size; i++) {
                vals[i] += vec[i];
            }
            return *this;
        }

        // vector subtraction
        Point operator-(const Vector& vec) const {
            if(size != vec.getSize()) {
                throwSizeMismatchError(*this, vec);
            }
            Point result(size);
            for(int i = 0; i < size; i++){
                result.vals[i] = vals[i] - vec[i];
            }
            return result;
        }
        Point& operator-=(const Vector& vec) {
            if(size != vec.getSize()) {
                throwSizeMismatchError(*this, vec);
            }
            for(int i = 0; i < size; i++) {
                vals[i] -= vec[i];
            }
            return *this;
        }
        Point& subtract(const Vector& vec) {
            if(size != vec.getSize()) {
                throwSizeMismatchError(*this, vec);
            }
            for(int i = 0; i < size; i++) {
                vals[i] -= vec[i];
            }
            return *this;
        }

        // point operations
        Point& abs() {
            for(int i = 0; i < size; i++) {
                if(vals[i] < 0) {
                    vals[i] *= -1;
                }
            }
            return *this;
        }
        Point& changeSize(const int newSize) {
            if(newSize <= 0) {
                destroy();
                return *this;
            }
            if(newSize > size) {
                extend(newSize);
            } else {
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
        Point& fillPoint(const double fillValue) {
            for(int i = 0; i < size; i++) {
                vals[i] = fillValue;
            }
            return *this;
        }
        Point& zero() {
            return fillPoint(0);
        }
        Point& toPoint4() {
            // This is meant to be used for graphics purposes. most transforms need the 3d point to have a fourth
            // dimension with a value of 1. so... this does that.
            changeSize(4);
            vals[3] = 1;
            return *this;
        }
        Point& toPoint3() {
            // This is here to undo the toPoint4 command
            changeSize(3);
            return *this;
        }
        void clear() {
            destroy();
        }

    };
}

#endif //GRAPHICSUTILITIES_POINT_H
