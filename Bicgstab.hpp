/**
 * @file Bicgstab.hpp
 * @author Chiamaka Okeke 
 * @brief Biconjugate gradient stabilized method
 * @version 0.1
 * @date 2020-12-13
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;
/**
 * @brief A class for matrices
 * 
 */
class matrix
{

public:
    /**
     * @brief Construct an uninitialized matrix object for both square and non-square matrices
     * 
     */
    matrix(const size_t &);
    matrix(const size_t &, const size_t &);

    /**
     * @brief Initializing matrix with flattened vector
     * 
     */
    matrix(const size_t &, const vector<double> &);
    matrix(const size_t &, const size_t &, const vector<double> &);

    /**
     * @brief Member function to obtatin the number of rows and columns
     * 
     * @return size 
     */
    size_t get_rows() const;
    size_t get_cols() const;

    /**
     * @brief Member function to determine vector norm
     * 
     * @return returns norm of a vector 
     */
    double l2_norm(vector<double> const &);

    /**
     * @brief Acessing the matrices with operator()
     * 
     * @return Returns individual elements of the matrix 
     */
    double &operator()(const size_t &, const size_t &);
    double operator()(const size_t &, const size_t &) const;
    double &operator()(size_t &i);

    /**
     * @brief Overloading operator* as a friend function to enable multiplication between matrices and vectors
     * 
     * @param A Matrix
     * @param x Vector
     * @return Returns a vector
     */

    friend vector<double> operator*(const matrix &A, const vector<double> &x);
    friend vector<double> operator*(const vector<double> &x, const matrix &A);
    /**
     * @brief Declaration of the biconjugate gradient stabilized solver as a friend function
     * 
     * @param A Matrix
     * @param b Vector
     * @param verbose True
     * @param TOL Tolerance
     * @param maxit Maximum iteration
     * @return Returns a vector as approximate value
     */
    //friend vector<double> bicgstab(const matrix &A, const vector<double> &b, bool verbose, double TOL, int maxit);
    friend vector<double> bicgstab(const matrix &A, const vector<double> &b, const vector<double> &x0, bool verbose, double TOL, int maxit);

    // Exception to be thrown no value was assigned to rows and columns
    class wrong_parameter
    {
    };
    // Exception to be thrown if the vector of elements provided to the constructor is of the wrong size
    class initializer_wrong_size
    {
    };
    // Exception to be thrown if two matrices of different sizes are added
    class incompatible_sizes_add
    {
    };
    // Exception to be thrown if two matrices of different sizes are subtracted
    class incompatible_sizes_subtract
    {
    };
    // Exception to be thrown if two vectors of incompartible sizes are utilized
    class incompatible_vectors
    {
    };
    // Exception to be thrown if two matrices of incompatible sizes are multiplied
    class incompatible_sizes_multiply
    {
    };
    // Exception to be thrown if matrix and vector of incompatible sizes are multiplied
    class incompatible_sizes_mat_vec_multiply
    {
    };
    // Exception to be thrown when matrix will be divided by zero
    class division_by_zero_not_allowed
    {
    };
    // Exception to be thrown when matrix will be divided by zero
    class condition_for_bicgstab_not_satisfied
    {
    };

    /**
     * @brief The private member consists of the number of rows, columns and the input data that was flattened to 1-dimension
     * 
     */
private:
    size_t rows{0};
    size_t cols{0};
    vector<double> elements;
};

/**
 * @brief Overloading of operators for the matrix
 * 
 */
ostream &operator<<(ostream &, const matrix &);
matrix operator+(const matrix &, const matrix &);
matrix operator+=(matrix &, const matrix &);
matrix operator-(const matrix &);
matrix operator-(const matrix &, const matrix &);
matrix operator-=(matrix &, const matrix &);
matrix operator/(const matrix &, const double &);
matrix operator*(const matrix &, const matrix &);
matrix operator*(const double &, const matrix &);
matrix operator*(const matrix &, const double &);

/**
 * @brief Declaration of the matrix-vector multiplication
 * 
 * @param A Matrix
 * @param x Vector
 * @return Returns a vector 
 */
vector<double> operator*(const matrix &A, const vector<double> &x);
vector<double> operator*(const vector<double> &x, const matrix &A);

/**
 * @brief Overloading of operators for the vector
 */
//Operator overload for vector
ostream &operator<<(ostream &, const vector<double> &);
vector<double> operator+(const vector<double> &, const vector<double> &);
vector<double> operator+=(vector<double> &, const vector<double> &);
vector<double> operator-(const vector<double> &);
vector<double> operator-(const vector<double> &, const vector<double> &);
vector<double> operator-=(vector<double> &, const vector<double> &);
double operator*(const vector<double> &, const vector<double> &);
vector<double> operator*(const double &, const vector<double> &);
double inner_product(const vector<double> &U, const vector<double> &V);
vector<double> operator*(const vector<double> &, const double &);
vector<double> operator/(const vector<double> &, const double &);

/**
 * @brief Actual declaration of the biconjugate gradient stabilized solver 
 * 
 * @param A Matrix
 * @param b Vector
 * @param verbose True
 * @param TOL Tolerance
 * @param maxit Maximun iteration
 * @return Returns a vector 
 */
vector<double> bicgstab(const matrix &A, const vector<double> &b, bool verbose = 0, double TOL = 1E-12, int maxit = 256);
vector<double> bicgstab(const matrix &A, const vector<double> &b, const vector<double> &x0, bool verbose = 0, double TOL = 1E-12, int maxit = 256);
