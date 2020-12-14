/**
 * @file Bicgstab.cpp
 * @author Chiamaka Okeke
 * @brief Biconjugate gradient stabilized method
 * @version 0.1
 * @date 2020-12-13
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "Bicgstab.hpp"

using namespace std;

/**
 * @brief Implementation of the constructor for a square matrix
 * 
 * @param size equals the number of rows and columns. If size is 3, then row = 3, col = 3.
 */
matrix::matrix(const size_t &size) : rows(size), cols(size)
{
    if (size == 0)
        throw wrong_parameter{};
    elements = vector<double>(rows * cols);
}
matrix::matrix(const size_t &size, const vector<double> &input_vector)
    : rows(size), cols(size), elements(input_vector)
{
    if (rows == 0 or cols == 0)
        throw wrong_parameter{};
    if (input_vector.size() != rows * cols)
        throw initializer_wrong_size{};
}
/**
 * @brief Implementation of the constructor for a non-sqaure matrix
 * 
 * @param row The number of rows
 * @param col the number of columns
 */
matrix::matrix(const size_t &row, const size_t &col)
    : rows(row), cols(col)
{
    if (rows == 0 or cols == 0)
        throw wrong_parameter{};
    elements = vector<double>(rows * cols);
}
matrix::matrix(const size_t &row, const size_t &col, const vector<double> &input_vector)
    : rows(row), cols(col), elements(input_vector)
{
    if (rows == 0 or cols == 0)
        throw wrong_parameter{};
    if (input_vector.size() != rows * cols)
        throw wrong_parameter{};
}
/**
 * @brief Returns the rows and columns of the matrix
 * 
 * @return size_t The value
 */
size_t matrix::get_rows() const
{
    return rows;
}
size_t matrix::get_cols() const
{
    return cols;
}

/**
 * @brief Calculates the norm of a vector
 * 
 * @param u The vector
 * @return double The norm
 */
//L2 norm of a vector
double l2_norm(vector<double> const &u)
{
    double result = 0.;
    for (size_t i = 0; i < u.size(); ++i)
    {
        result += u[i] * u[i];
    }
    return sqrt(result);
}
/**
 * @brief Implementation of the overloaded operator() for accessing the elements of a non-square matrix
 * 
 * @param row The number of rows
 * @param col The number of columns
 * @return double& 
 */
double &matrix::operator()(const size_t &row, const size_t &col)
{
    return elements[(cols * row) + col];
}

double matrix::operator()(const size_t &row, const size_t &col) const
{
    return elements[(cols * row) + col];
}
/**
 * @brief Implementation of the overloaded operator() for accessing the elements of a square matrix
 * 
 * @param i Row
 * @return double& A reference to the elements 
 */
double &matrix::operator()(size_t &i)
{
    i -= 1;
    size_t row, col;
    row = i % cols;
    col = i / cols;
    return elements[i * row * col];
}
/**
 * @brief Implementation of the overloaded operator << for printing to stream
 * 
 * @param output Matrix
 * @param A Matrix
 * @return ostream& Reference to the matrix
 */
ostream &operator<<(ostream &output, const matrix &A)
{
    size_t rows{A.get_rows()}, cols{A.get_cols()};
    for (size_t i{0}; i < rows; i++)
    {
        output << "( ";
        for (size_t j{0}; j < cols; j++)
            output << A(i, j) << '\t';
        output << ")\n";
    }
    output << '\n';
    return output;
}
/**
 * @brief Implementation of the overloaded operator + for matrix summation
 * 
 * @param A Matrix
 * @param B Matrix
 * @return matrix The summation of both matrices
 */
matrix operator+(const matrix &A, const matrix &B)
{
    size_t rows{A.get_rows()}, cols{A.get_cols()};
    if ((rows != B.get_rows()) or (cols != B.get_cols()))
        throw matrix::incompatible_sizes_add{};
    matrix c(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            c(i, j) = A(i, j) + B(i, j);
    return c;
}
matrix operator+=(matrix &A, const matrix &B)
{
    A = A + B;
    return A;
}
/**
 * @brief Implementation of the overloaded operator - for negating a matrix
 * 
 * @param B Matrix
 * @return matrix A negative matrix
 */
matrix operator-(const matrix &B)
{
    size_t rows{B.get_rows()}, cols{B.get_cols()};
    matrix c(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            c(i, j) = -B(i, j);
    return c;
}
/**
 * @brief Implementation of the overloaded operator - for subtraction between matrices
 * 
 * @param A Matrix
 * @param B Matrix
 * @return matrix The result from the subtraction 
 */
matrix operator-(const matrix &A, const matrix &B)
{
    size_t rows{A.get_rows()}, cols{B.get_cols()};
    if ((rows != B.get_rows()) or (cols != B.get_cols()))
        throw matrix::incompatible_sizes_add{};
    matrix c(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            c(i, j) = A(i, j) - B(i, j);
    return c;
}
matrix operator-=(matrix &A, const matrix &B)
{
    A = A - B;
    return A;
}
/**
 * @brief Implementation of the overloaded operator* for matrix multiplication
 * 
 * @param A Matrix
 * @param B Matrix
 * @return matrix The result of the multiplication betwen matrices
 */
matrix operator*(const matrix &A, const matrix &B)
{
    size_t a_rows{A.get_rows()}, a_cols{A.get_cols()};
    size_t b_rows{B.get_rows()}, b_cols{B.get_cols()};
    if (a_cols != b_rows)
        throw matrix::incompatible_sizes_multiply{};
    matrix c(a_rows, b_cols);
    for (size_t i{0}; i < a_rows; i++)
        for (size_t j{0}; j < b_cols; j++)
            for (size_t k{0}; k < a_cols; k++)
                c(i, j) += A(i, k) * B(k, j);
    return c;
}
/**
 * @brief Implementation of the overloaded operator* for matrix-scalar  multiplication
 * 
 * @param c scalar
 * @param A Matrix
 * @return matrix The result of the multiplication betwen matrix-scalar
 */
matrix operator*(const double &c, const matrix &A)
{
    size_t rows{A.get_rows()}, cols{A.get_cols()};
    matrix s(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            s(i, j) = c * A(i, j);
    return s;
}
matrix operator*(const matrix &m, const double &c)
{
    return c * m;
}
/**
 * @brief Implementation of the overloaded operator* for matrix-vector multiplication
 * 
 * @param A Matrix
 * @param x Vector
 * @return Vector The result of the multiplication betwen matrix-vector
 */
vector<double> operator*(const matrix &A, const vector<double> &x)
{
    if (A.cols != x.size())
    {
        throw matrix::wrong_parameter{};
    }
    vector<double> output(A.rows);
    for (size_t i = 0; i < A.rows; i++)
    {
        for (size_t j = 0; j < x.size(); j++)
        {
            output[i] += A(i, j) * x[j];
        }
    }
    return output;
}
/**
 * @brief 
 * 
 * @param x Vector
 * @param A Matrix
 * @return Vector The result of the multiplication betwen vector-matrix
 */
vector<double> operator*(const vector<double> &x, const matrix &A)
{
    if (x.size() != A.rows)
    {
        throw matrix::wrong_parameter{};
    }
    vector<double> output(A.cols);

    for (size_t i = 0; i < A.cols; i++)
    {
        for (size_t j = 0; j < A.rows; j++)
        {
            output[i] += x[j] * A(i, j);
        }
    }
    return output;
}
/**
 * @brief Implementation of the overloaded operator/ for matrix-scalar division
 * 
 * @param A Matrix
 * @param c Scalar
 * @return matrix The result of the division betwen vector-matrix
 */
matrix operator/(const matrix &A, const double &c)
{
    if (c == 0.0)
    {
        throw matrix::wrong_parameter{};
    }
    double d = 1.0 / c;

    return A * d;
}
/**
 * @brief Implementation of the overloaded operator<< for printing a vector to stream
 * 
 * @param output Vector
 * @param v Vector
 * @return ostream& Reference to the vector
 */
ostream &operator<<(ostream &output, const vector<double> &v)
{
    size_t s{v.size() - 1};
    output << "(";
    for (size_t i{0}; i < s; i++)
        output << v[i] << ", ";
    output << v[s] << ")\n";
    return output;
}
/**
 * @brief Implementation of the overloaded operator+ for addition between vectors
 * 
 * @param v1 Vector
 * @param v2 Vector
 * @return vector<double> The result of the addition between vectors
 */
vector<double> operator+(const vector<double> &v1, const vector<double> &v2)
{
    size_t s{v1.size()};
    if (s != v2.size())
        throw matrix::incompatible_vectors{};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = v1[i] + v2[i];
    return u;
}
vector<double> operator+=(vector<double> &v1, const vector<double> &v2)
{
    v1 = v1 + v2;
    return v1;
}
/**
 * @brief Implementation of the overloaded operator- for negating a vector
 * 
 * @param v Vector
 * @return Vector A negative vector
 */
vector<double> operator-(const vector<double> &v)
{
    size_t s{v.size()};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = -v[i];
    return u;
}
/**
 * @brief Implementation of the overloaded operator- for subtraction between vectors
 * 
 * @param v1 Vector
 * @param v2 Vector
 * @return Vector The result of the subtraction between vectors
 */
vector<double> operator-(const vector<double> &v1, const vector<double> &v2)
{
    size_t s{v1.size()};
    if (s != v2.size())
        throw matrix::incompatible_vectors{};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = v1[i] - v2[i];
    return u;
}
vector<double> operator-=(vector<double> &v1, const vector<double> &v2)
{
    v1 = v1 - v2;
    return v1;
}
/**
 * @brief Implementation of the overloaded operator* for multiplication between vectors
 * 
 * @param v1 Vector
 * @param v2 Vector
 * @return double The result of the multiplication between vectors
 */
double operator*(const vector<double> &v1, const vector<double> &v2)
{
    size_t s{v1.size()};
    if (s != v2.size())
        throw matrix::incompatible_vectors{};
    double p{0};
    for (size_t i{0}; i < s; i++)
        p += v1[i] * v2[i];
    return p;
}
/**
 * @brief Implementation of the overloaded operator* for multiplication between scalar-vectors
 * @param x Scalar
 * @param v Vector
 * @return vector<double> The result of the multiplication between scalar-vectors
 */
vector<double> operator*(const double &x, const vector<double> &v)
{
    size_t s{v.size()};
    vector<double> u(s);
    for (size_t i{0}; i < s; i++)
        u[i] = x * v[i];
    return u;
}
vector<double> operator*(const vector<double> &v, const double &x)
{
    return x * v;
}
/**
 * @brief Implementation of the overloaded operator/ for division between vector-scalar
 * 
 * @param v vector
 * @param x scalar
 * @return vector<double> The result of the multiplication between scalar-vectors
 */
vector<double> operator/(const vector<double> &v, const double &x)
{
    if (x == 0)
    {
        throw matrix::incompatible_vectors{};
    }
    size_t s{v.size()};
    vector<double> u(s);
    double d = 1.0 / x;
    for (size_t i{0}; i < s; i++)
        u[i] = v[i] * d;
    return u;
}
/**
 * @brief Implementation of the biconjugate gardient stabilized method
 * 
 * @param A Matrix
 * @param b vector
 * @param verbose 
 * @param TOL tolerance
 * @param maxit Maximum iteration
 * @return vector<double> 
 */
vector<double> bicgstab(const matrix &A, const vector<double> &b, bool verbose, double TOL, int maxit)
{
    vector<double> x0(b.size());
    return bicgstab(A, b, x0, verbose, TOL, maxit);
}

vector<double> bicgstab(const matrix &A, const vector<double> &b, const vector<double> &x0, bool verbose, double TOL, int maxit)
{

    if (A.rows != A.cols)
    {
        throw matrix::wrong_parameter{};
    }
    if (A.cols != b.size())
    {
        throw matrix::wrong_parameter{};
    }
    if (maxit < 1)
    {
        throw matrix::wrong_parameter{};
    }
    if (TOL <= 0)
    {
        throw matrix::wrong_parameter{};
    }
    /**
 * @brief Definition of all parameters used in the algorithm
 */
    int n = 1;
    double alpha, beta, rho, rho_new, omega, omega_new, rtildeV, residual1, residual2, residual;
    vector<double> xn(x0), x_new(b.size());
    vector<double> r_initial(b.size()), r_new(b.size());
    vector<double> P_new(b.size()), P_old(b.size());
    vector<double> V_new(b.size());
    vector<double> h(b.size()), S(b.size()), t(b.size());
    /**
 * @brief Initial guesses
 */
    r_initial = b - A * xn;
    vector<double> rtilde(r_initial);
    vector<double> P_initial(xn);
    vector<double> V_initial(xn);
    rho = 1;
    alpha = 1;
    omega = 1;
    /**
 * @brief While loop
 */
    while (n < maxit)
    {
        rho_new = rtilde * r_initial;
        beta = (rho_new / rho) * (alpha / omega);
        P_new = r_initial + beta * (P_initial - (omega * V_initial));
        V_new = A * P_new;
        rtildeV = rtilde * V_new;
        alpha = rho_new / rtildeV;
        h = xn + (alpha * P_new);
        residual1 = (l2_norm(b - (A * xn))) / l2_norm(b);
        if (residual1 < TOL)
        {
            x_new = h;
            break;
        }
        S = r_initial - (alpha * V_new);
        t = A * S;
        omega_new = (t * S) / (t * t);
        x_new = h + (omega_new * S);
        residual2 = l2_norm(S) / l2_norm(b);
        if (residual2 < TOL)
        {

            break;
        }
        r_new = S - (omega_new * t);
        residual = (r_new * r_new);
        /**
 * @brief Iterates updates
 */
        n++;
        xn = x_new;
        r_initial = r_new;
        rho = rho_new;
        V_initial = V_new;
        P_initial = P_new;
        omega = omega_new;
    }
    /**
 * @brief Displays information about the results from the solver
 */
    if (verbose)
    {
        if (n > maxit)
        {
            cout << "BICGStab terminated after reaching the maximum iteration,\n";
            cout << "Residual: " << residual << "\n";
            cout << "Approximate value:" << x_new << "\n";
        }
        else
        {
            cout << "BICGStab converged to tolerance after " << n << " iteration(s) ";
            cout << "Residual: " << residual << "\n";
            cout << "Approximate value: " << x_new << "\n";
        }
    }
    return x_new;
}
