
#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define MATRIX_DUMP(mat, file) (Matrix_dump(mat, file, #mat))
#define MATRIX_ERROR(msg) do {                                                                    \
    set_color(COLOR_RED, stderr);                                                                 \
    fprintf(stderr, "ERROR: %s in function: `%s()` in file: `%s` in line number: %d.\n",         \
        msg, __func__, __FILE__, __LINE__);                                                       \
    set_color(COLOR_RESET, stderr);                                                               \
    exit(1);                                                                                      \
} while(false)     
#define INF (1.0/0.0)               

typedef struct Matrix_t Matrix;

// Printing and debuging
void Matrix_dump(Matrix* mat, FILE* file, const char* name);

// Memory memory management
Matrix* Matrix_create(double* carr, size_t rows, size_t colls);
void Matrix_destroy(Matrix* mat);
Matrix* Matrix_copy(Matrix* mat);

// Creating special matrices
Matrix* Matrix_I(size_t n);
Matrix* Matrix_scalar(size_t n, double alpha);
Matrix* Matrix_dig(Matrix* vec);
Matrix* Matrix_noise(size_t rows, size_t colls, double from, double to);

// Get's and set's
size_t Matrix_getRows(Matrix* mat);
size_t Matrix_getColls(Matrix* mat);
double* Matrix_at(Matrix* mat, size_t i, size_t j);
Matrix* Matrix_getDig(Matrix* mat);

// Row operations
void Matrix_swap_rows(Matrix* mat, size_t i, size_t j);
void Matrix_scale_row(Matrix* mat, size_t row, double alpha);
void Matrix_add_row_to_row(Matrix* mat, size_t from, size_t to);

// matrix operations
void Matrix_scale_self(Matrix* mat, double alpha);

// Basic operations and arithmetic
Matrix* Matrix_add(Matrix* a, Matrix* b);
Matrix* Matrix_sub(Matrix* a, Matrix* b);
Matrix* Matrix_scale(Matrix* mat, double alpha);
Matrix* Matrix_dot(Matrix* a, Matrix* b);
Matrix* Matrix_transpose(Matrix* mat);
Matrix* Matrix_minor(Matrix* mat, size_t i, size_t j);
double Matrix_frob(Matrix* mat);

// Advanced operations
Matrix* Matrix_adj(Matrix* mat);
double Matrix_det(Matrix* mat);
Matrix* Matrix_inv(Matrix* mat);
double Matrix_gaussian_elimination(Matrix* mat);

// Eigenvalues and eigenvectors
double Matrix_smallest_eigenvalue(Matrix* mat);
double Matrix_greatest_eigenvalue(Matrix* mat);
Matrix* Matrix_eigenvalues(Matrix* mat, double lambda_min); // TODO: find `lambda_min` using: `Matrix_smallest_eigenvalue` 

#endif // MATRIX_H_
