
#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define MATRIX_DUMP(mat, file) (Matrix_dump(mat, file, #mat))

typedef struct Matrix_t Matrix;

void Matrix_error(const char* msg);

Matrix* Matrix_create(double* carr, size_t rows, size_t colls);
void Matrix_destroy(Matrix* mat);
Matrix* Matrix_copy(Matrix* mat);

Matrix* Matrix_I(size_t n);
Matrix* Matrix_scalar(size_t n, double alpha);
Matrix* Matrix_noise(size_t rows, size_t colls, double from, double to);

size_t Matrix_getRows(Matrix* mat);
size_t Matrix_getColls(Matrix* mat);
double* Matrix_at(Matrix* mat, size_t i, size_t j);
void Matrix_dump(Matrix* mat, FILE* file, const char* name);

Matrix* Matrix_add(Matrix* a, Matrix* b);
Matrix* Matrix_sub(Matrix* a, Matrix* b);
Matrix* Matrix_scale(Matrix* mat, double alpha);
Matrix* Matrix_dot(Matrix* a, Matrix* b);
Matrix* Matrix_transpose(Matrix* mat);
Matrix* Matrix_minor(Matrix* mat, size_t i, size_t j);

Matrix* Matrix_adj(Matrix* mat);
double Matrix_det(Matrix* mat);
Matrix* Matrix_inv(Matrix* mat);

#endif // MATRIX_H_
