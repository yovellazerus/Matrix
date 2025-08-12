
#include "Matrix.h"

typedef enum {
    COLOR_RESET = 0,
    COLOR_RED,
    COLOR_GREEN,
    COLOR_YELLOW,
    COLOR_BLUE,
    COLOR_MAGENTA,
    COLOR_CYAN,
    COLOR_WHITE
} Color;

static void set_color(Color color, FILE* file) {
    if(!file){
        file = stdout;
    }
    switch (color) {
        case COLOR_RED:     fprintf(file, "\033[0;31m"); break;
        case COLOR_GREEN:   fprintf(file, "\033[0;32m"); break;
        case COLOR_YELLOW:  fprintf(file, "\033[0;33m"); break;
        case COLOR_BLUE:    fprintf(file, "\033[0;34m"); break;
        case COLOR_MAGENTA: fprintf(file, "\033[0;35m"); break;
        case COLOR_CYAN:    fprintf(file, "\033[0;36m"); break;
        case COLOR_WHITE:   fprintf(file, "\033[0;37m"); break;
        case COLOR_RESET:
        default:            fprintf(file, "\033[0m");    break;
    }
}

static double randomInRange(double from, double to){
    if(from >= to){
        return 0.0;
    }
    double seg = (double) rand() / RAND_MAX;
    if(rand() % 2 == 0){
        seg *= to;
    }
    else{
        seg *= from;
    }
    return seg;
}

struct Matrix_t {
    double* data;
    size_t rows;
    size_t colls;
};

Matrix *Matrix_create(double *carr, size_t rows, size_t colls)
{
    if(rows == 0 || colls == 0){
        return NULL;
    }

    Matrix* res = (Matrix*)malloc(sizeof(*res));
    if(!res){
        return NULL;
    }

    double* new_arr = (double*)malloc(sizeof(*new_arr) * rows * colls);
    if(!new_arr){
        free(res);
        return NULL;
    }

    if(!carr){
        for(int i = 0; i < rows*colls; i++){
            new_arr[i] = 0.0;
        }
    }
    else{
        for(int i = 0; i < rows*colls; i++){
            new_arr[i] = carr[i];
        }
    }

    res->data = new_arr;
    res->rows = rows;
    res->colls = colls;

    return res;
}

void Matrix_destroy(Matrix *mat)
{
    if(!mat){
        return;
    }
    free(mat->data);
    free(mat);
}

Matrix *Matrix_copy(Matrix *mat)
{
    if(!mat){
        return NULL;
    }

    Matrix* res = (Matrix*)malloc(sizeof(*res));
    if(!res){
        return NULL;
    }

    double* new_arr = (double*)malloc(sizeof(*new_arr) * mat->rows * mat->colls);
    if(!new_arr){
        free(res);
        return NULL;
    }

    
    for(int i = 0; i < mat->rows*mat->colls; i++){
        new_arr[i] = mat->data[i];
    }
    
    res->data = new_arr;
    res->rows = mat->rows;
    res->colls = mat->colls;

    return res;
}

Matrix *Matrix_I(size_t n)
{
    return Matrix_scalar(n, 1.0);
}

Matrix *Matrix_scalar(size_t n, double alpha)
{
    if(n == 0){
        return NULL;
    }
    Matrix* res = Matrix_create(NULL, n, n);
    if(!res){
        return NULL;
    }
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            if(i == j) *Matrix_at(res, i, j) = alpha;
        }
    }
    return res;
}

Matrix *Matrix_dig(Matrix *vec)
{
    if(!vec){
        return NULL;
    }
    if(vec->rows != 1 && vec->colls != 1){
        MATRIX_ERROR("can not make a diagonal matrix from a non flat-matrix (non vector)");
        return NULL;
    }
    size_t n = 0;
    bool row_vec = false;
    if(vec->rows != 1){
        n = vec->rows;
    }
    else{
        n = vec->colls;
        row_vec = true;
    }
    Matrix* res = Matrix_create(NULL, n, n);
    if(!res){
        return NULL;
    }
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            if(i == j){
                if(row_vec) {
                    *Matrix_at(res, i, j) = *Matrix_at(vec, 1, j);
                }
                else{
                    *Matrix_at(res, i, j) = *Matrix_at(vec, i, 1);
                }
            } 
        }
    }
    return res;
}

Matrix *Matrix_noise(size_t rows, size_t colls, double from, double to)
{
    srand(time(NULL));
    if(rows == 0 || colls == 0 || from >= to){
        return NULL;
    }
    Matrix* res = Matrix_create(NULL, rows, colls);
    if(!res){
        return NULL;
    }
    for(int i = 1; i <= rows; i++){
        for(int j = 1; j <= colls; j++){
            *Matrix_at(res, i, j) = randomInRange(from, to);
        }
    }
    return res;
}

size_t Matrix_getRows(Matrix *mat)
{
    if(!mat){
        return 0;
    }
    return mat->rows;
}

size_t Matrix_getColls(Matrix *mat)
{
    if(!mat){
        return 0;
    }
    return mat->colls;
}

double *Matrix_at(Matrix *mat, size_t i, size_t j)
{
    if(i == 0 || i > mat->rows || j == 0 || j > mat-> colls){
        MATRIX_ERROR("Indexing error, index out of range");
        return NULL;
    }
    if(!mat){
        return NULL;
    }
    return &(mat->data[(i - 1) * (mat->colls) + (j - 1)]);
}

Matrix *Matrix_getDig(Matrix *mat)
{
    if(!mat){
        return NULL;
    }
    if(mat->rows != mat->colls){
        MATRIX_ERROR("can not get diagonal of non square matrix");
        return NULL;
    }
    Matrix* res = Matrix_create(NULL, 1, mat->colls);
    if(!res){
        return NULL;
    }
    for(int i = 1; i <= mat->rows; i++){
        for(int j = 1; j <= mat->colls; j++){
            if(i == j) *Matrix_at(res, 1, j) = *Matrix_at(mat, i, j);
        }
    }
    return res;
}

void Matrix_swap_rows(Matrix *mat, size_t i, size_t j)
{
    if(!mat){
        return;
    }
    if(i <= 0 || i > mat->rows || j <= 0 || j > mat->rows){
        MATRIX_ERROR("invalid selection of rows to replace");
        return;
    }

    double* row_i = (double*)malloc(sizeof(*row_i) * mat->colls);
    if(!row_i){
        return;
    }
    for(int k = 0; k < mat->colls; k++){
        row_i[k] = *Matrix_at(mat, i, k+1);
    }
    double* row_j = (double*)malloc(sizeof(*row_j) * mat->colls);
    if(!row_j){
        free(row_i);
        return;
    }
    for(int k = 0; k < mat->colls; k++){
        row_j[k] = *Matrix_at(mat, j, k+1);
    }
    for(int k = 0; k < mat->colls; k++){
        *Matrix_at(mat, j, k+1) = row_i[k];
    }
    for(int k = 0; k < mat->colls; k++){
        *Matrix_at(mat, i, k+1) = row_j[k];
    }
    free(row_i);
    free(row_j);
}

void Matrix_scale_row(Matrix *mat, size_t row, double alpha)
{
    if(!mat){
        return;
    }
    if(row <= 0 || row > mat->rows){
        MATRIX_ERROR("attempting to change a row that does not exist in the matrix");
        return;
    }
    for(int j = 0; j < mat->colls; j++){
        *Matrix_at(mat, row, j+1) *= alpha; 
    }
}

void Matrix_add_row_to_row(Matrix *mat, size_t from, size_t to)
{
    if(!mat){
        return;
    }
    if(from <= 0 || from > mat->rows || to <= 0 || to > mat->rows){
        MATRIX_ERROR("invalid selection of rows to add to each other");
        return;
    }
    for(int j = 1; j <= mat->colls; j++){
        *Matrix_at(mat, to, j) += *Matrix_at(mat, from, j);
    }
}

void Matrix_dump(Matrix *mat, FILE *file, const char *name)
{
    if(!mat){
        return;
    }
    if(!file){
        file = stdout;
    }
    
    if(name){
        fprintf(file, "%s = [", name);
    }
    else{
       fprintf(file, "["); 
    }

    for(int i = 1; i <= mat->rows; i++){
        fprintf(file, "\n");
        for(int j = 1; j <= mat->colls; j++){
            fprintf(file, "\t%.2lf", *Matrix_at(mat, i, j));
        }
    }

    fprintf(file, "\n]\n");

}

Matrix *Matrix_add(Matrix *a, Matrix *b)
{
    if(!a || !b){
        return NULL;
    }
    if(a->rows != b->rows || a->colls != b->colls){
        MATRIX_ERROR("Dimension error, Matrix dimensions not suitable for adding and subtracting");
        return NULL;
    }
    Matrix* res = Matrix_create(NULL, a->rows, a->colls);
    if(!res){
        return NULL;
    }
    for(int i = 1; i <= a->rows; i++){
        for(int j = 1; j <= a->colls; j++){
            *Matrix_at(res, i, j) = *Matrix_at(a, i, j) + *Matrix_at(b, i, j);
        }
    }
    return res;
}

Matrix *Matrix_sub(Matrix *a, Matrix *b)
{
    return Matrix_add(a, Matrix_scale(b, -1.0));
}

Matrix *Matrix_scale(Matrix *mat, double alpha)
{
    if(!mat){
        return NULL;
    }
    Matrix* res = Matrix_copy(mat);
    if(!res){
        return NULL;
    }
    for(int i = 1; i <= res->rows; i++){
        for(int j = 1; j <= res->colls; j++){
            *Matrix_at(res, i, j) *= alpha;
        }
    }
    return res;
}

void Matrix_scale_self(Matrix *mat, double alpha)
{
    if(!mat){
        return;
    }
    for(int i = 1; i <= mat->rows; i++){
        for(int j = 1; j <= mat->colls; j++){
            *Matrix_at(mat, i, j) *= alpha;
        }
    }
}

Matrix *Matrix_dot(Matrix *a, Matrix *b)
{
    if(!a || !b){
        return NULL;
    }
    if(a->colls != b->rows){
        MATRIX_ERROR("Dimension error, Matrix dimensions not suitable for dot multiplication");
        return NULL;
    }
    size_t n = a->rows;
    size_t m = a->colls;
    size_t p = b->colls;
    Matrix* ab = Matrix_create(NULL, n, p);
    if(!ab){
        return NULL;
    }
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= p; j++){
            for(int k = 1; k <= m; k++){
                *Matrix_at(ab, i, j) += *Matrix_at(a, i, k) * *Matrix_at(b, k, j);
            }
        }
    }
    return ab;
}

Matrix *Matrix_transpose(Matrix *mat)
{
    if(!mat){
        return NULL;
    }
    Matrix* res = Matrix_create(NULL, mat->colls, mat->rows);
    if(!res){
        return NULL;
    }
    for(int i = 1; i <= mat->rows; i++){
        for(int j = 1; j <= mat->colls; j++){
            *Matrix_at(res, j, i) = *Matrix_at(mat, i, j);
        }
    }
    return res;
}

Matrix *Matrix_minor(Matrix *mat, size_t i, size_t j)
{
    if(!mat){
        return NULL;
    }
    if(i == 0 || i > mat->rows || j == 0 || j > mat-> colls){
        MATRIX_ERROR("Indexing error, index out of range");
        return NULL;
    }
    if(mat->rows != mat->colls){
        MATRIX_ERROR("minor is not define for non square matrixes");
        return NULL;
    }
    size_t n = mat->rows;
    Matrix* res = Matrix_create(NULL, n - 1, n - 1);
    if(!res){
        return NULL;
    }

    size_t i_res = 0;
    size_t j_res = 0;
    for(int i_in = 1; i_in <= n; i_in++){
        i_res++;
        if(i_in == i) i_res--;
        for(int j_in = 1; j_in <= n; j_in++){
            j_res++;
            if(j_in == j) j_res--;
            if(i_in != i && j_in != j){
                *Matrix_at(res, i_res, j_res) = *Matrix_at(mat, i_in, j_in);
            }
        }
        j_res = 0;
    }

    return res;
}

double Matrix_frob(Matrix *mat)
{
    if(!mat){
        return 0.0;
    }
    double res = 0.0;
    for(int i = 1; i <= mat->rows; i++){
        for(int j = 1; j <= mat->colls; j++){
            res += (*Matrix_at(mat, i, j)) * (*Matrix_at(mat, i, j));
        }
    }
    return sqrt(res); // TODO: write myself
}

Matrix *Matrix_adj(Matrix *mat)
{
    if(!mat){
        return NULL;
    }
    if(mat->rows != mat-> rows){
        MATRIX_ERROR("adj is not define for non square matrixes");
        return NULL;
    }
    size_t n = mat->rows;
    Matrix* res = Matrix_create(NULL, n, n);
    if(!res){
        return NULL;
    }
    double sign = 1.0;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            Matrix* minor = Matrix_minor(mat, j, i);
            if((j + i) % 2 == 0){
                sign = 1.0;
            }
            else{
                sign = -1.0;
            }
            *Matrix_at(res, i, j) = sign * Matrix_det(minor);
            Matrix_destroy(minor);
        }
    }
    return res;
}

double Matrix_det(Matrix *mat)
{
    if(!mat){
        return 0.0;
    }
    if(mat->rows != mat->colls){
        MATRIX_ERROR("determinant is not define for non square matrixes");
        return 0.0;
    }

    // Gaussian method, time complexity: O(n^3)
    Matrix* copy = Matrix_copy(mat);
    if(!copy){
        return 0.0;
    }
    double det_factor = Matrix_gaussian_elimination(copy);
    for(size_t i = 1; i <= copy->rows; i++){
        for(size_t j = 1; j <= copy->colls; j++){
            if(i == j) det_factor *= *Matrix_at(copy, i, j);
        } 
    }
    Matrix_destroy(copy);
    return det_factor;

    // Laplace's method, time complexity: O(n!)
    size_t n = mat->rows;
    if(n == 1){
        return *Matrix_at(mat, 1, 1);
    }

    double res = 0.0;
    size_t i = 1; // for the first row
    double sign = 1.0;
    for(size_t j = 1; j <= n; j++){
        if((j + i) % 2 == 0){
            sign = 1.0;
        }
        else{
            sign = -1.0;
        }
        Matrix* minor = Matrix_minor(mat, i, j);
        res += sign * *Matrix_at(mat, i, j) * Matrix_det(minor);
        Matrix_destroy(minor);
    }
    return res;
}

// Very inefficient! TODO: to change by solving a system of equations. Now possible with matrix ranking option 
Matrix *Matrix_inv(Matrix *mat)
{
    if(!mat){
        return NULL;
    }
    if(mat->rows != mat->colls){
        MATRIX_ERROR("inverse is not define for non square matrixes");
        return NULL;
    }
    double det = Matrix_det(mat);
    if(det == 0.0){
        MATRIX_ERROR("mathematical error, attempting to invert an invertible matrix");
        return NULL;
    }
    Matrix* adj = Matrix_adj(mat);
    Matrix* res = Matrix_scale(adj, 1.0 / det);
    Matrix_destroy(adj);
    return res;
}

// TODO: Consider adding an option to rank for canonical form
double Matrix_gaussian_elimination(Matrix *mat)
{
    if(!mat){
        return 0.0;
    }
    double det_factor = 1.0;
    for(size_t coll = 1; coll < Matrix_getColls(mat); coll++){
            for(size_t row = coll+1; row <= Matrix_getRows(mat); row++){
            if(*Matrix_at(mat, row, coll) == 0.0) continue;
            // promotion
            if(*Matrix_at(mat, coll, coll) == 0.0){
                Matrix_swap_rows(mat, coll, row);
                det_factor *= -1;
                continue;
            }
            double alpha = (*Matrix_at(mat, row, coll) / *Matrix_at(mat, coll, coll));
            Matrix_scale_row(mat, coll, -alpha);
            Matrix_add_row_to_row(mat, coll, row);
            det_factor *= -1 / alpha;
        }
    }
    return det_factor;
}

// TODO:not working...
double Matrix_smallest_eigenvalue(Matrix *mat)
{
    if(!mat){
        return -INF;
    }
    if(mat->colls != mat->rows){
        MATRIX_ERROR("Eigenvalues are not defined for a non-square matrix");
        return -INF;
    }

    double mat_det = Matrix_det(mat);
    Matrix* mat_shifted = NULL;

    if(fabs(mat_det) < 1e-12){  // consider near-singular also
        double sigma = 1e-3;
        Matrix* I_sigma = Matrix_scalar(mat->rows, sigma);
        if(!I_sigma){
            return -INF;
        }
        mat_shifted = Matrix_sub(mat, I_sigma);
        Matrix_destroy(I_sigma);
        if(!mat_shifted){
            return -INF;
        }
    } else {
        mat_shifted = mat;
    }

    Matrix* mat_inv = Matrix_inv(mat_shifted);
    if(mat_shifted != mat){
        Matrix_destroy(mat_shifted); // free only if shifted
    }

    if(!mat_inv){
        return -INF;
    }

    double lambda_min_inv = Matrix_greatest_eigenvalue(mat_inv);
    Matrix_destroy(mat_inv);
    if(lambda_min_inv == INF){
        return -INF;
    }

    return 1.0 / lambda_min_inv;
}

double Matrix_greatest_eigenvalue(Matrix *mat)
{
    if(!mat){
        return INF;
    }
    if(mat->colls != mat->rows){
        MATRIX_ERROR("Eigenvalues are not defined for a non-square matrix");
        return INF;
    }

    /*

    do max_iterations times:
    v <- Av / ||Av||

    */

    Matrix* A = mat;
    size_t max_iterations = 100;
    Matrix* v = Matrix_noise(mat->rows, 1, -1.0, 1.0); // TODO: `from` and `to` and `max_iterations` need to be generalize
    for(size_t i = 0; i < 100; i++){
        Matrix* Av = Matrix_dot(A, v);
        if(!Av){
            Matrix_destroy(v);
            return INF;
        }
        double norm = Matrix_frob(Av);
        Matrix_scale_self(Av, 1.0 / norm);
        Matrix_destroy(v);
        v = Av;
    }

    /*
    
    lambda_max = vt_A_v / vt_v
    
    */

    if(!v){
        return INF;
    }
    Matrix* vt = Matrix_transpose(v);
    if(!vt){
        Matrix_destroy(v);
        return INF;
    }
    double lambda_max = INF;
    Matrix* Av = Matrix_dot(A, v);
    if(!Av){
        Matrix_destroy(v);
        Matrix_destroy(vt);
        return INF;
    }
    Matrix* vt_A = Matrix_dot(vt, A);
    if(!vt_A){
        Matrix_destroy(v);
        Matrix_destroy(vt);
        Matrix_destroy(Av);
        return INF;
    }
    Matrix* vt_A_v = Matrix_dot(vt_A, v);
    if(!vt_A_v){
        Matrix_destroy(v);
        Matrix_destroy(vt);
        Matrix_destroy(Av);
        Matrix_destroy(vt_A);
        return INF;
    }
    Matrix* vt_v = Matrix_dot(vt, v);
    if(!vt_v){
        Matrix_destroy(v);
        Matrix_destroy(vt);
        Matrix_destroy(Av);
        Matrix_destroy(vt_A);
        Matrix_destroy(vt_A_v);
        return INF;
    }
    lambda_max = (*Matrix_at(vt_A_v, 1, 1) / *Matrix_at(vt_v, 1, 1));
    Matrix_destroy(v);
    Matrix_destroy(vt);
    Matrix_destroy(Av);
    Matrix_destroy(vt_A);
    Matrix_destroy(vt_A_v);
    Matrix_destroy(vt_v);
    return lambda_max;
}

Matrix *Matrix_eigenvalues(Matrix *mat, double lambda_min)
{
    if(!mat){
        return NULL;
    }
    if(mat->colls != mat->rows){
        MATRIX_ERROR("Eigenvalues are not defined for a non-square matrix");
        return NULL;
    }

    Matrix* res = Matrix_create(NULL, mat->rows, mat->colls);
    if(!res){
        return NULL;
    }

    size_t i = 1;
    double delta = 1e-3;
    double from = lambda_min; // Matrix_smallest_eigenvalue(mat) - delta;
    double to = Matrix_greatest_eigenvalue(mat) + delta;
    double characteristic_polynomial_val_prev = 0.0;
    double characteristic_polynomial_val_curr = 0.0;
    for(double lambda = from; lambda <= to; lambda += delta){
        Matrix* ILambda = Matrix_scalar(mat->rows, lambda);
        if(!ILambda){
            return NULL;
        }
        Matrix* A_minus_ILambda = Matrix_sub(mat, ILambda);
        if(!A_minus_ILambda){
            Matrix_destroy(ILambda);
            return NULL;
        }
        characteristic_polynomial_val_curr = Matrix_det(A_minus_ILambda);
        if(characteristic_polynomial_val_prev * characteristic_polynomial_val_curr < 0.0){
            *Matrix_at(res, i, i) = lambda;
            i++;
        }
        characteristic_polynomial_val_prev = characteristic_polynomial_val_curr;
        Matrix_destroy(ILambda);
        Matrix_destroy(A_minus_ILambda);
    }
    return res;
}
