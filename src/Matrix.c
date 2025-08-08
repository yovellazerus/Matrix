
#include "Matrix.h"

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

void Matrix_error(const char *msg)
{
    if(!msg){
        msg = "Unknown error";
    }
    fprintf(stderr, "ERROR: %s.\n", msg);
    exit(1);
}

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
        Matrix_error("Indexing error, index out of range");
        return NULL;
    }
    if(!mat){
        return NULL;
    }
    return &(mat->data[(i - 1) * (mat->colls) + (j - 1)]);
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
            fprintf(file, "\t%.4lf", *Matrix_at(mat, i, j));
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
        Matrix_error("Dimension error, Matrix dimensions not suitable for adding and subtracting");
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

Matrix *Matrix_dot(Matrix *a, Matrix *b)
{
    if(!a || !b){
        return NULL;
    }
    if(a->colls != b->rows){
        Matrix_error("Dimension error, Matrix dimensions not suitable for dot multiplication");
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
        Matrix_error("Indexing error, index out of range");
        return NULL;
    }
    if(mat->rows != mat->colls){
        Matrix_error("minor is not define for non square matrixes");
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

Matrix *Matrix_adj(Matrix *mat)
{
    if(!mat){
        return NULL;
    }
    if(mat->rows != mat-> rows){
        Matrix_error("adj is not define for non square matrixes");
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
        Matrix_error("determinant is not define for non square matrixes");
        return 0.0;
    }
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

// Very inefficient!
Matrix *Matrix_inv(Matrix *mat)
{
    if(!mat){
        return NULL;
    }
    if(mat->rows != mat->colls){
        Matrix_error("inverse is not define for non square matrixes");
        return NULL;
    }
    double det = Matrix_det(mat);
    if(det == 0.0){
        Matrix_error("mathematical error, attempting to invert an invertible matrix");
        return NULL;
    }
    Matrix* adj = Matrix_adj(mat);
    Matrix* res = Matrix_scale(adj, 1.0 / det);
    Matrix_destroy(adj);
    return res;
}
