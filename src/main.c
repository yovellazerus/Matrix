
#include "Matrix.h"

int main(int argc, char* argv[]){

    double carr1[] = {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
    };
    Matrix* m1 = Matrix_create(carr1, 3, 3);
    Matrix* m1_copy = Matrix_copy(m1);
    Matrix* m1_t = Matrix_transpose(m1);
    MATRIX_DUMP(m1, NULL);
    MATRIX_DUMP(m1_copy, NULL);
    MATRIX_DUMP(m1_t, NULL);
    Matrix_destroy(m1);
    Matrix_destroy(m1_copy);
    Matrix_destroy(m1_t);

    double carr_row_vec[] = {
        1, 2, 3, 4,
    };
    Matrix* row_vec = Matrix_create(carr_row_vec, 1, 4);
    Matrix* row_vec_t = Matrix_transpose(row_vec);
    MATRIX_DUMP(row_vec, NULL);
    MATRIX_DUMP(row_vec_t, NULL);
    Matrix_destroy(row_vec);
    Matrix_destroy(row_vec_t);

    double carr_coll_vec[] = {
        1,
        2,
        3,
        4,
    };
    Matrix* coll_vec = Matrix_create(carr_coll_vec, 4, 1);
    Matrix* coll_vec_t = Matrix_transpose(coll_vec);
    MATRIX_DUMP(coll_vec, NULL);
    MATRIX_DUMP(coll_vec_t, NULL);
    Matrix_destroy(coll_vec);
    Matrix_destroy(coll_vec_t);

    Matrix* I = Matrix_I(4);
    MATRIX_DUMP(I, NULL);
    Matrix_destroy(I);

    Matrix* m_scalar = Matrix_scalar(4, 3.14);
    MATRIX_DUMP(m_scalar, NULL);
    Matrix_destroy(m_scalar);

    Matrix* random = Matrix_noise(3, 4, -3.0, 4.0);
    MATRIX_DUMP(random, NULL);
    Matrix_destroy(random);

    double aarr[] = {
        1, 2,
        3, 4,
    };
    double barr[] = {
        1, 4,
        -1, 7,
    };
    Matrix* a = Matrix_create(aarr, 2, 2);
    Matrix* b = Matrix_create(barr, 2, 2);
    Matrix* a_plus_b = Matrix_add(a, b);
    Matrix* a_minos_b = Matrix_sub(a, b);
    Matrix* ab = Matrix_dot(a, b);
    MATRIX_DUMP(a_plus_b, NULL);
    MATRIX_DUMP(a_minos_b, NULL);
    MATRIX_DUMP(ab, NULL);
    Matrix_destroy(a);
    Matrix_destroy(b);
    Matrix_destroy(a_plus_b);
    Matrix_destroy(a_minos_b);
    Matrix_destroy(ab);

    double marr[] = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16, 
    };
    Matrix* mm = Matrix_create(marr, 4, 4);
    Matrix* minor23 = Matrix_minor(mm, 2, 3);
    MATRIX_DUMP(minor23, NULL);
    Matrix_destroy(mm);
    Matrix_destroy(minor23);

    double arr_det[] = {
        1, 2, 1,
        1, 1, 2,
        2, 1, 1,
    };
    Matrix* m_det = Matrix_create(arr_det, 3, 3);
    double det = Matrix_det(m_det);
    printf("\n%lf\n", det);
    Matrix* m_det_inv = Matrix_inv(m_det);
    Matrix* m_det_again = Matrix_dot(m_det, m_det_inv);
    MATRIX_DUMP(m_det_inv, NULL);
    MATRIX_DUMP(m_det_again, NULL);
    Matrix_destroy(m_det);
    Matrix_destroy(m_det_inv);
    Matrix_destroy(m_det_again);
    
    return 0;
}