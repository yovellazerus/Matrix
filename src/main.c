
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

    double arr4[] = {
        1, 2, 3, 4,
        1, 0.5, 3, 4,
        1, 2, 3.14, 4,
        1, 2, 3, 4.76,
    };
    Matrix* m4 = Matrix_create(arr4, 4, 4);
    Matrix* dig = Matrix_getDig(m4);
    MATRIX_DUMP(dig, NULL);
    Matrix_destroy(m4);
    Matrix_destroy(dig);

    double arr5[] = {
        1, 2, 3, 4, 
    };
    Matrix* vec1 = Matrix_create(arr5, 1, 4);
    Matrix* m5 = Matrix_dig(vec1);
    MATRIX_DUMP(m5, NULL);
    Matrix_destroy(m5);
    Matrix_destroy(vec1);

    double arr6[] = {
        4,
        3,
        2,
        1, 
    };
    Matrix* vec2 = Matrix_create(arr6, 4, 1);
    Matrix* m6 = Matrix_dig(vec2);
    MATRIX_DUMP(m6, NULL);
    Matrix_destroy(m6);
    Matrix_destroy(vec2);

    Matrix* noise1 = Matrix_noise(5, 5, -2, 2);
    Matrix* noise1_inv = Matrix_inv(noise1);
    Matrix* noise1_again = Matrix_dot(noise1, noise1_inv);
    MATRIX_DUMP(noise1, NULL);
    MATRIX_DUMP(noise1_inv, NULL);
    MATRIX_DUMP(noise1_again, NULL);
    Matrix_destroy(noise1);
    Matrix_destroy(noise1_inv);
    Matrix_destroy(noise1_again);

    double arr7[] = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
    };
    Matrix* m_forSwap = Matrix_create(arr7, 3, 4);
    Matrix_scale_row(m_forSwap, 1, 5.0);
    Matrix_swap_rows(m_forSwap, 2, 3);
    Matrix_add_row_to_row(m_forSwap, 1, 3);
    MATRIX_DUMP(m_forSwap, NULL);
    Matrix_destroy(m_forSwap);

    Matrix* m_for_g = Matrix_noise(4, 4, 2.0, 4.0);
    Matrix* m_for_g_origen = Matrix_copy(m_for_g);
    MATRIX_DUMP(m_for_g_origen, NULL);
    Matrix_gaussian_elimination(m_for_g);
    MATRIX_DUMP(m_for_g, NULL);
    Matrix_destroy(m_for_g);
    Matrix_destroy(m_for_g_origen);

    double arr_p[] = {
        0, 1, 2,
        0, 0, 3,
        2, 3, 4,
    };
    Matrix* mp = Matrix_create(arr_p, 3, 3);
    Matrix_gaussian_elimination(mp);
    MATRIX_DUMP(mp, NULL);
    Matrix_destroy(mp);

    double arr_eig[] = {
    2.0, -1.0,  0.0,  0.0,
   -1.0,  2.0, -1.0,  0.0,
    0.0, -1.0,  2.0, -1.0,
    0.0,  0.0, -1.0,  2.0
};
    Matrix* m_eig = Matrix_create(arr_eig, 4, 4);
    Matrix* m_eig_dig = Matrix_eigenvalues(m_eig);
    MATRIX_DUMP(m_eig, NULL);
    MATRIX_DUMP(m_eig_dig, NULL);
    Matrix_destroy(m_eig_dig);
    printf("%lf\n", Matrix_greatest_eigenvalue(m_eig));
    printf("%lf\n", Matrix_smallest_eigenvalue(m_eig));
    Matrix_destroy(m_eig);
    
    return 0;
}












