#include <stdio.h>
#include <gcrypt.h>
#include "point/point.h"
#include "jacobian_curve/jacobian_curve.h"
#include "params/params.h"
#include "utils/utils.h"

void print_delimiter() {
    printf("-----------------------------------------------------------------------------\n\n");
}

void print_equality(const point *p1, const point *p2, const jacobian_curve *curve) {
    printf(is_points_equals(p1, p2, curve) ? "YES\n" : "NO\n");
}

void print_is_on_curve(const point *p, const jacobian_curve *curve) {
    printf("Does the point lie on the curve? ");
    printf(is_point_on_curve(p, curve) ? "YES\n" : "NO\n");
}

void test() {
    int scalar, i;
    gcry_mpi_t p, q, t, a, x_base, y_base, one, two, q_plus_one, q_minus_one, mpi_scalar, k1, k2, k1_plus_k2;
    jacobian_curve curve;
    point point_e, point_b, point_e_plus_e, point_e_plus_b, point_b_plus_e, point_doubled_b, point_two_mult_b;
    point point_q_plus_one_mult_b, point_q_mult_b, point_minus_b, point_q_minus_one_mult_b;
    point point_scalar_mult_b_add, point_scalar_mult_b_montgomery;
    point point_k1_mult_b, point_k2_mult_b, point_k1_plus_k2_mult_b, point_k1_mult_b_plus_k2_mult_b;

    gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, p_str, 0, NULL);
    gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, q_str, 0, NULL);
    gcry_mpi_scan(&t, GCRYMPI_FMT_HEX, t_str, 0, NULL);
    gcry_mpi_scan(&a, GCRYMPI_FMT_HEX, a_str, 0, NULL);
    gcry_mpi_scan(&x_base, GCRYMPI_FMT_HEX, x_base_str, 0, NULL);
    gcry_mpi_scan(&y_base, GCRYMPI_FMT_HEX, y_base_str, 0, NULL);
    gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, NULL);
    gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, NULL);

    q_plus_one = gcry_mpi_new(0);
    q_minus_one = gcry_mpi_new(0);
    k1_plus_k2 = gcry_mpi_new(0);

    init_jacobian_curve(&curve, p, q, t, a, x_base, y_base);
    init_point_as_neutral(&point_e);
    init_point_by_mpi_numbers(&point_b, curve.x_base, curve.y_base, curve.z_base);

    gcry_mpi_addm(q_plus_one, curve.q, one, curve.p);
    gcry_mpi_subm(q_minus_one, curve.q, one, curve.p);

    add_point(&point_e_plus_e, &point_e, &point_e, &curve);
    add_point(&point_e_plus_b, &point_e, &point_b, &curve);
    add_point(&point_b_plus_e, &point_b, &point_e, &curve);
    double_point(&point_doubled_b, &point_b, &curve);
    mult_point_montgomery(&point_two_mult_b, &point_b, two, &curve);
    mult_point_montgomery(&point_q_plus_one_mult_b, &point_b, q_plus_one, &curve);
    mult_point_montgomery(&point_q_mult_b, &point_b, curve.q, &curve);
    negative_point(&point_minus_b, &point_b);
    mult_point_montgomery(&point_q_minus_one_mult_b, &point_b, q_minus_one, &curve);

    scalar = 957;
    gcry_mpi_scan(&mpi_scalar, GCRYMPI_FMT_HEX, "3BD", 0, NULL);

    init_point_as_neutral(&point_scalar_mult_b_add);
    for (i = 0; i < scalar; ++i) {
        add_point(&point_scalar_mult_b_add, &point_scalar_mult_b_add, &point_b, &curve);
    }
    mult_point_montgomery(&point_scalar_mult_b_montgomery, &point_b, mpi_scalar, &curve);

    gcry_mpi_scan(&k1, GCRYMPI_FMT_HEX, "AC024DE", 0, NULL);
    gcry_mpi_scan(&k2, GCRYMPI_FMT_HEX, "23924EFA", 0, NULL);
    gcry_mpi_addm(k1_plus_k2, k1, k2, curve.p);

    mult_point_montgomery(&point_k1_mult_b, &point_b, k1, &curve);
    mult_point_montgomery(&point_k2_mult_b, &point_b, k2, &curve);
    add_point(&point_k1_mult_b_plus_k2_mult_b, &point_k1_mult_b, &point_k2_mult_b, &curve);
    mult_point_montgomery(&point_k1_plus_k2_mult_b, &point_b, k1_plus_k2, &curve);


    print_delimiter();
    printf("Test 1\n\n");
    printf("Point E coordinates:\n");
    print_point_in_affine_coordinates(&point_e, &curve);
    print_is_on_curve(&point_e, &curve);
    printf("\n");

    printf("Point E + E coordinates:\n");
    print_point_in_affine_coordinates(&point_e_plus_e, &curve);
    print_is_on_curve(&point_e_plus_e, &curve);
    printf("\n");

    printf("Are points E and E + E equal? ");
    print_equality(&point_e, &point_e_plus_e, &curve);
    printf("\n");
    print_delimiter();

    printf("Test 2\n\n");
    printf("Point B coordinates:\n");
    print_point_in_affine_coordinates(&point_b, &curve);
    print_is_on_curve(&point_b, &curve);
    printf("\n");

    printf("Point E + B coordinates:\n");
    print_point_in_affine_coordinates(&point_e_plus_b, &curve);
    print_is_on_curve(&point_e_plus_b, &curve);
    printf("\n");

    printf("Are points B and E + B equal? ");
    print_equality(&point_b, &point_e_plus_b, &curve);
    printf("\n");

    printf("Point B + E coordinates:\n");
    print_point_in_affine_coordinates(&point_b_plus_e, &curve);
    print_is_on_curve(&point_b_plus_e, &curve);
    printf("\n");

    printf("Are points B and B + E equal? ");
    print_equality(&point_b, &point_b_plus_e, &curve);
    printf("\n");
    print_delimiter();

    printf("Test 3\n\n");
    printf("Point B + B coordinates:\n");
    print_point_in_affine_coordinates(&point_doubled_b, &curve);
    print_is_on_curve(&point_doubled_b, &curve);
    printf("\n");

    printf("Point 2*B coordinates:\n");
    print_point_in_affine_coordinates(&point_two_mult_b, &curve);
    print_is_on_curve(&point_two_mult_b, &curve);
    printf("\n");

    printf("Are points B + B and 2*B equal? ");
    print_equality(&point_doubled_b, &point_two_mult_b, &curve);
    printf("\n");
    print_delimiter();

    printf("Test 4\n\n");
    printf("Point (q + 1)*B coordinates:\n");
    print_point_in_affine_coordinates(&point_q_plus_one_mult_b, &curve);
    print_is_on_curve(&point_q_plus_one_mult_b, &curve);
    printf("\n");

    printf("Are points (q + 1)*B and B equal? ");
    print_equality(&point_q_plus_one_mult_b, &point_b, &curve);
    printf("\n");
    print_delimiter();

    printf("Test 5\n\n");
    printf("Point q*B coordinates:\n");
    print_point_in_affine_coordinates(&point_q_mult_b, &curve);
    print_is_on_curve(&point_q_mult_b, &curve);
    printf("\n");

    printf("Are points q*B and E equal? ");
    print_equality(&point_q_mult_b, &point_e, &curve);
    printf("\n");
    print_delimiter();

    printf("Test 6\n\n");
    printf("Point (q - 1)*B coordinates:\n");
    print_point_in_affine_coordinates(&point_q_minus_one_mult_b, &curve);
    print_is_on_curve(&point_q_minus_one_mult_b, &curve);
    printf("\n");

    printf("Point -B coordinates:\n");
    print_point_in_affine_coordinates(&point_minus_b, &curve);
    print_is_on_curve(&point_minus_b, &curve);
    printf("\n");

    printf("Are points (q - 1)*B and -B equal? ");
    print_equality(&point_e, &point_e_plus_e, &curve);
    printf("\n");
    print_delimiter();

    printf("Test 7\n\n");
    printf("Point B + B + ... + B (scalar times) coordinates:\n");
    print_point_in_affine_coordinates(&point_scalar_mult_b_add, &curve);
    print_is_on_curve(&point_scalar_mult_b_add, &curve);
    printf("\n");

    printf("Point scalar*B coordinates:\n");
    print_point_in_affine_coordinates(&point_scalar_mult_b_montgomery, &curve);
    print_is_on_curve(&point_scalar_mult_b_montgomery, &curve);
    printf("\n");

    printf("Are points B + B + ... + B (957 times) and 957*B equal? ");
    print_equality(&point_scalar_mult_b_add, &point_scalar_mult_b_montgomery, &curve);
    printf("\n");
    print_delimiter();

    printf("Test 8\n\n");
    printf("Point k1*B + k2*B coordinates:\n");
    print_point_in_affine_coordinates(&point_k1_mult_b_plus_k2_mult_b, &curve);
    print_is_on_curve(&point_k1_mult_b_plus_k2_mult_b, &curve);
    printf("\n");

    printf("Point (k1 + k2)*B:\n");
    print_point_in_affine_coordinates(&point_k1_plus_k2_mult_b, &curve);
    print_is_on_curve(&point_k1_plus_k2_mult_b, &curve);
    printf("\n");

    printf("Are points k1*B + k2*B and (k1 + k2)*B equal? ");
    print_equality(&point_k1_mult_b_plus_k2_mult_b, &point_k1_plus_k2_mult_b, &curve);
    printf("\n");
    print_delimiter();

    free_point(&point_e);
    free_point(&point_b);
    free_point(&point_e_plus_e);
    free_point(&point_e_plus_b);
    free_point(&point_b_plus_e);
    free_point(&point_doubled_b);
    free_point(&point_two_mult_b);
    free_point(&point_q_plus_one_mult_b);
    free_point(&point_q_mult_b);
    free_point(&point_minus_b);
    free_point(&point_q_minus_one_mult_b);
    free_point(&point_scalar_mult_b_add);
    free_point(&point_scalar_mult_b_montgomery);
    free_point(&point_k1_mult_b);
    free_point(&point_k2_mult_b);
    free_point(&point_k1_plus_k2_mult_b);
    free_point(&point_k1_mult_b_plus_k2_mult_b);
    free_jacobian_curve(&curve);

    gcry_mpi_release(a);
    gcry_mpi_release(x_base);
    gcry_mpi_release(y_base);
    gcry_mpi_release(one);
    gcry_mpi_release(two);
    gcry_mpi_release(q_plus_one);
    gcry_mpi_release(q_minus_one);
    gcry_mpi_release(mpi_scalar);
    gcry_mpi_release(k1);
    gcry_mpi_release(k2);
    gcry_mpi_release(k1_plus_k2);
}

int main() {
    test();
}
