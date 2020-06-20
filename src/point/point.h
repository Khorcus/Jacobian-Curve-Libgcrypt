#ifndef JACOBI_POINT_H
#define JACOBI_POINT_H

#include "../jacobian_curve/jacobian_curve.h"
#include <stdio.h>
#include <gcrypt.h>

typedef struct point point;

struct point {
    gcry_mpi_t x, y, z;
};


void init_point_as_neutral(point *p);

void init_point_by_mpi_numbers(point *p, gcry_mpi_t x, gcry_mpi_t y, gcry_mpi_t z);

void init_point_by_chars(point *p, const char *x_str, const char *y_str, const char *z_str);

void copy_point(point *p_res, const point *p);

void add_point(point *p_res, const point *p1, const point *p2, const jacobian_curve *curve);

void double_point(point *p_res, const point *p, const jacobian_curve *curve);

void mult_point_montgomery(point *p_res, const point *p, gcry_mpi_t scalar, const jacobian_curve *curve);

void point_in_affine_coordinates(point *p_res, const point *p, const jacobian_curve *curve);

void negative_point(point *p_res, const point *p);

int is_points_equals(const point *p1, const point *p2, const jacobian_curve *curve);

int is_point_on_curve(const point *p, const jacobian_curve *curve);

void print_point(const point *p);

void print_point_in_affine_coordinates(const point *p, const jacobian_curve *curve);

void free_point(point *p);

#endif //JACOBI_POINT_H
