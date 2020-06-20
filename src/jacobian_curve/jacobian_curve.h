#ifndef JACOBI_JACOBIAN_CURVE_H
#define JACOBI_JACOBIAN_CURVE_H

#include <gcrypt.h>

typedef struct jacobian_curve jacobian_curve;

struct jacobian_curve {
    gcry_mpi_t e, d, p, q, t, x_base, y_base, z_base;
};

void init_jacobian_curve(jacobian_curve *curve, gcry_mpi_t p, gcry_mpi_t q, gcry_mpi_t t,
                         gcry_mpi_t a, gcry_mpi_t x_base, gcry_mpi_t y_base);

void free_jacobian_curve(jacobian_curve *curve);

#endif //JACOBI_JACOBIAN_CURVE_H
