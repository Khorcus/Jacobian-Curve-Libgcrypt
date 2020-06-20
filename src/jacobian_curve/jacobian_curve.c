#include "jacobian_curve.h"

void init_jacobian_curve(jacobian_curve *curve, gcry_mpi_t p, gcry_mpi_t q, gcry_mpi_t t,
                         gcry_mpi_t a, gcry_mpi_t x_base, gcry_mpi_t y_base) {
    // TODO: optimize number of buffers
    gcry_mpi_t zero, two, three, four, sixteen, inv_four, inv_sixteen;
    gcry_mpi_t buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8, buf9, buf10, buf11, buf12, buf13, buf14, buf15, buf16;

    zero = gcry_mpi_new(0);
    inv_four = gcry_mpi_new(0);
    inv_sixteen = gcry_mpi_new(0);
    buf1 = gcry_mpi_new(0);
    buf2 = gcry_mpi_new(0);
    buf3 = gcry_mpi_new(0);
    buf4 = gcry_mpi_new(0);
    buf5 = gcry_mpi_new(0);
    buf6 = gcry_mpi_new(0);
    buf7 = gcry_mpi_new(0);
    buf8 = gcry_mpi_new(0);
    buf9 = gcry_mpi_new(0);
    buf10 = gcry_mpi_new(0);
    buf11 = gcry_mpi_new(0);
    buf12 = gcry_mpi_new(0);
    buf13 = gcry_mpi_new(0);
    buf14 = gcry_mpi_new(0);
    buf15 = gcry_mpi_new(0);
    buf16 = gcry_mpi_new(0);

    gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, NULL);
    gcry_mpi_scan(&three, GCRYMPI_FMT_HEX, "3", 0, NULL);
    gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, NULL);
    gcry_mpi_scan(&sixteen, GCRYMPI_FMT_HEX, "10", 0, NULL);

    gcry_mpi_invm(inv_four, four, p);
    gcry_mpi_invm(inv_sixteen, sixteen, p);

    // Initialize e
    gcry_mpi_mulm(buf1, t, t, p); // buf1 = t^2
    gcry_mpi_mulm(buf2, buf1, three, p); // buf2 = 3*(t^2)
    gcry_mpi_mulm(buf3, a, four, p); // buf3 = 4*a
    gcry_mpi_addm(buf4, buf2, buf3, p); // buf4 = 3*(t^2) + 4*a
    gcry_mpi_subm(buf5, zero, buf4, p); // buf5 =  -(3*(t^2) + 4*a)
    gcry_mpi_mulm(buf6, buf5, inv_sixteen, p); // buf6 = (-(3*(t^2) + 4*a)) / 16

    // initialize d
    gcry_mpi_mulm(buf7, t, three, p); // buf7 = 3*t
    gcry_mpi_mulm(buf8, buf7, inv_four, p); // buf8 = (3*t)/4

    // Initialize x_base
    gcry_mpi_subm(buf9, x_base, t, p); // buf9 = x_base - t
    gcry_mpi_mulm(buf10, buf9, two, p); // buf10 = 2*(x_base - t)

    // Initialize y_base
    gcry_mpi_mulm(buf11, x_base, two, p); // buf11 = 2*x_base
    gcry_mpi_addm(buf12, buf11, t, p); // buf12 = 2*x_base + t
    gcry_mpi_mulm(buf13, buf9, buf9, p); // buf13 = (x_base - t)^2
    gcry_mpi_mulm(buf14, buf12, buf13, p); // buf14 = (2*x_base + t)*(x_base - t)^2
    gcry_mpi_mulm(buf15, y_base, y_base, p); // buf15 = y_base^2
    gcry_mpi_subm(buf16, buf14, buf15, p); // buf16 = (2*x_base + t)*(x_base - t)^2 - (y_base^2)

    curve->e = buf6;
    curve->d = buf8;
    curve->x_base = buf10;
    curve->y_base = buf16;
    curve->z_base = mpi_copy(y_base);


    curve->p = p;
    curve->q = q;
    curve->t = t;

    gcry_mpi_release(zero);
    gcry_mpi_release(two);
    gcry_mpi_release(three);
    gcry_mpi_release(four);
    gcry_mpi_release(sixteen);
    gcry_mpi_release(inv_four);
    gcry_mpi_release(inv_sixteen);
    gcry_mpi_release(buf1);
    gcry_mpi_release(buf2);
    gcry_mpi_release(buf3);
    gcry_mpi_release(buf4);
    gcry_mpi_release(buf5);
    gcry_mpi_release(buf7);
    gcry_mpi_release(buf9);
    gcry_mpi_release(buf11);
    gcry_mpi_release(buf12);
    gcry_mpi_release(buf13);
    gcry_mpi_release(buf14);
    gcry_mpi_release(buf15);
}

void free_jacobian_curve(jacobian_curve *curve) {
    gcry_mpi_release(curve->e);
    gcry_mpi_release(curve->d);
    gcry_mpi_release(curve->p);
    gcry_mpi_release(curve->q);
    gcry_mpi_release(curve->x_base);
    gcry_mpi_release(curve->y_base);
    gcry_mpi_release(curve->z_base);
}
