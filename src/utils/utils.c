#include "utils.h"

void print_mpi(int count, ...) {
    int i;
    va_list ap;
    gcry_error_t err;
    unsigned char *buffer;

    va_start(ap, count);

    for (i = 0; i < count; i++) {
        err = gcry_mpi_aprint(GCRYMPI_FMT_HEX, &buffer, NULL, va_arg(ap, gcry_mpi_t));
        if (err != 0) {
            printf(" error: %d\n", err);
            return;
        }
        printf("%d) %s\n",i,  buffer);
        gcry_free(buffer);
    }
    va_end(ap);
}
