cmake_minimum_required(VERSION 3.15)
project(Jacobian-Curve-Libgcrypt C)

include_directories(/usr/local/Cellar/libgcrypt/1.8.5/include)
link_directories(/usr/local/Cellar/libgcrypt/1.8.5/lib)
link_libraries(gcrypt)

include_directories(/usr/local/Cellar/libgpg-error/1.38/include)
link_directories(/usr/local/Cellar/libgpg-error/1.38/lib)
link_libraries(gpg-error)

add_executable(Jacobian-Curve-Libgcrypt
        src/main.c
        src/point/point.c
        src/point/point.h
        src/params/params.h
        src/jacobian_curve/jacobian_curve.c
        src/jacobian_curve/jacobian_curve.h
        src/utils/utils.c
        src/utils/utils.h
        )