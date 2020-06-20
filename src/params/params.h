#ifndef JACOBI_PARAMS_H
#define JACOBI_PARAMS_H

/**
 * p_str – Характеристика простого поля, над которым определяется эллиптическая кривая (модуль)
 * q_str – Порядок подгруппы простого порядка группы точек эллиптической кривой
 * t_str - (t,0) - точка второго порядка в форме Вейерштрасса
 * a_str - Коэффициент эллиптической кривой в форме Вейерштрасса
 * x_base_str, y_base_str - Координаты точки B
 * (порождающего элемента подгруппы простого порядка) на кривой в форме Вейерштрасса
 */


#define p_str      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97"
#define q_str      "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67"
#define t_str      "100FE73F595FF158E974B44D478D9588744FE5C192AC47EA63075DCE7A14AAA"
#define a_str      "C2173F1513981673AF4892C23035A27CE25E2013BF95AA33B22C656F277E7335"
#define x_base_str "91E38443A5E82C0D880923425712B2BB658B9196932E02C78B2582FE742DAA28"
#define y_base_str "32879423AB1A0375895786C4BB46E9565FDE0B5344766740AF268ADB32322E5C"

#endif //JACOBI_PARAMS_H
