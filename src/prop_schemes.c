#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "prop_schemes.h"


// Expects wavefunction in local basis
void propagate_NISE(
    t_non *non,
    float *H,
    float *e, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci
) {
    float f;
    int index, N;
    int i;

    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;

    matrix_on_vector(H, cr, ci, N);
    // Exponentiate [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(e[i] * f);
        im_U[i] = -sin(e[i] * f);
    }

    // Multiply with matrix exponent
    vector_on_vector(re_U, im_U, cr, ci, N);

    trans_matrix_on_vector(H, cr, ci, N);
}

// Expects wavefunction in adiabatic or local basis
void propagate_nise_dba(
    t_non *non, 
    float *H_old, 
    float *H_new, 
    float *e, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci
) {
    float f;
    int index, N;
    float qc_ij, qc_ji, ediff;
    float* abs;
    int i, j, k;
    
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;
    abs = (float *)calloc(N, sizeof(float));
    
    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Average")) {
        matrix_on_vector(H_old, cr, ci, N);
    }

    // Compute absolute values of the wavefunction coeffs
    for (i = 0; i < N; i++) {
        abs[i] = sqrt(cr[i]*cr[i] + ci[i]*ci[i]);
    }

    // Exponentiate [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(e[i] * f);
        im_U[i] = -sin(e[i] * f);
    }

    // Compute unadjusted non-adiabatic coupling
    matrix_on_matrix(H_new, H_old, N);
    // Compute temperature adjustments
    thermal_correction(non, H_old, e, abs);
    // Exponentiate the non-adiabatic couplings
    matrix_exp(H_old, N);
    // Multiply with (real) non-adiabatic propagator
    trans_matrix_on_vector(H_old, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Average")) {
        trans_matrix_on_vector(H_new, cr, ci, N);
    }

    free(abs);
}

// Expects wavefunction in average eigenbasis or local basis
void propagate_nise_dbb(
    t_non *non, 
    float *H_avg,
    float *H_new, 
    float *e_avg,
    float *e,
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci
) {
    float f;
    int index, N, N2;
    float *ecopy;
    float *Hcopy;
    float *Hcc;
    float *Urcopy;
    float *Uicopy;
    float qc_ij, qc_ji, ediff;
    float* abs;
    int i, j, k;
    
    N = non->singles;
    N2 = N*N;
    f = non->deltat * icm2ifs * twoPi;

    abs = (float *) calloc(N, sizeof(float));
    Hcopy = (float *) calloc(N2, sizeof(float));
    Hcc = (float *) calloc(N2, sizeof(float));
    ecopy = (float *) calloc(N, sizeof(float));
    Urcopy = (float *) calloc(N, sizeof(float));
    Uicopy = (float *) calloc(N, sizeof(float));

    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Adiabatic")) {
        // Transfer local -> average basis
        trans_matrix_on_vector(H_avg, cr, ci, N);
    }

    // Compute absolute values of the wavefunction coeffs
    for (i = 0; i < N; i++) {
        abs[i] = sqrt(cr[i]*cr[i] + ci[i]*ci[i]);
    }
    
    clearvec(Hcc, N2);
    // Find site basis perturbation
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            // Diagonals
            Hcc[i + N*i] += e[j] * H_new[i + N*j] * H_new[i + N*j];
            // Off-diagonals
            for (k = 0; k < N; k++) {
                Hcc[i + N*j] += e[k] * H_new[i + N*k] * H_new[j + N*k];
            }
        }
    }

    // Find average eigenbasis perturbation
    copyvec(H_avg, Hcopy, N2);
    matrix_on_matrix(H_avg, Hcc, N);
    matrix_on_matrix(Hcc, Hcopy, N);
    free(Hcc);

    // Exponentiate diagonal terms [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(Hcopy[i + N*i] * f);
        im_U[i] = -sin(Hcopy[i + N*i] * f);
    }

    // Compute temperature adjustments
    thermal_correction(non, Hcopy, e, abs);

    // Exponentiate the couplings
    diagonalizeLPD(Hcopy, ecopy, N);
    // Multiply with adjusted couplings.
    matrix_on_vector(Hcopy, cr, ci, N);
    for (i = 0; i < N; i++) {
        Urcopy[i] = cos(ecopy[i] * f);
        Uicopy[i] = -sin(ecopy[i] * f);
    }
    vector_on_vector(Urcopy, Uicopy, cr, ci, N);
    trans_matrix_on_vector(Hcopy, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Adiabatic")) {
        // Convert average -> local basis
        matrix_on_vector(H_avg, cr, ci, N);
    }

    free(abs);
    free(Hcopy);
    free(ecopy);
    free(Urcopy);
    free(Uicopy);
}

// Expects wavefunction in adiabatic basis
void propagate_nise_dbc(
    t_non *non, 
    float *H_old, 
    float *H_new,
    float *e_old,
    float *e_new, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci
) {
    float f;
    int index, N, N2;
    float *ecopy;
    float *Hcopy;
    float *Hcc;
    float *Urcopy;
    float *Uicopy;
    float qc_ij, qc_ji, ediff;
    float *abs;
    int i, j, k;
    
    N = non->singles;
    N2 = N*N;
    f = non->deltat * icm2ifs * twoPi;

    abs = (float *) calloc(N, sizeof(float));
    Hcopy = (float *) calloc(N2, sizeof(float));
    Hcc = (float *) calloc(N2, sizeof(float));
    ecopy = (float *) calloc(N, sizeof(float));
    Urcopy = (float *) calloc(N, sizeof(float));
    Uicopy = (float *) calloc(N, sizeof(float));

    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Average")) {
        matrix_on_vector(H_old, cr, ci, N);
    }

    // Compute absolute values of the wavefunction coeffs
    for (i = 0; i < N; i++) {
        abs[i] = sqrt(cr[i]*cr[i] + ci[i]*ci[i]);
    }
    
    clearvec(Hcc, N2);
    // Find site basis perturbation
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            // Diagonals
            Hcc[i + N*i] += e_new[j] * H_new[i + N*j] * H_new[i + N*j];
            // Off-diagonals
            for (k = 0; k < N; k++) {
                Hcc[i + N*j] += e_new[k] * H_new[i + N*k] * H_new[j + N*k];
            }
        }
    }

    // Find average eigenbasis perturbation
    copyvec(H_old, Hcopy, N2);
    matrix_on_matrix(H_old, Hcc, N);
    matrix_on_matrix(Hcc, Hcopy, N);
    free(Hcc);

    // Exponentiate diagonal terms [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(Hcopy[i + N*i] * f);
        im_U[i] = -sin(Hcopy[i + N*i] * f);
    }

    // Compute temperature adjustments
    thermal_correction(non, Hcopy, e_new, abs);
    // Exponentiate the couplings
    diagonalizeLPD(Hcopy, ecopy, N);
    // Multiply with adjusted couplings.
    matrix_on_vector(Hcopy, cr, ci, N);
    for (i = 0; i < N; i++) {
        Urcopy[i] = cos(ecopy[i] * f);
        Uicopy[i] = -sin(ecopy[i] * f);
    }
    vector_on_vector(Urcopy, Uicopy, cr, ci, N);
    trans_matrix_on_vector(Hcopy, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Average")) {
        trans_matrix_on_vector(H_new, cr, ci, N);
    }

    free(abs);
    free(Hcopy);
    free(ecopy);
    free(Urcopy);
    free(Uicopy);
}

// Expects wavefunction in adiabatic basis
void propagate_tnise(
    t_non *non, 
    float *H_old, 
    float *H_new, 
    float *e, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci
) {
    float f;
    int index, N;
    float qc_ij, qc_ji, ediff;
    int i, j, k;
    float norm2;
    
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;
    float kBT=non->temperature*0.695; // Kelvin to cm-1

    // Exponentiate [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(e[i] * f);
        im_U[i] = -sin(e[i] * f);
    }

    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Average")) {
        matrix_on_vector(H_old, cr, ci, N);
    }

    // Compute unadjusted non-adiabatic coupling
    matrix_on_matrix(H_new, H_old, N);
    // Compute temperature adjustments
    for (i = 0; i < N; i++) {
        norm2 = 0;
        for (j = 0; j < N; j++) {
            if (j != i) {
                ediff = e[i] - e[j];
                H_old[i + N*j] *= exp(-0.25 * ediff / kBT);
                norm2 += H_old[i + N*j] * H_old[i + N*j];
            }
        }
        // Renormalise
        H_old[i + N*i] = sqrt(1 - norm2);
    }

    // Multiply with (real) non-adiabatic propagator
    trans_matrix_on_vector(H_old, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    if (!strcmp(non->basis, "Local") || !strcmp(non->basis, "Average")) {
        trans_matrix_on_vector(H_new, cr, ci, N);
    }

    norm2 = 0;
    for (i = 0; i < N; i++) {
        norm2 += cr[i]*cr[i] + ci[i]*ci[i];
    }
    for (i = 0; i < N; i++) {
        cr[i] /= sqrt(norm2);
        ci[i] /= sqrt(norm2);
    }
}

void thermal_correction(
    t_non* non, 
    float *H_old, 
    float *e, 
    float *abs
) {
    int N;
    float qc_ij, qc_ji, ediff;
    int i, j;
    float boltz;
    
    N = non->singles;
    float kBT=non->temperature*0.695; // Kelvin to cm-1

    for (i = 0; i < N; i++) {
        // Reset diagonal
        H_old[i + N*i] = 0;
        for (j = 0; j < i; j++) {
            ediff = e[i] - e[j];
            boltz = exp(ediff / kBT);
            qc_ij = sqrt(2 / (1 + boltz));
            qc_ji = qc_ij * sqrt(boltz);
            // Correction by Kleinekathoefer
            if (abs[i] - abs[j] != 0) {
                H_old[i + N*j] *= abs[j]*qc_ij - abs[i]*qc_ji;
                H_old[i + N*j] /= abs[j] - abs[i];
            } else {
                H_old[i + N*j] *= (qc_ij - qc_ji) / 2;
            }
            H_old[j + N*i] = -H_old[i + N*j];
        }
    }
}