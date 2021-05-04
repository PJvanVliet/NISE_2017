#include <stdio.h>
#include <stdlib.h>
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
#include "pop_single_t2.h"
#include <stdarg.h>

// Print population results to file
void pop_print(char* filename, float* pop_nise, float* pop_prezhdo, float* pop_alt, t_non* non, int sampleCount) {
    FILE* out = fopen(filename, "w");
    for (int t2 = 0; t2 < non->tmax2 + 1; t2++) {
            pop_nise[t2] /= sampleCount;
            pop_prezhdo[t2] /= sampleCount;
            pop_alt[t2] /= sampleCount;
            fprintf(out, "%f %e %e %e\n", t2 * non->deltat, pop_nise[t2], pop_prezhdo[t2], pop_alt[t2]);
    }
    fclose(out);
}

void reset_wavefn(float* cr_nise, float* ci_nise, float* cr_prezhdo, float* ci_prezhdo, float* cr_alt, float* ci_alt, int N) {
    clearvec(cr_nise, N);
    clearvec(ci_nise, N);
    clearvec(cr_prezhdo, N);
    clearvec(ci_prezhdo, N);
    clearvec(cr_alt, N);
    clearvec(ci_alt, N);
    cr_nise[0] = cr_prezhdo[0] = cr_alt[0] = 1.0;
    ci_nise[0] = ci_prezhdo[0] = ci_alt[0] = 0.0;
}

void propagate_NISE(
    t_non *non,
    float *H, 
    float *e, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci,
    float *pop_nise,
    float *pop_nise_ad,
    int t2
) {
    float f;
    int index, N;
    float re, im;
    int i;

    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;

    // Exponentiate [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(e[i] * f);
        im_U[i] = -sin(e[i] * f);
    }

    // Transfer to eigen basis
    matrix_on_vector(H,cr,ci,N);

    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    float pop = cr[0]*cr[0] + ci[0]*ci[0];
    pop_nise_ad[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
    
    // Transfer back to site basis
    trans_matrix_on_vector(H,cr,ci,N);

    pop = cr[0]*cr[0] + ci[0]*ci[0];
    pop_nise[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
}

void propagate_prezhdo(
    t_non *non, 
    float *H_old, 
    float *H_new, 
    float *e, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci,
    float *pop_prezhdo,
    float *pop_prezhdo_ad,
    int t2
) {
    float f;
    int index, N;
    float qc_ij, qc_ji, ediff;
    float re, im;
    float* abs;
    int i, j, k;
    
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;
    float kBT=non->temperature*0.695; // Kelvin to cm-1

    abs = (float *)calloc(N, sizeof(float));

    // Compute absolute values of the wavefunction coeffs
    for (i = 0; i < N; i++) {
        abs[i] = sqrt(cr[i]*cr[i] + ci[i]*ci[i]);
    }

    // Exponentiate [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(e[i] * f);
        im_U[i] = -sin(e[i] * f);
    }

    // Convert to adiabatic basis
    matrix_on_vector(H_old, cr, ci, N);

    // Compute unadjusted non-adiabatic coupling
    matrix_on_matrix(H_new, H_old, N);

    // Compute temperature adjustments
    for (i = 0; i < N; i++) {
        // Reset diagonal
        H_old[i + N*i] = 0;
        for (j = 0; j < i; j++) {
            ediff = e[i] - e[j];
            qc_ij = sqrt(2 / (1 + exp(ediff / kBT)));
            qc_ji = sqrt(2 / (1 + exp(-ediff / kBT)));
            H_old[i + N*j] *= abs[j]*qc_ij + abs[i]*qc_ji;
            // Correction by Kleinekathofer
            H_old[i + N*j] /= abs[i] + abs[j];
            H_old[j + N*i] = -H_old[i + N*j];
        }
    }

    // Exponentiate the non-adiabatic couplings
    matrix_exp(H_old, N);

    // Multiply with (real) non-adiabatic propagator
    trans_matrix_on_vector(H_old, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    float pop = cr[0]*cr[0] + ci[0]*ci[0];
    pop_prezhdo_ad[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
    
    // Transform back to site basis
    trans_matrix_on_vector(H_new, cr, ci, N);

    pop = cr[0]*cr[0] + ci[0]*ci[0];
    pop_prezhdo[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];

    free(abs);
}

void propagate_alt(
    t_non *non, 
    float *H_old, 
    float *H_new, 
    float *e, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci,
    float *pop_alt,
    float *pop_alt_ad,
    int t2
) {
    float f;
    int index, N;
    float qc_ij, qc_ji, ediff;
    float re, im;
    float* abs;
    int i, j, k;
    
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;
    float kBT=non->temperature*0.695; // Kelvin to cm-1

    abs = (float *)calloc(N, sizeof(float));

    // Compute absolute values of the wavefunction coeffs
    for (i = 0; i < N; i++) {
        abs[i] = sqrt(cr[i]*cr[i] + ci[i]*ci[i]);
    }

    // Exponentiate [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos(e[i] * f);
        im_U[i] = -sin(e[i] * f);
    }
    
    // Convert to adiabatic basis
    matrix_on_vector(H_old, cr, ci, N);

    // Compute unadjusted non-adiabatic coupling
    matrix_on_matrix(H_new, H_old, N);

    // Compute temperature adjustments
    for (i = 0; i < N; i++) {
        // Reset diagonals
        H_old[i + N*i] = 0;
        for (j = 0; j < i; j++) {
            ediff = e[i] - e[j];
            // TODO: make sure there are no overflows!
            if (ediff > 0) {
                qc_ij = exp(-ediff / kBT);
                qc_ji = 1;
            } else {
                qc_ij = 1;
                qc_ji = exp(ediff / kBT);
            }
            H_old[i + N*j] *= abs[j]*qc_ij + abs[i]*qc_ji;
            // Correction by Kleinekathofer
            H_old[i + N*j] /= abs[i] + abs[j];
            H_old[j + N*i] = -H_old[i + N*j];
        }
    }

    // Exponentiate the non-adiabatic couplings
    matrix_exp(H_old, N);

    // Multiply with (real) non-adiabatic propagator
    trans_matrix_on_vector(H_old, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    float pop = cr[0]*cr[0] + ci[0]*ci[0];
    pop_alt_ad[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
    
    // Transform back to site basis
    trans_matrix_on_vector(H_new, cr, ci, N);

    pop = cr[0]*cr[0] + ci[0]*ci[0];
    pop_alt[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];

    free(abs);
}

void row_swap(float* a, int row1, int row2, int N) {
    for (int i = 0; i < N; i++) {
        float temp = a[i + row1*N];
        a[i + row1*N] = a[i + row2*N];
        a[i + row2*N] = temp; 
    }
}

void swaps(float* H_new, float* H_old, int N) {
    int N2;
    float *Hcopy;
    float temp;
    int imax, jmax;
    int i, j, k;
    float min_diag, max_offdiag;
    
    N2 = N * N;
    
    Hcopy = (float *)calloc(N2, sizeof(float));
    copyvec(H_old, Hcopy, N2);
    matrix_on_matrix(H_new, Hcopy, N);

    // Determine maximum off-diagonal 
    // and minimum diagonal elements.
    imax = 0, jmax = 0;
    min_diag = 1, max_offdiag = 0;
    for (i = 0; i < N; i++) {
        temp = fabs(Hcopy[i + N*i]);
        if (temp < min_diag) {
            min_diag = temp;
        }
        for (j = 0; j < i; j++) {
            temp = fabs(Hcopy[j + N*i]);
            if (temp > max_offdiag) {
                imax = i, jmax = j;
                max_offdiag = temp;
            }
        }
    }

    // Check if rows must be swapped
    while (min_diag < max_offdiag) {
        row_swap(H_new, imax, jmax, N);

        // Recompute min_diag and max_offdiag
        copyvec(H_old, Hcopy, N2);
        matrix_on_matrix(H_new, Hcopy, N);

        // Determine maximum off-diagonal 
        // and minimum diagonal elements.
        imax = 0, jmax = 0;
        min_diag = 1, max_offdiag = 0;
        for (i = 0; i < N; i++) {
            temp = fabs(Hcopy[i + N*i]);
            if (temp < min_diag) {
                min_diag = temp;
            }
            for (j = 0; j < i; j++) {
                temp = fabs(Hcopy[j + N*i]);
                if (temp > max_offdiag) {
                    imax = i, jmax = j;
                    max_offdiag = temp;
                }
            }
        }
    }

    // Check whether we need to change the sign of the eigenvectors
    copyvec(H_old, Hcopy, N2);
    matrix_on_matrix(H_new, Hcopy, N);
    for (i = 0; i < N; i++) {
        if (Hcopy[i + N*i] < 0) {
            for (j = 0; j < N; j++) {
                H_new[i + N*j] = -H_new[i + N*j];
            }
        }
    }
    
    free(Hcopy);
}

void pop_single_t2(t_non* non) {
    // Initialise base variables
    int N = non->singles;
    int nn2 = N * (N + 1) / 2;
    int N2 = N * N;
    int ti, tm, samples;
    int sampleCount;
    int i, j, k;
    float *Hamil_i_e;
    float *H_new;
    float *H_old;
    float *Hcopy;
    float *re_U;
    float *im_U;
    float *e;
    float *cr_nise;
    float *ci_nise;
    float *cr_prezhdo;
    float *ci_prezhdo;
    float *cr_alt;
    float *ci_alt;
    float *pop_nise;
    float *pop_prezhdo;
    float *pop_alt;
    float *pop_nise_ad;
    float *pop_prezhdo_ad;
    float *pop_alt_ad;

    Hamil_i_e = (float *) calloc(nn2, sizeof(float));
    H_new = (float *) calloc(N2, sizeof(float));
    H_old = (float *) calloc(N2, sizeof(float));
    Hcopy = (float *)calloc(N2, sizeof(float));
    re_U = (float *) calloc(N, sizeof(float));
    im_U = (float *) calloc(N, sizeof(float));
    e = (float *) calloc(N, sizeof(float));
    cr_nise = (float *) calloc(N, sizeof(float));
    ci_nise = (float *) calloc(N, sizeof(float));
    cr_prezhdo = (float *) calloc(N, sizeof(float));
    ci_prezhdo = (float *) calloc(N, sizeof(float));
    cr_alt = (float *) calloc(N, sizeof(float));
    ci_alt = (float *) calloc(N, sizeof(float));
    pop_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_prezhdo = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_alt = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_nise_ad = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_prezhdo_ad = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_alt_ad = (float *) calloc(non->tmax2 + 1, sizeof(float));
    
    // Determine number of samples
    sampleCount = (non->length - non->tmax2 - 1) / non->sample + 1;
    printf("Total number of samples: %i\n", sampleCount);
    // Set the initial population
    pop_nise[0] = pop_prezhdo[0] = pop_alt[0] = 1.0 * sampleCount;

    // Read the Hamiltonian file
    FILE* H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        exit(1);
    }

    // Loop over samples
    for (samples = 0; samples < sampleCount; samples++) {
        ti = samples * non->sample;

        // Load first Hamiltonian
        if (read_He(non, Hamil_i_e, H_traj, ti) != 1) {
            printf("Hamiltonian trajectory file too short, could not fill buffer!\n");
            exit(1);
        }
        build_diag_H(Hamil_i_e, H_new, e, N);

        // Reset the wavefunctions
        reset_wavefn(cr_nise, ci_nise, cr_prezhdo, ci_prezhdo, cr_alt, ci_alt, N);

        // Start integrating the Schrödinger equation
        // NISE:
        for (int t2 = 0; t2 < non->tmax2; t2++) {
            int tm = ti + t2 + 1;
            // Copy old Hamiltonian
            copyvec(H_new, H_old, N2);
            // Load new Hamiltonian
            if (read_He(non, Hamil_i_e, H_traj, tm) != 1) {
                printf("Hamiltonian trajectory file too short, could not fill buffer!!!\n");
                exit(1);
            }
            build_diag_H(Hamil_i_e, H_new, e, N);

            // Check if we need to perform any swaps, to maximise overlap
            swaps(H_new, H_old, N);

            // Run NISE
            propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise, pop_nise, pop_nise_ad, t2);


            // Run Prezhdo
            copyvec(H_old, Hcopy, N2);
            propagate_prezhdo(non, Hcopy, H_new, e, re_U, im_U, cr_prezhdo, ci_prezhdo, pop_prezhdo, pop_prezhdo_ad, t2);
            pop = cr_prezhdo[0]*cr_prezhdo[0] + ci_prezhdo[0]*ci_prezhdo[0];

            pop_prezhdo[t2 + 1] += pop;
            // Run alt
            propagate_alt(non, H_old, H_new, e, re_U, im_U, cr_alt, ci_alt, pop_alt, pop_alt_ad, t2);
            pop = cr_alt[0] * cr_alt[0] + ci_alt[0] * ci_alt[0];

            pop_alt[t2 + 1] += pop;
        }
    }
    printf("\n");
    
    char* fn = "pop_t2.txt";
    pop_print(fn, pop_nise, pop_prezhdo, pop_alt, non, sampleCount);

    free(Hamil_i_e), free(H_new), free(H_old), free(Hcopy);
    free(re_U), free(im_U), free(e);
    free(cr_nise), free(ci_nise);
    free(cr_prezhdo), free(ci_prezhdo);
    free(cr_alt), free(ci_alt);
    free(pop_nise), free(pop_prezhdo), free(pop_alt);
}