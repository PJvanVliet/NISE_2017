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
void pop_print(char* filename, float* pop_nise, float* pop_standard, float* pop_harmonic, t_non* non, int sampleCount) {
    FILE* out = fopen(filename, "w");
    for (int t2 = 0; t2 < non->tmax2 + 1; t2++) {
            pop_nise[t2] /= sampleCount;
            pop_standard[t2] /= sampleCount;
            pop_harmonic[t2] /= sampleCount;
            fprintf(out, "%f %e %e %e\n", t2 * non->deltat, pop_nise[t2], pop_standard[t2], pop_harmonic[t2]);
    }
    fclose(out);
}

void reset_wavefn(float* cr_nise, float* ci_nise, float* cr_standard, float* ci_standard, float* cr_harmonic, float* ci_harmonic, int N) {
    clearvec(cr_nise, N);
    clearvec(ci_nise, N);
    clearvec(cr_standard, N);
    clearvec(ci_standard, N);
    clearvec(cr_harmonic, N);
    clearvec(ci_harmonic, N);
    cr_nise[0] = cr_standard[0] = cr_harmonic[0] = 1.0;
    ci_nise[0] = ci_standard[0] = ci_harmonic[0] = 0.0;
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

    pop_nise_ad[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
    
    // Transfer back to site basis
    trans_matrix_on_vector(H,cr,ci,N);

    pop_nise[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
}

void propagate_standard(
    t_non *non, 
    float *H_old, 
    float *H_new, 
    float *e, 
    float *re_U, 
    float *im_U, 
    float *cr, 
    float *ci,
    float *pop_standard,
    float *pop_standard_ad,
    int t2
) {
    float f;
    int index, N;
    float qc_ij, qc_ji, ediff;
    float re, im;
    float* abs;
    int i, j, k;
    float* cpf;
    
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;
    float kBT=non->temperature*0.695; // Kelvin to cm-1

    abs = (float *)calloc(N, sizeof(float));
    cpf = (float *)calloc(N*N, sizeof(float));
    clearvec(cpf, N*N);

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
            float boltz = exp(-ediff / kBT);
            qc_ij = sqrt(2 / (1 + boltz));
            qc_ji = qc_ij * sqrt(boltz);
            // Correction by Kleinekathofer
            if (abs[i] - abs[j] != 0) {
                H_old[i + N*j] *= abs[j]*qc_ij - abs[i]*qc_ji;
                H_old[i + N*j] /= abs[j] - abs[i];
            } else {
                H_old[i + N*j] *= (qc_ij - qc_ji) / 2;
            }
            H_old[j + N*i] = -H_old[i + N*j];
        }
    }

    // Exponentiate the non-adiabatic couplings
    matrix_exp(H_old, N);

    // Multiply with (real) non-adiabatic propagator
    trans_matrix_on_vector(H_old, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);
    
    // // Compute coherence penalty functional (CPF)
    // for (i = 0; i < N; i++) {
    //     for (j = 0; j < i; j++) {
    //         float a = cr[i] * cr[j] + ci[i] * ci[j];
    //         float b = ci[i] * cr[j] - cr[i] * ci[j];
    //         cpf[i + N*j] = -non->deltat / non->dephasing * (a*a + b*b);
    //         cpf[j + N*i] = cpf[i + N*j];
    //     }
    // }
    // // Exponentiate the CPF
    // matrix_exp(cpf, N);
    // // Multiply with CPF
    // matrix_on_vector(cpf, cr, ci, N);

    pop_standard_ad[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
    
    // Transform back to site basis
    trans_matrix_on_vector(H_new, cr, ci, N);

    pop_standard[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];

    free(abs);
}

void propagate_harmonic(
    t_non *non, 
    float *H_old, 
    float *H_new, 
    float *e, 
    float *re_U, 
    float *im_U,
    float *cr, 
    float *ci,
    float *pop_harmonic,
    float *pop_harmonic_ad,
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
            float boltz = exp(-ediff / kBT);
            qc_ij = sqrt(ediff / kBT / (1 - boltz));
            qc_ji = qc_ij * sqrt(boltz);
            // Correction by Kleinekathofer
            if (abs[i] - abs[j] != 0) {
                H_old[i + N*j] *= abs[j]*qc_ij - abs[i]*qc_ji;
                H_old[i + N*j] /= abs[j] - abs[i];
            } else {
                H_old[i + N*j] *= (qc_ij - qc_ji) / 2;
            }
            H_old[j + N*i] = -H_old[i + N*j];
        }
    }
    
    // Exponentiate the non-adiabatic couplings
    matrix_exp(H_old, N);
    
    // Multiply with (real) non-adiabatic propagator
    trans_matrix_on_vector(H_old, cr, ci, N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);

    pop_harmonic_ad[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];
    
    // Transform back to site basis
    trans_matrix_on_vector(H_new, cr, ci, N);

    pop_harmonic[t2 + 1] += cr[0]*cr[0] + ci[0]*ci[0];

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
    float *cr_standard;
    float *ci_standard;
    float *cr_harmonic;
    float *ci_harmonic;
    float *pop_nise;
    float *pop_standard;
    float *pop_harmonic;
    float *pop_nise_ad;
    float *pop_standard_ad;
    float *pop_harmonic_ad;
    float *cr_copy;
    float *ci_copy;

    Hamil_i_e = (float *) calloc(nn2, sizeof(float));
    H_new = (float *) calloc(N2, sizeof(float));
    H_old = (float *) calloc(N2, sizeof(float));
    Hcopy = (float *)calloc(N2, sizeof(float));
    re_U = (float *) calloc(N, sizeof(float));
    im_U = (float *) calloc(N, sizeof(float));
    e = (float *) calloc(N, sizeof(float));
    cr_nise = (float *) calloc(N, sizeof(float));
    ci_nise = (float *) calloc(N, sizeof(float));
    cr_standard = (float *) calloc(N, sizeof(float));
    ci_standard = (float *) calloc(N, sizeof(float));
    cr_harmonic = (float *) calloc(N, sizeof(float));
    ci_harmonic = (float *) calloc(N, sizeof(float));
    pop_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_standard = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_harmonic = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_nise_ad = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_standard_ad = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_harmonic_ad = (float *) calloc(non->tmax2 + 1, sizeof(float));
    cr_copy = (float *) calloc(N, sizeof(float));
    ci_copy = (float *) calloc(N, sizeof(float));
    
    // Determine number of samples
    sampleCount = (non->length - non->tmax2 - 1) / non->sample + 1;
    printf("Total number of samples: %i\n", sampleCount);
    // Set the initial population
    pop_nise[0] = pop_standard[0] = pop_harmonic[0] = 1.0 * sampleCount;

    // Read the Hamiltonian file
    FILE* H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        exit(1);
    }

    // Loop over samples
    for (samples = 0; samples < sampleCount; samples++) {
        if (samples % 100 == 0) {
            printf("samples = %i\n", samples);
        }
        ti = samples * non->sample;

        // Load first Hamiltonian
        if (read_He(non, Hamil_i_e, H_traj, ti) != 1) {
            printf("Hamiltonian trajectory file too short, could not fill buffer!\n");
            exit(1);
        }
        build_diag_H(Hamil_i_e, H_new, e, N);

        // Reset the wavefunctions
        reset_wavefn(cr_nise, ci_nise, cr_standard, ci_standard, cr_harmonic, ci_harmonic, N);
        copyvec(cr_nise, cr_copy, N);
        copyvec(ci_nise, ci_copy, N);
        matrix_on_vector(H_new, cr_copy, ci_copy, N);
        float pop_copy = cr_copy[0] * cr_copy[0] + ci_copy[0] * ci_copy[0];
        pop_nise_ad[0] += pop_copy;
        pop_standard_ad[0] += pop_copy;
        pop_harmonic_ad[0] += pop_copy;

        // Start integrating the Schrödinger equation
        // NISE:
        for (int t2 = 0; t2 < non->tmax2; t2++) {
            int tm = ti + t2 + 1;
            // Copy old Hamiltonian
            copyvec(H_new, H_old, N2);
            // Load new Hamiltonian
            if (read_He(non, Hamil_i_e, H_traj, tm) != 1) {
                exit(1);
            }
            build_diag_H(Hamil_i_e, H_new, e, N);

            // Check if we need to perform any swaps, to maximise overlap
            swaps(H_new, H_old, N);

            // Run NISE
            propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise, pop_nise, pop_nise_ad, t2);

            // Run Prezhdo
            copyvec(H_old, Hcopy, N2);
            propagate_standard(non, Hcopy, H_new, e, re_U, im_U, cr_standard, ci_standard, pop_standard, pop_standard_ad, t2);

            // Run alt
            propagate_harmonic(non, H_old, H_new, e, re_U, im_U, cr_harmonic, ci_harmonic, pop_harmonic, pop_harmonic_ad, t2);
        }
    }
    
    char* fn = "pop_t2.txt";
    pop_print(fn, pop_nise, pop_standard, pop_harmonic, non, sampleCount);
    fn = "pop_t2_ad.txt";
    pop_print(fn, pop_nise_ad, pop_standard_ad, pop_harmonic_ad, non, sampleCount);

    free(Hamil_i_e), free(H_new), free(H_old), free(Hcopy);
    free(re_U), free(im_U), free(e);
    free(cr_nise), free(ci_nise);
    free(cr_standard), free(ci_standard);
    free(cr_harmonic), free(ci_harmonic);
    free(pop_nise), free(pop_standard), free(pop_harmonic);
    free(pop_nise_ad), free(pop_standard_ad), free(pop_harmonic_ad);
    free(cr_copy), free(ci_copy);
}