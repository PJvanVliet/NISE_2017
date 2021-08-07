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
void pop_print(char* filename, float* pop_nise, float* pop_nise_dba, float* pop_nise_dbb, t_non* non, int sampleCount) {
    FILE* out = fopen(filename, "w");
    for (int t2 = 0; t2 < non->tmax2 + 1; t2++) {
            pop_nise[t2] /= sampleCount;
            pop_nise_dba[t2] /= sampleCount;
            pop_nise_dbb[t2] /= sampleCount;
            fprintf(out, "%f %e %e %e\n", t2 * non->deltat, pop_nise[t2], pop_nise_dba[t2], pop_nise_dbb[t2]);
    }
    fclose(out);
}

void coh_print(char* filename, float* cohr_nise, float* cohi_nise, float* cohr_nise_dba, float* cohi_nise_dba, float* cohr_nise_dbb, float* cohi_nise_dbb, t_non* non, int sampleCount) {
    FILE* out = fopen(filename, "w");
    for (int t2 = 0; t2 < non->tmax2 + 1; t2++) {
            cohr_nise[t2] /= sampleCount;
            cohr_nise_dba[t2] /= sampleCount;
            cohr_nise_dbb[t2] /= sampleCount;
            cohi_nise[t2] /= sampleCount;
            cohi_nise_dba[t2] /= sampleCount;
            cohi_nise_dbb[t2] /= sampleCount;
            cohr_nise[t2] = sqrt(cohr_nise[t2]*cohr_nise[t2] + cohi_nise[t2]*cohi_nise[t2]);
            cohr_nise_dba[t2] = sqrt(cohr_nise_dba[t2]*cohr_nise_dba[t2] + cohi_nise_dba[t2]*cohi_nise_dba[t2]);
            cohr_nise_dbb[t2] = sqrt(cohr_nise_dbb[t2]*cohr_nise_dbb[t2] + cohi_nise_dbb[t2]*cohi_nise_dbb[t2]);
            fprintf(out, "%f %e %e %e\n", t2 * non->deltat, cohr_nise[t2], cohr_nise_dba[t2], cohr_nise_dbb[t2]);
    }
    fclose(out);
}

void reset_wavefn(float* cr_nise, float* ci_nise, float* cr_nise_dba, float* ci_nise_dba, float* cr_nise_dbb, float* ci_nise_dbb, int N) {
    clearvec(cr_nise, N);
    clearvec(ci_nise, N);
    clearvec(cr_nise_dba, N);
    clearvec(ci_nise_dba, N);
    clearvec(cr_nise_dbb, N);
    clearvec(ci_nise_dbb, N);
    cr_nise[1] = cr_nise_dba[1] = cr_nise_dbb[1] = 1.0;
}

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
    float re, im;
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

    // Compute unadjusted non-adiabatic coupling
    matrix_on_matrix(H_new, H_old, N);

    // Compute temperature adjustments
    for (i = 0; i < N; i++) {
        // Reset diagonal
        H_old[i + N*i] = 0;
        for (j = 0; j < i; j++) {
            ediff = e[i] - e[j];
            float boltz = exp(ediff / kBT);
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

    free(abs);
}

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
    float *Urcopy;
    float *Uicopy;
    float qc_ij, qc_ji, ediff;
    float re, im;
    float* abs;
    int i, j, k;
    
    N = non->singles;
    N2 = N*N;
    f = non->deltat * icm2ifs * twoPi;
    float kBT=non->temperature*0.695; // Kelvin to cm-1

    abs = (float *)calloc(N, sizeof(float));
    Hcopy = (float *) calloc(N2, sizeof(float));
    ecopy = (float *) calloc(N, sizeof(float));
    Urcopy = (float *) calloc(N, sizeof(float));
    Uicopy = (float *) calloc(N, sizeof(float));
    clearvec(Hcopy, N2);
    clearvec(Urcopy, N);
    clearvec(Uicopy, N);

    // Compute absolute values of the wavefunction coeffs
    for (i = 0; i < N; i++) {
        abs[i] = sqrt(cr[i]*cr[i] + ci[i]*ci[i]);
    }

    // Find back site basis perturbation (only diagonals!)
    for (i = 0; i < N; i++) {
        for (k = 0; k < N; k++) {
            ecopy[i] += e[k] * H_new[i + N*k] * H_new[i + N*k];
            ecopy[i] -= e_avg[k] * H_avg[i + N*k] * H_avg[i + N*k];
        }
    }

    // Find back average eigenbasis perturbation
    for (i = 0; i < N; i++) {
        // Diagonals
        for (k = 0; k < N; k++) {
            Hcopy[i + N*i] += ecopy[k] * H_avg[i + N*k] * H_avg[i + N*k];
        }
        for (j = 0; j < i; j++) {
            for (k = 0; k < N; k++) {
                Hcopy[i + N*j] += ecopy[k] * H_avg[i + N*k] * H_avg[i + N*j];
            }
            Hcopy[j + N*i] = Hcopy[i + N*j];
        }
    }

    // Exponentiate average Hamiltonian [U=exp(-i/h H dt)]
    for (i = 0; i < N; i++) {
        re_U[i] = cos((e_avg[i] + Hcopy[i + N*i]) * f);
        im_U[i] = -sin((e_avg[i] + Hcopy[i + N*i]) * f);
    }

    // Compute temperature adjustments
    for (i = 0; i < N; i++) {
        // Reset diagonal
        Hcopy[i + N*i] = 0;
        for (j = 0; j < i; j++) {
            ediff = e[i] - e[j];
            float boltz = exp(ediff / kBT);
            qc_ij = sqrt(2 / (1 + boltz));
            qc_ji = qc_ij * sqrt(boltz);
            // Correction by Kleinekathofer
            if (abs[i] - abs[j] != 0) {
                Hcopy[i + N*j] *= abs[j]*qc_ij - abs[i]*qc_ji;
                Hcopy[i + N*j] /= abs[j] - abs[i];
            } else {
                Hcopy[i + N*j] *= (qc_ij - qc_ji) / 2;
            }
            Hcopy[j + N*i] = Hcopy[i + N*j];
        }
    }

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
    float *H_avg;
    float *e_avg;
    float *H_new;
    float *H_old;
    float *Hcopy;
    float *re_U;
    float *im_U;
    float *e;
    float *cr_nise;
    float *ci_nise;
    float *cr_nise_dba;
    float *ci_nise_dba;
    float *cr_nise_dbb;
    float *ci_nise_dbb;
    float *pop_nise;
    float *pop_nise_dba;
    float *pop_nise_dbb;
    float *cohr_nise;
    float *cohr_nise_dba;
    float *cohr_nise_dbb;
    float *cohi_nise;
    float *cohi_nise_dba;
    float *cohi_nise_dbb;

    Hamil_i_e = (float *) calloc(nn2, sizeof(float));
    H_avg = (float *) calloc(N2, sizeof(float));
    e_avg = (float *) calloc(N, sizeof(float));
    H_new = (float *) calloc(N2, sizeof(float));
    H_old = (float *) calloc(N2, sizeof(float));
    Hcopy = (float *)calloc(N2, sizeof(float));
    re_U = (float *) calloc(N, sizeof(float));
    im_U = (float *) calloc(N, sizeof(float));
    e = (float *) calloc(N, sizeof(float));
    cr_nise = (float *) calloc(N, sizeof(float));
    ci_nise = (float *) calloc(N, sizeof(float));
    cr_nise_dba = (float *) calloc(N, sizeof(float));
    ci_nise_dba = (float *) calloc(N, sizeof(float));
    cr_nise_dbb = (float *) calloc(N, sizeof(float));
    ci_nise_dbb = (float *) calloc(N, sizeof(float));
    pop_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_nise_dba = (float *) calloc(non->tmax2 + 1, sizeof(float));
    pop_nise_dbb = (float *) calloc(non->tmax2 + 1, sizeof(float));
    cohr_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
    cohr_nise_dba = (float *) calloc(non->tmax2 + 1, sizeof(float));
    cohr_nise_dbb = (float *) calloc(non->tmax2 + 1, sizeof(float));
    cohi_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
    cohi_nise_dba = (float *) calloc(non->tmax2 + 1, sizeof(float));
    cohi_nise_dbb = (float *) calloc(non->tmax2 + 1, sizeof(float));
    
    // Determine number of samples
    sampleCount = (non->length - non->tmax2 - 1) / non->sample + 1;
    printf("Total number of samples: %i\n", sampleCount);
    // Set the initial population
    pop_nise[0] = pop_nise_dba[0] = pop_nise_dbb[0] = 1.0 * sampleCount;

    // Read the Hamiltonian file
    FILE* H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        exit(1);
    }

    // Find average Hamiltonian
    clearvec(H_avg, N2);
    for (int t2 = 0; t2 < (non->length - non->tmax2 - 1); t2++) {
        if (read_He(non, Hamil_i_e, H_traj, t2) != 1) {
            printf("Hamiltonian trajectory file too short, could not fill buffer!\n");
            exit(1);
        }
        int a, b;
        // Build square Hamiltonian from triagonal matrix
        for (a = 0; a < N; a++) {
            H_avg[a + N * a] += Hamil_i_e[a + N * a - (a * (a + 1)) / 2]; // Diagonal
            for (b = a + 1; b < N; b++) {
                H_avg[a + N * b] += Hamil_i_e[b + N * a - (a * (a + 1)) / 2];
                H_avg[b + N * a] += Hamil_i_e[b + N * a - (a * (a + 1)) / 2];
            }
        }
    }

    for (a = 0; a < N; a++) {
        H_avg[a + N * a] /= (non->length - non->tmax2); // Diagonal
        for (b = a + 1; b < N; b++) {
            H_avg[a + N * b] /= (non->length - non->tmax2);
            H_avg[b + N * a] /= (non->length - non->tmax2);
        }
    }

    // Diagonalise Hamiltonian
    diagonalizeLPD(H_avg, e_avg, N);

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
        reset_wavefn(cr_nise, ci_nise, cr_nise_dba, ci_nise_dba, cr_nise_dbb, ci_nise_dbb, N);

        // Start NISE procedure
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

            if (!strcmp(non->basis, "Local")) {
                // NISE
                propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise);
                pop_nise[t2 + 1] += cr_nise[1] * cr_nise[1] + ci_nise[1] * ci_nise[1];
                cohr_nise[t2 + 1] += cr_nise[0] * cr_nise[1] + ci_nise[0] * ci_nise[1];
                cohi_nise[t2 + 1] += ci_nise[0] * ci_nise[1] - ci_nise[1] * cr_nise[0];
                
                // NISE-DBa
                // Transfer to adiabatic basis
                matrix_on_vector(H_old, cr_nise_dba, ci_nise_dba, N);
                // Propagate
                copyvec(H_old, Hcopy, N2);
                propagate_nise_dba(non, Hcopy, H_new, e, re_U, im_U, cr_nise_dba, ci_nise_dba);
                // Transfer back to site basis
                trans_matrix_on_vector(H_new, cr_nise_dba, ci_nise_dba, N);
                pop_nise_dba[t2 + 1] += cr_nise_dba[1]*cr_nise_dba[1] + ci_nise_dba[1]*ci_nise_dba[1];
                cohr_nise_dba[t2 + 1] += cr_nise_dba[0] * cr_nise_dba[1] + ci_nise_dba[0] * ci_nise_dba[1];
                cohi_nise_dba[t2 + 1] += ci_nise_dba[0] * cr_nise_dba[1] - ci_nise_dba[1] * cr_nise_dba[0];

                // NISE-DBb
                // Transform to average eigenbasis
                trans_matrix_on_vector(H_avg, cr_nise_dbb, ci_nise_dbb, N);
                // Propagate
                propagate_nise_dbb(non, H_avg, H_new, e_avg, e, re_U, im_U, cr_nise_dbb, ci_nise_dbb);
                // Transform back to site basis
                matrix_on_vector(H_avg, cr_nise_dbb, ci_nise_dbb, N);
                pop_nise_dbb[t2 + 1] += cr_nise_dbb[1]*cr_nise_dbb[1] + ci_nise_dbb[1]*ci_nise_dbb[1];
                cohr_nise_dbb[t2 + 1] += cr_nise_dbb[0] * cr_nise_dbb[1] + ci_nise_dbb[0] * ci_nise_dbb[1];
                cohi_nise_dbb[t2 + 1] += ci_nise_dbb[0] * cr_nise_dbb[1] - ci_nise_dbb[1] * cr_nise_dbb[0];

            } else if (!strcmp(non->basis, "Adiabatic")) {
                // NISE
                // Transfer to site
                matrix_on_vector(H_old, cr_nise, ci_nise, N);
                // Propagate
                propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise);
                // Transfer back to adiabatic basis
                trans_matrix_on_vector(H_new, cr_nise, ci_nise, N);
                pop_nise[t2 + 1] += cr_nise[1] * cr_nise[1] + ci_nise[1] * ci_nise[1];
                cohr_nise[t2 + 1] += cr_nise[0] * cr_nise[1] + ci_nise[0] * ci_nise[1];
                cohi_nise[t2 + 1] += ci_nise[0] * cr_nise[1] - ci_nise[1] * cr_nise[0];

                // NISE-DBa
                copyvec(H_old, Hcopy, N2);
                propagate_nise_dba(non, Hcopy, H_new, e, re_U, im_U, cr_nise_dba, ci_nise_dba);
                pop_nise_dba[t2 + 1] += cr_nise_dba[1]*cr_nise_dba[1] + ci_nise_dba[1]*ci_nise_dba[1];
                cohr_nise_dba[t2 + 1] += cr_nise_dba[0] * cr_nise_dba[1] + ci_nise_dba[0] * ci_nise_dba[1];
                cohi_nise_dba[t2 + 1] += ci_nise_dba[0] * cr_nise_dba[1] - ci_nise_dba[1] * cr_nise_dba[0];

                // NISE-DBb
                // Transfer adiabatic -> site basis
                trans_matrix_on_vector(H_new, cr_nise_dbb, ci_nise_dbb, N);
                // Transfer site -> average eigenbasis
                trans_matrix_on_vector(H_avg, cr_nise_dbb, ci_nise_dbb, N);
                // Propagate
                propagate_nise_dbb(non, H_avg, H_new, e_avg, e, re_U, im_U, cr_nise_dbb, ci_nise_dbb);
                // Transfer average -> site basis
                matrix_on_vector(H_avg, cr_nise_dbb, ci_nise_dbb, N);
                // Transfer site -> adiabatic basis
                matrix_on_vector(H_new, cr_nise_dbb, ci_nise_dbb, N);
                pop_nise_dbb[t2 + 1] += cr_nise_dbb[1]*cr_nise_dbb[1] + ci_nise_dbb[1]*ci_nise_dbb[1];
                cohr_nise_dbb[t2 + 1] += cr_nise_dbb[0] * cr_nise_dbb[1] + ci_nise_dbb[0] * ci_nise_dbb[1];
                cohi_nise_dbb[t2 + 1] += ci_nise_dbb[0] * cr_nise_dbb[1] - ci_nise_dbb[1] * cr_nise_dbb[0];

            } else if (!strcmp(non->basis, "Average")) {
                // NISE
                // Transfer to site basis
                trans_matrix_on_vector(H_avg, cr_nise, ci_nise, N);
                // Propagate
                propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise);
                // Transfer back to average eigenbasis
                trans_matrix_on_vector(H_avg, cr_nise, ci_nise, N);
                pop_nise[t2 + 1] += cr_nise[1] * cr_nise[1] + ci_nise[1] * ci_nise[1];
                cohr_nise[t2 + 1] += cr_nise[0] * cr_nise[1] + ci_nise[0] * ci_nise[1];
                cohi_nise[t2 + 1] += ci_nise[0] * cr_nise[1] - ci_nise[1] * cr_nise[0];

                // NISE-DBa
                copyvec(H_old, Hcopy, N2);
                // Transfer average -> site
                trans_matrix_on_vector(H_avg, cr_nise_dba, ci_nise_dba, N);
                // Transfer site -> adiabatic
                trans_matrix_on_vector(H_new, cr_nise_dba, ci_nise_dba, N);
                // Propagate
                propagate_nise_dba(non, Hcopy, H_new, e, re_U, im_U, cr_nise_dba, ci_nise_dba);
                // Transfer adiabatic -> site
                matrix_on_vector(H_new, cr_nise_dba, ci_nise_dba, N);
                // Transfer site -> average
                matrix_on_vector(H_avg, cr_nise_dba, ci_nise_dba, N);
                pop_nise_dba[t2 + 1] += cr_nise_dba[1]*cr_nise_dba[1] + ci_nise_dba[1]*ci_nise_dba[1];
                cohr_nise_dba[t2 + 1] += cr_nise_dba[0] * cr_nise_dba[1] + ci_nise_dba[0] * ci_nise_dba[1];
                cohi_nise_dba[t2 + 1] += ci_nise_dba[0] * cr_nise_dba[1] - ci_nise_dba[1] * cr_nise_dba[0];

                // NISE-DBb
                copyvec(H_old, Hcopy, N2);
                propagate_nise_dbb(non, H_avg, H_new, e_avg, e, re_U, im_U, cr_nise_dbb, ci_nise_dbb);
                pop_nise_dbb[t2 + 1] += cr_nise_dbb[1]*cr_nise_dbb[1] + ci_nise_dbb[1]*ci_nise_dbb[1];
                cohr_nise_dbb[t2 + 1] += cr_nise_dbb[0] * cr_nise_dbb[1] + ci_nise_dbb[0] * ci_nise_dbb[1];
                cohi_nise_dbb[t2 + 1] += ci_nise_dbb[0] * cr_nise_dbb[1] - ci_nise_dbb[1] * cr_nise_dbb[0];
            }
        }
    }
    
    char* fn = "pop_t2.txt";
    pop_print(fn, pop_nise, pop_nise_dba, pop_nise_dbb, non, sampleCount);
    fn = "coh_t2.txt";
    coh_print(fn, cohr_nise, cohi_nise, cohr_nise_dba, cohi_nise_dba, cohr_nise_dbb, cohi_nise_dbb, non, sampleCount);

    free(Hamil_i_e);
    free(H_avg);
    free(e_avg);
    free(H_new);
    free(H_old);
    free(Hcopy);
    free(re_U);
    free(im_U);
    free(e);
    free(cr_nise);
    free(ci_nise);
    free(cr_nise_dba);
    free(ci_nise_dba);
    free(cr_nise_dbb);
    free(ci_nise_dbb);
    free(pop_nise);
    free(pop_nise_dba);
    free(pop_nise_dbb);
    free(cohr_nise);
    free(cohr_nise_dba);
    free(cohr_nise_dbb);
    free(cohi_nise);
    free(cohi_nise_dba);
    free(cohi_nise_dbb);
}