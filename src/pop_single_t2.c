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
#include "pop_single_t2.h"
#include "prop_schemes.h"


// Print population results to file
void pop_print(char* filename, t_non* non, int sampleCount, int count, ...) {
    int j;
    float *pop;
    FILE* out = fopen(filename, "w");
    
    va_list ap;
    for (int t2 = 0; t2 < non->tmax2 + 1; t2++) {
        fprintf(out, "%f ", t2 * non->deltat);
        va_start(ap, count);
        for (j = 0; j < count; j++) {
            pop = va_arg(ap, float*);
            fprintf(out, "%e ", pop[t2]/sampleCount);
        }
        fprintf(out, "\n");
        va_end(ap);
    }
    fclose(out);
}

// Print coherence results to file
void coh_print(char* filename, t_non* non, int sampleCount, int count, ...) {
    int j;
    float *cohr;
    float *cohi;
    float cohabs;
    FILE* out = fopen(filename, "w");

    va_list ap;
    for (int t2 = 0; t2 < non->tmax2 + 1; t2++) {
        fprintf(out, "%f ", t2 * non->deltat);
        va_start(ap, count);
        for (j = 0; j < count/2; j++) {
            cohr = va_arg(ap, float*);
            cohi = va_arg(ap, float*);
            cohabs = sqrt(cohr[t2]*cohr[t2] + cohi[t2]*cohi[t2]) / sampleCount;
            fprintf(out, "%e ", cohabs);
        }
        fprintf(out, "\n");
        va_end(ap);
    }
    fclose(out);
}

void reset_wavefn(int N, int count, ...) {
    int j;
    float* cr;
    float* ci;
    va_list ap;

    va_start(ap, count);
    for (j = 0; j < count/2; j++) {
        cr = va_arg(ap, float*);
        ci = va_arg(ap, float*);
        clearvec(cr, N);
        clearvec(ci, N);
        cr[1] = 1.0;
    }
    va_end(ap);
}

void update_trajectories(int t2, int N, float* cr, float* ci, float* pop, float* cohr, float* cohi) {
    pop[t2 + 1] += cr[N-1] * cr[N-1] + ci[N-1] * ci[N-1];
    cohr[t2 + 1] += cr[N-2] * cr[N-1] + ci[N-2] * ci[N-1];
    cohi[t2 + 1] += ci[N-2] * ci[N-1] - ci[N-1] * cr[N-2];
}

void avg_hamil(t_non* non, FILE *H_traj, float* H_avg, float* e_avg, int N) {
    int N2, nn2;
    int a, b, c;
    int L;
    float *Hamil_i_e;
    
    N2 = N * N;
    nn2 = N * (N + 1) / 2;
    clearvec(H_avg, N2);
    Hamil_i_e = (float *)calloc(nn2, sizeof(float));
    L = (non->tmax2*10 < non->length) ? non->tmax*10 : non->length;

    for (int t2 = 0; t2 < L; t2++) {
        if (read_He(non, Hamil_i_e, H_traj, t2) != 1) {
            printf("Hamiltonian trajectory file too short, could not fill buffer!\n");
            exit(1);
        }
        // Find sum of Hamiltonians
        for (a = 0; a < N; a++) {
            H_avg[a + N * a] += Hamil_i_e[a + N * a - (a * (a + 1)) / 2]; // Diagonal
            for (b = a + 1; b < N; b++) {
                H_avg[a + N * b] += Hamil_i_e[b + N * a - (a * (a + 1)) / 2];
                H_avg[b + N * a] += Hamil_i_e[b + N * a - (a * (a + 1)) / 2];
            }
        }
    }
    
    free(Hamil_i_e);
    // Divide by total number of samples
    for (a = 0; a < N; a++) {
        H_avg[a + N * a] /= non->length; // Diagonal
        for (b = a + 1; b < N; b++) {
            H_avg[a + N * b] /= non->length;
            H_avg[b + N * a] /= non->length;
        }
    }
    // Diagonalise Hamiltonian
    diagonalizeLPD(H_avg, e_avg, N);
}

void col_swap(float* a, int col1, int col2, int N) {
    float temp;
    for (int i = 0; i < N; i++) {
        temp = a[col1 + i*N];
        a[col1 + i*N] = a[col2 + i*N];
        a[col2 + i*N] = temp; 
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
    
    // Set initial values
    imax = 0, jmax = 0;
    min_diag = 0, max_offdiag = 1;

    // Check if columns must be swapped
    while (min_diag < max_offdiag) {
        if (imax != jmax) {
            col_swap(H_new, imax, jmax, N);
        }
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
    float *e_old;
    float *H_new;
    float *H_old;
    float *Hcopy;
    float *re_U;
    float *im_U;
    float *e;
    // Flags
    int nise;
    int nise_dba;
    int nise_dbb;
    int nise_dbc;
    int tnise;
    // Variables specific to NISE
    float *cr_nise;
    float *ci_nise;
    float *pop_nise;
    float *cohr_nise;
    float *cohi_nise;
    // Variables specific to NISE-DBa
    float *cr_nise_dba;
    float *ci_nise_dba;
    float *pop_nise_dba;
    float *cohr_nise_dba;
    float *cohi_nise_dba;
    // Variables specific to NISE-DBb
    float *cr_nise_dbb;
    float *ci_nise_dbb;
    float *pop_nise_dbb;
    float *cohr_nise_dbb;
    float *cohi_nise_dbb;
    // Variables specific to tNISE
    float *cr_tnise;
    float *ci_tnise;
    float *pop_tnise;
    float *cohr_tnise;
    float *cohi_tnise;

    // Set all methods to run by default
    nise = 1;
    nise_dba = 1;
    nise_dbb = 1;
    tnise = 1;

    Hamil_i_e = (float *) calloc(nn2, sizeof(float));
    H_avg = (float *) calloc(N2, sizeof(float));
    e_avg = (float *) calloc(N, sizeof(float));
    e_old = (float *) calloc(N, sizeof(float));
    H_new = (float *) calloc(N2, sizeof(float));
    H_old = (float *) calloc(N2, sizeof(float));
    Hcopy = (float *)calloc(N2, sizeof(float));
    re_U = (float *) calloc(N, sizeof(float));
    im_U = (float *) calloc(N, sizeof(float));
    e = (float *) calloc(N, sizeof(float));

    // Decide whether to use NISE
    if (nise == 1) {
        cr_nise = (float *) calloc(N, sizeof(float));
        ci_nise = (float *) calloc(N, sizeof(float));
        pop_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
        cohr_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
        cohi_nise = (float *) calloc(non->tmax2 + 1, sizeof(float));
    }
    // Decide whether to use NISE-DBa
    if (nise_dba == 1) {
        cr_nise_dba = (float *) calloc(N, sizeof(float));
        ci_nise_dba = (float *) calloc(N, sizeof(float));
        pop_nise_dba = (float *) calloc(non->tmax2 + 1, sizeof(float));
        cohr_nise_dba = (float *) calloc(non->tmax2 + 1, sizeof(float));
        cohi_nise_dba = (float *) calloc(non->tmax2 + 1, sizeof(float));
    }
    // Decide whether to use NISE-DBb
    if (nise_dbb == 1) {
        cr_nise_dbb = (float *) calloc(N, sizeof(float));
        ci_nise_dbb = (float *) calloc(N, sizeof(float));
        pop_nise_dbb = (float *) calloc(non->tmax2 + 1, sizeof(float));
        cohr_nise_dbb = (float *) calloc(non->tmax2 + 1, sizeof(float)); 
        cohi_nise_dbb = (float *) calloc(non->tmax2 + 1, sizeof(float));
    }
    // Decide whether to use tNISE
    if (tnise == 1) {
        cr_tnise = (float *) calloc(N, sizeof(float));
        ci_tnise = (float *) calloc(N, sizeof(float));
        pop_tnise = (float *) calloc(non->tmax2 + 1, sizeof(float));
        cohr_tnise = (float *) calloc(non->tmax2 + 1, sizeof(float));
        cohi_tnise = (float *) calloc(non->tmax2 + 1, sizeof(float));
    }
    
    // Determine number of samples
    sampleCount = (non->length - non->tmax2 - 1) / non->sample + 1;
    printf("Total number of samples: %i\n", sampleCount);
    // Set the initial population
    pop_nise[0] = pop_nise_dba[0] = pop_nise_dbb[0] = pop_tnise[0] = 1.0 * sampleCount;

    // Read the Hamiltonian file
    FILE* H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        exit(1);
    }

    // Find average Hamiltonian
    avg_hamil(non, H_traj, H_avg, e_avg, N);

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
        reset_wavefn(
            N, 8, 
            cr_nise, ci_nise, 
            cr_nise_dba, ci_nise_dba, 
            cr_nise_dbb, ci_nise_dbb,
            cr_tnise, ci_tnise
        );

        // Start NISE procedure
        for (int t2 = 0; t2 < non->tmax2; t2++) {
            int tm = ti + t2 + 1;
            // Copy old Hamiltonian
            copyvec(H_new, H_old, N2);
            copyvec(e, e_old, N);
            // Load new Hamiltonian
            if (read_He(non, Hamil_i_e, H_traj, tm) != 1) {
                exit(1);
            }
            build_diag_H(Hamil_i_e, H_new, e, N);

            // Check if we need to perform any swaps, to maximise overlap
            swaps(H_new, H_old, N);
            
            if (nise == 1) {
                if (!strcmp(non->basis, "Adiabatic")) {
                    // Transfer adiabatic -> site basis
                    matrix_on_vector(H_old, cr_nise, ci_nise, N);
                    // Propagate
                    propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise);
                    // Transfer site -> adiabatic basis
                    trans_matrix_on_vector(H_new, cr_nise, ci_nise, N);
                } else if (!strcmp(non->basis, "Average")) {
                    // Transfer average -> site basis
                    matrix_on_vector(H_avg, cr_nise, ci_nise, N);
                    // Propagate
                    propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise);
                    // Transfer site -> average basis
                    trans_matrix_on_vector(H_avg, cr_nise, ci_nise, N);
                } else {
                    // Propagate
                    propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise);
                }
                update_trajectories(t2, N, cr_nise, ci_nise, pop_nise, cohr_nise, cohi_nise);
            }
            
            if (nise_dba == 1) {
                if (!strcmp(non->basis, "Average")) {
                    // Transfer average -> site basis
                    matrix_on_vector(H_avg, cr_nise_dba, ci_nise_dba, N);
                    // Propagate
                    copyvec(H_old, Hcopy, N2);
                    propagate_nise_dba(non, Hcopy, H_new, e, re_U, im_U, cr_nise_dba, ci_nise_dba);
                    // Transfer site -> average basis
                    trans_matrix_on_vector(H_avg, cr_nise_dba, ci_nise_dba, N);
                } else {
                    // Propagate
                    copyvec(H_old, Hcopy, N2);
                    propagate_nise_dba(non, Hcopy, H_new, e, re_U, im_U, cr_nise_dba, ci_nise_dba);
                }
                update_trajectories(t2, N, cr_nise_dba, ci_nise_dba, pop_nise_dba, cohr_nise_dba, cohi_nise_dba);
            }

            if (nise_dbb == 1) {
                if (!strcmp(non->basis, "Adiabatic")) {
                    // Transfer adiabatic -> site basis
                    matrix_on_vector(H_old, cr_nise_dbb, ci_nise_dbb, N);
                    // Propagate
                    propagate_nise_dbb(non, H_avg, H_new, e_avg, e, re_U, im_U, cr_nise_dbb, ci_nise_dbb);
                    // Transfer site -> adiabatic basis
                    trans_matrix_on_vector(H_new, cr_nise_dbb, ci_nise_dbb, N);
                } else {
                    // Propagate
                    propagate_nise_dbb(non, H_avg, H_new, e_avg, e, re_U, im_U, cr_nise_dbb, ci_nise_dbb);
                }
                update_trajectories(t2, N, cr_nise_dbb, ci_nise_dbb, pop_nise_dbb, cohr_nise_dbb, cohi_nise_dbb);
            }
            
            if (tnise == 1) {
                if (!strcmp(non->basis, "Average")) {
                    // Transfer average -> site basis
                    matrix_on_vector(H_avg, cr_tnise, ci_tnise, N);
                    // Propagate
                    copyvec(H_old, Hcopy, N2);
                    propagate_tnise(non, Hcopy, H_new, e_old, e, re_U, im_U, cr_tnise, ci_tnise);
                    // Transfer site -> average basis
                    trans_matrix_on_vector(H_avg, cr_tnise, ci_tnise, N);
                } else {
                    // Propagate
                    copyvec(H_old, Hcopy, N2);
                    propagate_tnise(non, Hcopy, H_new, e_old, e, re_U, im_U, cr_tnise, ci_tnise);
                }
                update_trajectories(t2, N, cr_tnise, ci_tnise, pop_tnise, cohr_tnise, cohi_tnise);
            }
        }
    }
    
    char* fn_pop;
    char* fn_coh;
    // Export population and coherence data
    if (!strcmp(non->basis, "Local")) {
        fn_pop = "pop_t2_local.txt";
        fn_coh = "coh_t2_local.txt";
    } else if (!strcmp(non->basis, "Adiabatic")) {
        fn_pop = "pop_t2_adiabatic.txt";
        fn_coh = "coh_t2_adiabatic.txt";
    } else if (!strcmp(non->basis, "Average")) {
        fn_pop = "pop_t2_average.txt";
        fn_coh = "coh_t2_average.txt";
    }

    // TODO: Get rid of "magic" numbers
    pop_print(
        fn_pop, non, sampleCount, 4, 
        pop_nise, pop_nise_dba, pop_nise_dbb, pop_tnise
    );
    coh_print(
        fn_coh, non, sampleCount, 8, 
        cohr_nise, cohi_nise, 
        cohr_nise_dba, cohi_nise_dba, 
        cohr_nise_dbb, cohi_nise_dbb,
        cohr_tnise, cohi_tnise
    );

    free(Hamil_i_e);
    free(H_avg);
    free(e_avg);
    free(e_old);
    free(H_new);
    free(H_old);
    free(Hcopy);
    free(re_U);
    free(im_U);
    free(e);
    
    if (nise == 1) {
        free(cr_nise);
        free(ci_nise);
        free(pop_nise);
        free(cohr_nise);
        free(cohi_nise);
    }
    if (nise_dba == 1) {
        free(cr_nise_dba);
        free(ci_nise_dba);
        free(pop_nise_dba);
        free(cohr_nise_dba);
        free(cohi_nise_dba);
    }
    if (nise_dbb == 1) {
        free(cr_nise_dbb);
        free(ci_nise_dbb);
        free(pop_nise_dbb);
        free(cohr_nise_dbb);
        free(cohi_nise_dbb);
    }
    if (tnise == 1) {
        free(cr_tnise);
        free(ci_tnise);
        free(pop_tnise);
        free(cohr_tnise);
        free(cohi_tnise);
    }
}