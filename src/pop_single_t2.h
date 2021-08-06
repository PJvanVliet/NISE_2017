#ifndef _POP_T2_
#define _POP_T2_

void pop_print(char* filename, float* pop_nise, float* pop_prezhdo, float* pop_alt, t_non* non, int sampleCount);
void coh_print(char* filename, float* cohr_nise, float* cohi_nise, float* cohr_standard, float* cohi_standard, float* cohr_harmonic, float* cohi_harmonic, t_non* non, int sampleCount);
void reset_wavefn(float* cr_nise, float* ci_nise, float* cr_prezhdo, float* ci_prezhdo, float* cr_alt, float* ci_alt, int N);
void propagate_NISE(t_non* non, float *H, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_standard(t_non* non, float *H_old, float *H_new, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_harmonic(t_non* non, float *H_avg, float *H_new, float *e_avg, float *e, float *re_U, float *im_U, float *cr, float *ci);
void row_swap(float *a, int row1, int row2, int N);
void swaps(float *H_new, float *H_old, int N);
void pop_single_t2(t_non* non);

#endif // _POP_T2_