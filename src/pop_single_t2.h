#ifndef _POP_T2_
#define _POP_T2_

void pop_print(char* filename, float* pop_nise, float* pop_prezhdo, float* pop_alt, t_non* non, int sampleCount);
void reset_wavefn(float* cr_nise, float* ci_nise, float* cr_prezhdo, float* ci_prezhdo, float* cr_alt, float* ci_alt, int N);
void propagate_NISE(t_non* non, float *H, float *e, float *re_U, float *im_U, float *cr, float *ci, float *pop_nise, float *pop_nise_ad, int t2);
void propagate_prezhdo(t_non* non, float *H_old, float *H_new, float *e, float *re_U, float *im_U, float *cr, float *ci, float *pop_prezhdo, float *pop_prezhdo_ad, int t2);
void propagate_alt(t_non* non, float *H_old, float *H_new, float *e, float *re_U, float *im_U, float *cr, float *ci, float *pop_alt, float *pop_alt_ad, int t2);
void row_swap(float *a, int row1, int row2, int N);
void swaps(float *H_new, float *H_old, int N);
void pop_single_t2(t_non* non);

#endif // _POP_T2_