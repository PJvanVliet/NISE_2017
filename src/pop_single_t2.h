#ifndef _POP_T2_
#define _POP_T2_

void pop_print(char* filename, float* pop_nise, float* pop_nise_dba, float* pop_nise_dbb, t_non* non, int sampleCount);
void coh_print(char* filename, float* cohr_nise, float* cohi_nise, float* cohr_nise_dba, float* cohi_nise_dba, float* cohr_nise_dbb, float* cohi_nise_dbb, t_non* non, int sampleCount);
void reset_wavefn(float* cr_nise, float* ci_nise, float* cr_nise_dba, float* ci_nise_dba, float* cr_nise_dbb, float* ci_nise_dbb, int N);
void propagate_NISE(t_non* non, float *H, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_nise_dba(t_non* non, float *H_old, float *H_new, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_nise_dbb(t_non* non, float *H_avg, float *H_new, float *e_avg, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_nise_dbc(t_non* non, float *H_old, float *H_new, float *e_old, float *e_new, float *re_U, float *im_U, float *cr, float *ci);
void propagate_tnise(t_non* non, float *H_old, float *H_new, float *e, float *re_U, float *im_U, float *cr, float *ci);
void thermal_correction(t_non* non, float *H_old, float *e, float *abs, float kBT);
void row_swap(float *a, int row1, int row2, int N);
void swaps(float *H_new, float *H_old, int N);
void pop_single_t2(t_non* non);

#endif // _POP_T2_