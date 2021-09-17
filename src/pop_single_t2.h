#ifndef _POP_T2_
#define _POP_T2_

void pop_print(char* filename, t_non* non, int sampleCount, int count, ...);
void coh_print(char* filename, t_non* non, int sampleCount, int count, ...);
void reset_wavefn(t_non* non, int N, int count, ...);
void update_trajectories(t_non* non, int t2, int N, float *cr, float *ci, float *pop, float *cohr, float *cohi);
void avg_hamil(t_non* non, FILE* H_traj, float *H_avg, float *e_avg, int N);
void col_swap(float *a, int row1, int row2, int N);
void swaps(float *H_new, float *H_old, int N);
void pop_single_t2(t_non* non);

#endif // _POP_T2_
