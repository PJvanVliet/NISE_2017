#ifndef _PROP_SCHEMES_
#define _PROP_SCHEMES_

void propagate_NISE(t_non* non, float *H, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_nise_dba(t_non* non, float *H_old, float *H_new, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_nise_dbb(t_non* non, float *H_avg, float *H_new, float *e_avg, float *e, float *re_U, float *im_U, float *cr, float *ci);
void propagate_tnise(t_non* non, float *H_old, float *H_new, float *e_old, float *e_new, float *re_U, float *im_U, float *cr, float *ci);
void thermal_correction(t_non* non, float *H_old, float *e, float *abs);

#endif // _PROP_SCHEMES_