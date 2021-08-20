/* part of github.com:rqtl/qtl2scan */
/* 1-d optimization by Brent's method */
# ifndef BRENT_FMIN_H
# define BRENT_FMIN_H

/***********************************************************************
 * Taken from R ver 3.2.2
 * R-3.2.2/src/library/stats/src/optimize.c
 ***********************************************************************/

/* see brent_fmin.c */

double Brent_fmin(double ax, double bx, double guessx, double (*f)(double, void *), \
    void *fargs, double tol, unsigned int max_iter, int *info_ptr);

# endif /* BRENT_FMIN_H */
