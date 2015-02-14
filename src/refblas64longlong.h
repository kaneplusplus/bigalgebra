extern "C"
{
  void int8_dgemm (char *, char *, long long *, long long *,
	       long long *, double *, double *, long long *,
	       double *, long long *, double *, double *,
	       long long *);

  void int8_daxpy (long long *, double *, double *, long long *,
	       double *, long long *);
/*
  void dcopy (long long *, double *, long long *, double *,
	       long long *);
  void int8_dscal (long long *, double *, double *, long long *);
  void int8_dgeev (char *, char *, long long *, double *, long long *, 
               double *, double *, double *, long long *, double *,
               long long *, double *, long long *, long long *);
*/
  void int8_dpotrf (char *, long long *, double *, long long *, long long *);
/*
  void int8_dgeqrf (long long *, long long *, double *, long long *, double *,
                double *, long long *, long long *);
  void int8_dgesdd (char *, long long *, long long *, double *, long long *,
                double *, double *, long long *, double *,
                long long *, double *, long long *, long long *,
                long long *);
*/
}
