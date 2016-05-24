rjd.fortran <- function(X, weight=NULL, maxiter=100, eps=1e-06, na.action = na.fail)
{
  msg <- "'rjd_fortran' is deprecated.\n Use 'frjd' instead.\n See
'help(rjd)' for details."
  .Deprecated("frjd", package="JADE", msg=msg)
  frjd(X, weight, maxiter, eps, na.action)
}


