## NOTA BENE: In Matlab, there's the function balance(.) which
## calls LAPACK's  dgebal  *AND* which has the option to also return the
## transformation *diagonal* matrix D , not just the transformed matrix.
balance <- function(A, job = c("B","N", "P","S"))
    .Call("R_dgebal", A, match.arg(job))
dgebal <- balance
