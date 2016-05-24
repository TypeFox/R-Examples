#  *****************************************************************************
#   File : invgen.R
#         ************************************************************
#   Description :
#       General inverse for a matrix
#   Version : 1.01
#   Date : 2007-10-01
#         ************************************************************
#   Author : Julien Damon <julien.damon@gmail.com>
#   License : LGPL
#   URL: https://github.com/Looping027/far
#  *****************************************************************************

#  *****************************************************************************
#   Title : invgen
#         ************************************************************
#   Description :
#       General inverse for a matrix
#   Version : 1.01
#   Date : 2007-10-01
#  *****************************************************************************
invgen <- function (a, tol = sqrt(.Machine$double.eps))
{
    if (length(dim(a)) > 2 || !(is.numeric(a) || is.complex(a)))
        stop("a must be a numeric or complex matrix")
    if (!is.matrix(a))
        a <- as.matrix(a)
    asvd <- La.svd(a)
    if (is.complex(a))
        {
            asvd$u <- Conj(asvd$u)
            asvd$v <-t(Conj(asvd$vt))
        } else {
            asvd$v <-t(asvd$vt)
        }
    Positive <- asvd$d > max(tol * asvd$d[1], 0)
    if (!any(Positive))
        array(0, dim(a)[2:1])
    else asvd$v[,Positive] %*% ((1/asvd$d[Positive]) * t(asvd$u[,Positive]))
}
