##
## INPUT:
## H: Hessian matrix
## fixed: logical indicator for which elements are fixed
##
## RETURN:
## logical for whether 'H' is negative definite
## 
checkLocalID <- function(H, fixed)
{
    H <- H[!fixed, !fixed]

    ## 'chol' returns an error for non-positive definite matrices
    ans <- tryCatch(chol(-H), error = identity)
    ans <- !inherits(ans, "error")
    return(ans)
}
