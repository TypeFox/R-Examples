##
## Remove all intercepts from a formula with a single left-hand side and one or
## more right-hand sides.
##
## ARGUMENTS:
##   formula.: an object of class "Formula"
##   mf: the model frame (included to deal with the case where 'formula.' is of
##       the form y ~ .)
##
## RETURN:
##   'formula.', with each rhs term including "- 1"
##
removeIntercepts <- function(formula., mf)
{
    nrhs <- length(formula.)[2]
    flist <- vector("list", nrhs)

    ## Keep LHS for first term [the first line is to ensure no error if the
    ## formula is of the form y ~ .]
    ft <- terms(formula(formula., lhs = 1, rhs = 1), data = mf)
    flist[[1]] <- update(formula(ft), . ~ . - 1)

    ## Don't use LHS in subsequent terms (if any)
    for (i in seq_len(nrhs-1)) {
        ft <- terms(formula(formula., lhs = 0, rhs = i+1), data = mf)
        flist[[i+1]] <- update(formula(ft), ~ . - 1)
    }

    ans <- do.call(as.Formula, flist)
    return(ans)
}
