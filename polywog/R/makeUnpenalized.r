##
## Find which columns of the expanded design matrix should be excluded from
## the adaptive LASSO penalty
##
makeUnpenalized <- function(unpenalized, polyTerms)
{
    ## Check that all of the specified terms are actually in the model, and
    ## issue a warning if not
    notInModel <- !(unpenalized %in% rownames(polyTerms))
    if (any(notInModel)) {
        warning("The following elements of 'unpenalized' are not terms in the model: ",
                paste(unpenalized[notInModel], collapse = ", "))
    }

    ans <- rownames(polyTerms) %in% unpenalized
    names(ans) <- rownames(polyTerms)
    ans
}
