condition.table <-
function (x, variables, condition = NULL, condition.value = NULL) 
{
    if (!is.null(condition.value) && length(condition) != length(condition.value)) 
        stop("condition and condition.value must have same length")
    if (length(intersect(variables, condition)) > 0) 
        stop("margin and condition must be disjoint")
    k = length(variables)
    marg = marginTable(x, c(variables, condition))
    if (length(condition) == 0) 
        return(marg/sum(marg))
    variables <- seq_len(k)
    condition <- k + seq_along(condition)
    cond <- propTable(marg, condition)
    if (is.null(condition.value)) {
        out = cond
        dim(out) = dim(x)[c(variables, condition)]
    }
    else if (is.list(condition.value)) {
        out = subtable(cond, condition, condition.value)
        dim(out) = c(dim(x)[variables], fsapply(condition.value, 
            length))
    }
    else {
        out = subtable(cond, condition, condition.value)
        dim(out) = dim(x)[variables]
    }
    return(out)
}
