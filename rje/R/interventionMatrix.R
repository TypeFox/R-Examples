interventionMatrix <-
function (x, variables, condition, dim = NULL, incols = FALSE) 
{
    d = 2 - incols
    if (length(variables) == 0) 
        return(x)
    if (is.null(dim)) 
        dim = rep(2, log2(dim(x)[d]))
    tmp = conditionMatrix(x, variables, condition, dim = dim, 
        incols = incols)
    vars = sort(c(variables, condition))
    patt = patternRepeat(seq_len(ncol(tmp)), vars, dim)
    if (incols) 
        tmp = tmp[patt, ]
    else tmp = tmp[, patt]
    return(x/c(tmp))
}
