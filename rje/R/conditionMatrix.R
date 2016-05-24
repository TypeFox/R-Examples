conditionMatrix <-
function (x, variables, condition = NULL, condition.value = NULL, 
    dim = NULL, incols = FALSE) 
{
    d = 2 - incols
    if (is.null(dim)) 
        dim = rep(2, log2(dim(x)[d]))
    if (prod(dim) != dim(x)[d]) 
        stop("Dimensions do not match")
    if (is.null(condition) || length(condition) == 0) 
        return(marginMatrix(x, variables, dim, incols = incols))
    if (!is.null(condition.value)) {
        if (length(condition.value) != length(condition)) 
            stop("Condition length does not match number of conditioning variables")
        if (any(condition.value > dim[condition])) 
            stop("Invalid conditioning values chosen for specified dimension")
        ord = order(condition)
        condition.value = condition.value[ord]
        condition = condition[ord]
        rest = seq_along(dim)[-condition]
        wh = sum((c(condition.value, 1) - 1) * c(1, cumprod(dim[condition]))) + 
            1
        patt = patternRepeat(seq_len(prod(dim[condition])), condition, 
            dim)
        x = x[, patt == wh]
        dim[condition] = 1
    }
    mar = sort(union(variables, condition))
    mtch = match(condition, mar)
    jt = marginMatrix(x, margin = mar, dim = dim, incols = incols)
    cond = marginMatrix(jt, margin = mtch, dim = dim[mar], incols = incols)
    patt = patternRepeat(seq_len(ncol(cond)), mtch, dim[mar])
    if (incols) {
        out = jt/c(cond[patt, ])
    }
    else {
        out = jt/c(cond[, patt])
    }
    out
}
