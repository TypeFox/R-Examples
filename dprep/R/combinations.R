combinations <-
function (numcol) 
{
    combine = rep(0, 0)
    variables = seq(1, numcol)
    n = variables[1]
    combine = n
    for (k in 1:(numcol - 1)) {
        ntemp = (n + (((-1)^(k + 1)) * k))
        if ((ntemp == 0) | (ntemp == numcol)) 
            n = numcol
        else n = (n + (((-1)^(k + 1)) * k))%%numcol
        combine = c(combine, n)
    }
    combinations = combine
    repet = ceiling((numcol - 1)/2) - 1
    m = seq(1, repet)
    for (j in m) {
        combine = (combine + 1)%%numcol
        combine[combine == 0] = numcol
        combinations = cbind(combinations, combine)
    }
    colnames(combinations) = NULL
    return(combinations)
}
