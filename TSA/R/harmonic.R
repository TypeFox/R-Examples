`harmonic` <-
function (x, m = 1) 
{
    if (!is.ts(x) || (2 * m) > frequency(x)) 
        stop("x need to be a time series with 2m <= frequency(x)")
    y = outer(2 * pi * time(x), 1:m)
    cosy = apply(y, 2, cos)
    siny = apply(y, 2, sin)
    mult = 2 * (1:m)
    colnames(cosy) = paste(paste("cos(", mult, sep = ""), "*pi*t)", 
        sep = "")
    colnames(siny) = paste(paste("sin(", mult, sep = ""), "*pi*t)", 
        sep = "")
    out = cbind(cosy, siny)
    colnames(out) = c(colnames(cosy), colnames(siny))
    if ((2 * m) == frequency(x)) 
        out = out[, -(2 * m)]
    invisible(out)
}
