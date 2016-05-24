`peaks` <-
function(series, span = 3, do.pad = TRUE)
{
    if ((span = as.integer(span))%%2 != 1)
        stop("'span' must be odd")
    s1 = 1:1 + (s = span%/%2)
    if (span == 1)
        return(rep.int(TRUE, length(series)))
    z = embed(series, span)
    v = apply(z[, s1] > z[, -s1, drop = FALSE], 1, all)
    if (do.pad) {
        pad = rep.int(FALSE, s)
        return(c(pad, v, pad))
    }
    else return(v)

    
}

