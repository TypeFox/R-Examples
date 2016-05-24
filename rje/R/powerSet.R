powerSet <-
function (x, rev = FALSE) 
{
    out = list(x[c()])
    if (length(x) == 1) 
        return(c(out, list(x)))
    for (i in seq_along(x)) {
        if (rev) 
            out = c(lapply(out, function(y) c(y, x[i])), out)
        else out = c(out, lapply(out, function(y) c(y, x[i])))
    }
    out
}
