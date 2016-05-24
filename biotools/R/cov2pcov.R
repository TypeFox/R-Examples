cov2pcov <- 
function(m, vars1, vars2 = seq(1, ncol(m))[-vars1])
{
    p <- (d <- dim(m))[1L]
    if (!is.numeric(m) || length(d) != 2L || p != d[2L]) 
        stop("'m' is not a square numeric matrix")
    Is <- sqrt(1/diag(m))
    if (any(!is.finite(Is))) 
        warning("diag(.) had 0 or NA entries; non-finite result is doubtful")
    pcov <- m[vars1, vars1] - 
       m[vars1, vars2] %*% solve(m[vars2, vars2], m[vars2, vars1])
    return(pcov)
}
