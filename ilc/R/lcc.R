lcc <-
function(ax, bxt, ltx, bx, kt){
    n <- length(kt)
    k <- length(ax)
    # A few simple error checks
    if (length(bxt) > 1 && length(bxt) != k)
        stop('mismatch in the b0(x) parameter length\n')
    if (length(ltx) == 1) ltx <- rep(ltx, n+k-1)
    if (length(ltx) != n+k-1) 
        stop('mismatch in the cohort-specific parameter length\n')
    if (length(bx) > 1 && length(bx) != k) 
        stop('mismatch in the b1(x) parameter length\n')
    # factor period:
    ft <- gl(n, k) 
    # factor cohort:
    ftx <- factor(as.numeric(ft)+seq(k-1, 0))
    matrix(ax + bxt*ltx[ftx] + bx*kt[ft], nrow=k, ncol=n, 
           dimnames=list(names(ax), names(kt)))
}
