contr.cum <-
function(n) 
{
    if (is.numeric(n) && length(n) == 1) 
        levs <- 1:n
    else {
        levs <- n
        n <- length(n)
    }
    contr <- array(0, c(n, n), list(levs, levs))
    contr[col(contr) <= row(contr)] <- 1
    if (n < 2) stop(paste("Contrasts not defined for", n - 1, "degrees of freedom"))
    contr <- contr[, -1, drop = FALSE]
    contr
}
