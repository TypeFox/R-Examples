##
##  c h i n e s e .R
##


chinese <- function(a, m) {
    stopifnot(is.numeric(a), is.numeric(m))
    if (any(!isNatural(m)) || any(floor(a) != ceiling(a)))
        stop("Arguments 'a', 'm' must be vectors of integers resp. natural numbers.")
    n <- length(m)
    if (length(a) != n)
        stop("Arguments 'a' an 'm' must be vectors of equal length.")

    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            if (GCD(m[i], m[j]) > 1) {
                stop("Elements in argument 'm' must be pairwise coprime.")
            }
        }
    }

    M <- prod(m)
    x <- 0
    for (i in 1:n) {
        Mmi <- prod(m[-i])  # or: Mmi <- M/m[i]
        mmi <- modinv(Mmi, m[i])
        x <- x + a[i] * Mmi * mmi
    }
    return(x %% M)
}
