

`nsamp` <- function (n, k, replace = FALSE, ordered = FALSE){
    m <- length(n)
    if (length(k) != m) 
        stop("number of urns doesn't equal number of sample sizes")
    if (length(replace) != m) {
        replace <- rep(replace, length.out = m)
    }
    if (length(ordered) != m) {
        ordered <- rep(ordered, length.out = m)
    }
    res <- c()
    for (i in 1:m) if (isTRUE(replace[i])) {
        if (isTRUE(ordered[i])) {
            res[i] <- n[i]^k[i]
        }
        else {
            res[i] <- choose(n[i] - 1 + k[i], k[i])
        }
    }
    else {
        if (isTRUE(ordered[i])) {
            res[i] <- factorial(n[i])/factorial(n[i] - k[i])
        }
        else {
            res[i] <- choose(n[i], k[i])
        }
    }
    return(res)
}



`permsn` <- function (x, m)
{

    # require(combinat)
    if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x)

        x <- seq(x)
    temp <- combn(x, m)
    if ( isTRUE(all.equal(m,1)) ) {

        P <- temp
    } else if (isTRUE(all.equal(m, length(x)))) {

        temp <- matrix(x, ncol = 1)
        P <- array(unlist(permn(temp[, 1])), dim = c(m, factorial(m)))
    } else {
        k <- dim(temp)[1]
        n <- dim(temp)[2]
        P <- array(unlist(permn(temp[, 1])), dim = c(k, factorial(k)))
        for (i in 2:n) {
            a <- temp[, i]
            perms <- array(unlist(permn(a)), dim = c(k, factorial(k)))
            P <- cbind(P, perms)
        }


    }
    return(P)
}
