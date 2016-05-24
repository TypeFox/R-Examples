##
##  d r o p l e t . R  Droplet Functions for e and pi
##


dropletE <- function(n) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 1)
    if (n > 1000)
        cat("Warning: Calculating n > 1000 digits will take a long time!\n")

    n <- n + 1
    a <- numeric(n+1) + 1.0

    E <- "2."
    for (i in 1:(n-1)) {
        a <- 10 * a
        for (j in (n+1):2) {
            a[j-1] <- a[j-1] + a[j] %/% (j+1)
            a[j]   <- a[j] %%  (j+1)
        }
        E <- paste(E, a[1] %/% 2, sep="")
        a[1] <- a[1] %% 2
    }
    return(E)
}


dropletPi <- function(n) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 1)
    if (n > 1000)
        cat("Warning: Calculating n > 1000 digits will take a long time!\n")

    Pi <- ""
    p <- ""
    no_nines <- 0
    d <- n + 2
    N <- floor(10*d/3 + 1)
    a <- rep(2, N+1)

    for (l in 1:d) {
        a <- 10 * a
        for (i in (N+1):2) {
            j <- 2*i-1
            q <- a[i] %/% j; r <- a[i] %% j
            a[i] <- r
            a[i-1] <- a[i-1] + q*(i-1)
        }
        q <- a[1] %/% 10; r <- a[1] %% 10
        a[1] <- r
        if (q < 9) {
            nines <- paste(rep("9", no_nines), collapse="")
            Pi <- paste(Pi, p, nines, sep="")
            p <- q  # paste(q)
            no_nines <- 0
        } else if (q == 9) {
            no_nines <- no_nines + 1
        } else if (q == 10) {
            p <- p + 1
            nines <- paste(rep("0", no_nines), collapse="")
            Pi <- paste(Pi, p, nines, sep="")
            p <- 0
            no_nines <- 0
        } else {
            stop("dropletPi: Internal eror occurred: inform the maintainer.")
        }
    }
    Pi <- paste(substr(Pi,1,1), ".", substr(Pi,2,nchar(Pi)), sep="")
    return(Pi)
}

