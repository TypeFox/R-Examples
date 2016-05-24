rationalfit <- function(x, y, d1 = 5, d2 = 5) {
    stopifnot(is.numeric(x), is.numeric(y))
    n <- length(x)
    if (n <= 2)
        stop("Length of arguments 'x' and 'y' must be greater than 2.")
    if (length(y) != n)
        stop("Arguments 'x' ans 'y' must be of the same length.")

    if (is.unsorted(x))
        stop("Argument 'x' must be a sorted vector")

    p <- finds(!is.finite(y))

    dinf <- c()
    while (length(p) > 0) {
       y <- y * (x - x[p[1]])       # adjust remaining y values
       y <- y[-p[1]]                # remove bad y value, now a NaN
       dinf <- c(dinf, x[p[1]])     # remember where pole was   
       x <- x[-p[1]]                # now remove that x value too
       if (d2 > 0) d2 <- d2 - 1     # reduce expected order of den
       p <- finds(!is.finite(y))     # have all Inf values been removed yet?
    }

    yy <- length(y)                 # x and y have a new length
    an <- outer(x, d1:0, "^")       # vandermonde matrix
    ad <- outer(x, d2:0, "^")
    for (k in 1:yy)
        ad[k, ] <- y[k] * ad[k, ]

    # A is basically N-y*D
    A <- cbind(an, -ad)             #  LS solution is in the null space of A
    
    V <- svd(A, nv = ncol(A))$v     # [u,s,v]=svd(A); % null space is in the cols of V
    ND <- V[, ncol(A)]              # use the "most null" vector
    
    N <- ND[1:(d1+1)]
    D <- ND[(d1+2):length(ND)]
        
    D1 <- D[1]
    if (D1 == 0) D1 <- 1
    N <- N/D1
    D <- D/D1
    D <- Poly(c(dinf, roots(D)))    # and then add the removed +/- Inf poles back in

    eps <- .Machine$double.eps      # remove small imaginary parts
    if (all(Im(D) < eps)) D <- Re(D)
    
    maxd <- max(abs(D))             # normalize max Den value to be +/- 1
    if (maxd == 0) maxd <- 1
    N <- N/maxd
    D <- D/maxd

    return(list(p1 = N, p2 = D))
}
