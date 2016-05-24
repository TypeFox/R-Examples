##
##  l e b e s g u e . R  Lebesgue Coefficient
##


lebesgue <- function(x, refine = 4, plotting = FALSE) {
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector.")
    if (! refine %in% c(2,3,4,5,6,7,8,9,10))
        stop("Argument 'refine' must be one of 2,3,4,5,6,7,8,9,10.")

    n <- length(x)
    N <- 2^refine * n + 1
    X <- matrix(rep(x, times = n), nrow = n, ncol = n)

    # weights
    w <- 1 / apply(X - t(X) + diag(n), 1, prod)

    # refine grid points
    xp <- seq(min(x), max(x), length.out = N)
    xdiff <- matrix(rep(xp, times = n), n, N, byrow = TRUE) -
             matrix(rep(x,  times = N), n, N)

    inds <- (xdiff == 0)
    lfun <- apply(xdiff, 2, prod)
    xdiff[inds] <- .Machine$double.eps

    # compute Lebesgue function
    Y <- abs((diag(w) %*% matrix(rep(lfun, times = n), n, N, byrow = TRUE)) / xdiff)
    lebfun <- apply(Y, 2, sum)
    if (plotting) {
        plot(xp, lebfun, type = "l", col = "blue", lty = 2, lwd = 2,
             xlab="Grid points", ylab="Coefficients", main = "Lebesgue Function")
        grid()
    }

    # return Lebesgue coefficient
    return(max(lebfun))
}
