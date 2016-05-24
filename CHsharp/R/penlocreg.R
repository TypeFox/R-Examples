penlocreg <- function(x, y, xgrid, degree = 0, h, lambda, L, ...) {
    m <- length(xgrid)
    y <- y[order(x)];  x <- sort(x); n <- length(x)
    Smoother <-function(x, degree, h, xgrid) {
    # smoother matrix -- A in the paper
    # x = original design
        A <- matrix(0, ncol=length(xgrid), nrow=length(x))
        for (j in 1:length(xgrid)){
            X <- outer(x, 0:degree, function(x,y) (x-xgrid[j])^y/factorial(y))
            K <- diag(dnorm(x-xgrid[j], sd=h))
            A[,j] <- (solve(t(X)%*%K%*%X)%*%t(X)%*%K)[1,]
        }
    A
    }
    AInterpolation <- function(A, xgrid) {
# Rows of A matrix interpolated as Functions of x gridpoints
        n <- nrow(A)  # the Smoother Matrix is A-transpose
        m <- ncol(A)
        AFun <- vector(n, mode="list")
        for (k in 1:n) {
            x <- xgrid
            y <- A[k,]  # y values evaluated at x gridpoints
            y <- c(rep(y[1], 5), y, rep(y[m], 5))
            x <- c(x[1]-(5:1), x, x[m]+1:5)
            AFun[[k]] <- splinefun(x,y)
        }
    AFun  # a list of interpolated functions for use in numerical derivative
            # calculation
    }
    penaltyMatrix <- function(n, xgrid, A, L, ...) {
        m <- length(xgrid)
        B <- matrix(0, ncol = m, nrow = n)
        AFun <- AInterpolation(A, xgrid)
        for (i in 1:n) {
            B[i,] <- L(xgrid, AFun[[i]], ...)
        }
    B
    }
    A <- Smoother(x, degree, h, xgrid)
    B <- penaltyMatrix(n, xgrid, A, L, ...)
    B <- B[, -c(1, m)]
    ysharp <- sharpen(x, y, lambda, B)
    list(x=x, y=ysharp, A=A, B=B)
}


