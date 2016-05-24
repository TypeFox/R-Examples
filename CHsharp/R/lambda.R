lambda <- function(x, y, h, d, xgrid, A, B, niterations=2) {
    lambdahat <- numeric(niterations)
    lambdas <- c(0, 1e-8*2^seq(0, 30, length = 51))
    ysharp <- y
    for (j in 1:niterations) {
        MISEhat <- NULL
        g.lp <- locpoly(x, ysharp, bandwidth = h, degree = d)
        ghat <- function(x) approx(g.lp$x, g.lp$y,  xout = x)$y
        if (j==1) {
            sigmahat2 <- (var(y - ghat(x)) + (var(y) - var(ghat(x))))
        } else {
            sigmahat2 <- var(y-ghat(x))
        }
        for (lbd in lambdas) {
            MISEhat <- c(MISEhat, MISE(x, xgrid, sigmahat2, lbd, h, ghat, 
                A, B)[3])
        }
        lambdahat[j] <- max(lambdas[which.min(MISEhat)])
        if (lambdahat[j]!=0) lambdas <- lambdahat[j]*seq(2^(-3), 2^3, length=51)       
        ysharp <- sharpen(x, y, lambdahat[j], B)
    } 
    lambdahat
}
