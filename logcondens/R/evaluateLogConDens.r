evaluateLogConDens <- function(xs, res, which = 1:5, gam = NULL, print = FALSE){

    x <- res$x
    w <- res$w
    phi <- res$phi
    Fhat <- res$Fhat
    IsKnot <- res$IsKnot
    n <- length(x)
    values.mat <- matrix(NA, ncol = 6, nrow = length(xs))

    ## quantities for smooth estimates
    if (max(which) > 3){
        if (is.null(gam) == TRUE){  
            VarFn <- LocalVariance(x = x, w = w, phi = phi)
            gam <- sqrt(res$sig ^ 2 - VarFn)
        }

        js <- 2:n
        xj <- x[js - 1]
        xj1 <- x[js]
        f <- exp(phi)
        a <- c(NA, diff(phi) / diff(x))[js]
    }


## now compute values of functions for each x0 in xs
for (i in 1:length(xs)){

    values <- rep(NA, 5)
    x0 <- xs[i]

    if (x0 < x[1]){values[1:3] <- c(-Inf, 0, 0)}
    if (x0 == x[1]){values[1:3] <- c(phi[1], exp(phi[1]), 0)}
    if (x0 > x[n]){values[1:3] <- c(-Inf, 0, 1)}
    
    if (x0 > x[1] && x0 <= x[n]){
        x.knot <- res$knots
        phi.knot <- phi[IsKnot > 0]
        k <- length(x.knot[x.knot < x0])
        if ((1 %in% which) | (2 %in% which)){
            phi.x0 <- (1 - (x0 - x.knot[k]) / (x.knot[k + 1] - x.knot[k])) * phi.knot[k] + (x0 - x.knot[k]) / (x.knot[k + 1] - x.knot[k]) * phi.knot[k + 1]
            f.x0 <- exp(phi.x0)
            values[1:2] <- c(phi.x0, f.x0)
            }
        
        if (3 %in% which){
            j <- length(x[x < x0])
            Fhat.x0 <- Fhat[j] + (x[j + 1] - x[j]) * J00(phi[j], phi[j + 1], (x0 - x[j]) / (x[j + 1] - x[j]))
            values[3] <- Fhat.x0
            }
    }

    ## compute quantities that appear in smooth PDF and CDF only   
    if (max(which) > 3){qs <- Q00(x = x0, a = a, u = xj, v = xj1, gamma = gam, QFhat = (5 %in% which))}
    
    ## smooth PDF
    if (4 %in% which){values[4] <- sum(f[js - 1] * qs$q)}

    ## smooth CDF
    if (5 %in% which){values[5] <- sum(f[js - 1] * qs$Q)}

    values.mat[i, ] <- c(x0, values)
    
    i10 <- i / length(xs) * 10
    if ((round(i10) == i10) & (print == TRUE)){print(paste(i10 * 10, "% of computation of smooth estimates done", sep = ""))}
}
    
    values.mat[, c(FALSE, (1:5 %in% which) == FALSE)] <- NA
    colnames(values.mat) <- c("xs", "log-density", "density", "CDF", "smooth.density", "smooth.CDF")   
    return(values.mat)
}
