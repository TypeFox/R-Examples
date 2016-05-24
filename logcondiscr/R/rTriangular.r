## function to generate random numbers
rTriangular <- function(n, a, b, c, d, e){

    ## generate random numbers according to Devroye's algorithm
    ## from the triangular density

    x <- seq(a, d, by = 1)
    dens <- dTriangular(a, b, c, d, e)
    logdens <- log(dens)

    rand <- rep(NA, n)
    m <- b     ## location of the mode
    pm <- dens[x == m]
    w <- 1 + pm / 2

    for (i in 1:n){

        cond <- TRUE

        while (identical(cond, TRUE)){
            U <- runif(1)
            W <- runif(1)
            S <- 2 * round(runif(1)) - 1
            if (U <= (w / (1 + w))){Y <- runif(1) * w / pm} else {Y <- (w + rexp(1)) / pm}
            X <- S * round(Y)

            ## now choose the correct value of the density once m + X is in {a, ..., d}
            densXm <- 0
            if (((m + X) >= a) & ((m + X) <= d)){densXm <- dens[x == (m + X)]}
            cond <- (W * min(1, exp(w - pm * Y))) > (densXm / pm)
        }

        rand[i] <- m + X
        #if (i / 1000 == round(i / 1000)){print(i)}
    }

res <- list("rand" = rand, "x" = x, "dens" = dens)
return(res)
}
