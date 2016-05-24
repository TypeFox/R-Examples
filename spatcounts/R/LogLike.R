`LogLike` <-
function (Yin, fm.X, region, regmodel, thinning = 1, burnin = 1){
    phi <- regmodel$phi
    omega <- regmodel$omega
    r <- regmodel$r
    be <- regmodel$beta
    ga <- regmodel$gamma
    t.i <- regmodel$t.i
    Coef <- regmodel$Coefficients
    if (fm.X == ~1) 
        xb <- matrix(1, length(Yin), 1)
    else xb <- model.matrix(fm.X)
    y <- as.matrix(Yin)
    fm.X <- cbind(region, xb)
    regindex <- fm.X[, 1]
    if (is.null(r) & is.null(omega) == TRUE) {
        m <- length(phi)
	phi <- phi[,burnin:m]
	ga <- ga[,burnin:m]
	be <- be[,burnin:m]
	n <- length(phi)
        nthin <- floor(n/thinning)
        phi <- phi[(1:nthin) * thinning]
        ga <- ga[, (1:nthin) * thinning]
        be <- be[, (1:nthin) * thinning]
        ll <- matrix(NA, dim(y)[1], nthin)
        for (i in 1:nthin) {
            mu.i <- t.i * exp(xb[, 1:dim(be)[1]] %*% be[, i] + 
                ga[regindex, i])
            ll[, i] <- log(mu.i) + (y - 1) * log(mu.i + (phi[i] - 
                1) * y) - y * log(phi[i]) - 1/phi[i] * (mu.i + 
                (phi[i] - 1) * y) - lgamma(y + 1)
        }
    }
    if (is.null(r) == FALSE) {
        n <- length(r)
	r <- r[,burnin:n]
	ga <- ga[,burnin:n]
	be <- be[,burnin:n]
        n <- length(r)
        nthin <- floor(n/thinning)
        r <- r[(1:nthin) * thinning]
        ga <- ga[, (1:nthin) * thinning]
        be <- be[, (1:nthin) * thinning]
        ll <- matrix(NA, dim(y)[1], nthin)
        for (i in 1:nthin) {
            mu.i <- t.i * exp(xb[, 1:dim(be)[1]] %*% be[, i] + 
                ga[regindex, i])
            ll[, i] <- lgamma(r[i] + y) - lgamma(r[i]) - lgamma(y + 
                1) + r[i] * log(r[i]) + y * log(mu.i) - (r[i] + 
                y) * log(r[i] + mu.i)
        }
    }
    if (is.null(omega) == FALSE) {
        n <- length(phi)
	phi <- phi[,burnin:n]
	omega <- omega[,burnin:n]
	ga <- ga[,burnin:n]
	be <- be[,burnin:n]
        n <- length(phi)
        nthin <- floor(n/thinning)
        phi <- phi[(1:nthin) * thinning]
        omega <- omega[(1:nthin) * thinning]
        ga <- ga[, (1:nthin) * thinning]
        be <- be[, (1:nthin) * thinning]
        ll <- matrix(NA, dim(y)[1], nthin)
        for (i in 1:nthin) {
            mu.i <- t.i * exp(xb[, 1:dim(be)[1]] %*% be[, i] + 
                ga[regindex, i])
            ll[, i] <- ifelse(y == 0, 1, 0) * (log(omega[i] + 
                (1 - omega[i]) * exp(-1/phi[i] * mu.i))) + ifelse(y > 
                0, 1, 0) * (log(1 - omega[i]) + log(mu.i) + (y - 
                1) * log(mu.i + (phi[i] - 1) * y) - y * log(phi[i]) - 
                1/phi[i] * (mu.i + (phi[i] - 1) * y) - lgamma(y + 
                1))
        }
    }
    LogLike.r <- list(ll = ll, Coef = Coef)
    return(LogLike.r)
}

