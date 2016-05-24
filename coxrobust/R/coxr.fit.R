
coxr.fit <- function(x, y, trunc, init, f.weight, singular.ok) {

    nrow   <- nrow(x)
    ncol   <- ncol(x)
    sorted <- order(y[,1])

    time    <- y[sorted,1]
    status  <- y[sorted,2]
    z       <- as.matrix(x[sorted,])

    skip <- integer(0)
    qrres <- qr(z)

    if ( qrres$rank < ncol ) {
       skip <- qrres$pivot[(qrres$rank+1):length(qrres$pivot)]
       notskip  <- qrres$pivot[1:qrres$rank]
    } else {
       notskip <- 1:ncol
    }

    maxabs <- double(0)
    for ( i in notskip ) {

        if (nlevels(as.factor(z[,i])) == 1) {

           skip <- c(skip, i)

        } else {

           temp <- max( abs(z[,i]) )
           if ( temp != 0 ) {
              z[,i] <- z[,i] / temp
              maxabs <- c(maxabs, temp)
           } else {
              maxabs <- c(maxabs, 1.0)
           }

        }

    }

	if ( length(skip) > 0 ) {
        
        skip <- sort(skip)

        if ( length(skip) == ncol ) {
            return(list(skip=skip))
        }

        if ( singular.ok ) {
            notskip <- 1:ncol
            notskip <- notskip[!1:ncol %in% skip]
            z       <- as.matrix(z[,notskip])
            ncol    <- ncol(z)
            init    <- init[notskip]
        } else {
            return(list(skip = skip))
        }

    }

    res      <- coxr.ple(rep(0, ncol), time, status, z, nrow, ncol)
    beta     <- res$beta
    beta.ple <- beta

    cv.ple   <- diag(1/maxabs, ncol = ncol, nrow=ncol) %*% res$hessinv %*%
                diag(1/maxabs, ncol = ncol, nrow=ncol)
    
    tmp <- as.matrix(beta-init*maxabs)
    wald.test <- t(tmp) %*% res$hess %*% tmp

    lmb <- time

    prev_ezbeta <- exp(z %*% beta)
    M <- quantile(lmb * prev_ezbeta, trunc, names = FALSE)
    lmb.ple <- as.double(cumsum( status / rev( cumsum( rev(prev_ezbeta) ) ) ))
    
    res  <- coxr.re(rep(0, ncol), lmb, status, z, prev_ezbeta, M, nrow, ncol, f.weight)
    beta <- res$beta

    for (i in 1:3) {

        lmb <- coxr.lambda(beta, lmb, status, z, prev_ezbeta, M, nrow, f.weight)

        prev_ezbeta <- exp(z %*% beta)
        M <- quantile(lmb * prev_ezbeta, trunc, names = FALSE)

        res  <- coxr.re(rep(0, ncol), lmb, status, z, prev_ezbeta, M, nrow, ncol, f.weight)
        beta <- res$beta

    }

    lmb <- coxr.lambda(beta, lmb, status, z, prev_ezbeta, M, nrow, f.weight)
    res <- coxr.covar(beta, lmb, status, z, prev_ezbeta, M, nrow, ncol, f.weight, res$hessinv)

    cv <- diag(1/maxabs, ncol = ncol, nrow = ncol) %*% res$covar %*%
          diag(1/maxabs, ncol = ncol, nrow = ncol)

    tmp <- sqrt(nrow)*as.matrix(beta-init*maxabs)
    ewald.test <- t(tmp) %*% chol2inv(chol(res$covar*nrow)) %*% tmp

    coef <- ple.coef <- double(ncol(x))
    if ( length(skip) > 0 ) {

        coef[notskip] <- beta/maxabs
        ple.coef[notskip] <- beta.ple/maxabs

        coef[skip] <- NA
        ple.coef[skip] <- NA

        var <- matrix(0, ncol(x), ncol(x))
        var[notskip, notskip] <- cv

        var.ple <- matrix(0, ncol(x), ncol(x))
        var.ple[notskip, notskip] <- cv.ple

    } else {

        coef <- beta/maxabs
        ple.coef <- beta.ple/maxabs
        var <- cv
        var.ple <- cv.ple

    }

    names(coef) <- names(ple.coef) <- dimnames(x)[[2]]

    list(
        coefficients = coef,
        ple.coefficients = ple.coef,
        lambda = lmb,
        lambda.ple = lmb.ple,
        var  = var,
        var.ple = var.ple,
        wald.test = as.double(wald.test),
        ewald.test = as.double(ewald.test),
        skip = skip)

}


coxr.ple <- function(beta, time, status, z, nrow, ncol) {

    eps <- 1e-06
    iter.max <- 100
    l <- 0.8;

    for (i in 1:iter.max) {

       res <- .C("ple", as.double(beta), as.double(time),
                as.integer(status), as.double(z), as.integer(nrow),
                as.integer(ncol), res = double(1), gradient = double(ncol),
                hessian = double(ncol*ncol))

        hess <- matrix(res$hessian, ncol, ncol)
        hessinv <- chol2inv(chol(hess))
        beta_new <- beta - l*(res$gradient %*% hessinv)
        minimum_new <- res$res

        norm <- sqrt(sum((beta - beta_new)^2))

        if ( norm <= eps ) {
            beta <- beta_new
            res <- .C("ple", as.double(beta), as.double(time),
                    as.integer(status), as.double(z), as.integer(nrow),
                    as.integer(ncol), res = double(1), gradient = double(ncol),
                    hessian = double(ncol*ncol), PACKAGE = "coxrobust")
            
            hess <- matrix(res$hessian, ncol, ncol)
            hessinv <- chol2inv(chol(hess))
            break
        } else if ( i>1 && minimum_new > minimum ) {
            beta <- (beta + beta_new) / 2
        } else {
            beta <- beta_new
            minimum <- minimum_new
        }

    }

    return(list(beta = as.double(beta), hess = hess, hessinv = hessinv))

}

coxr.re <- function(beta, time, status, z, prev_ezbeta, M, nrow, ncol, f.weight) {

    eps <- 1e-06
    iter.max <- 100
    l <- 0.8;
    
    for (i in 1:iter.max) {

        res <- .C("re", as.double(beta), as.double(time), as.integer(status),
                as.double(z), as.double(prev_ezbeta), as.double(M),
                as.integer(nrow),  as.integer(ncol), as.integer(f.weight),
                res = double(1),  gradient = double(ncol),
                hessian  = double(ncol*ncol))
        
        hess <- matrix(res$hessian, ncol, ncol)
        hessinv <- chol2inv(chol(hess))
        beta_new    <- beta - l*(res$gradient %*% hessinv)
        minimum_new <- res$res
        
        norm_re <- sqrt(sum((beta - beta_new)^2))

        if ( norm_re <= eps ) {
            beta <- beta_new
            res <- .C("re", as.double(beta), as.double(time), as.integer(status),
                    as.double(z), as.double(prev_ezbeta), as.double(M),
                    as.integer(nrow),  as.integer(ncol), as.integer(f.weight),
                    res = double(1),  gradient = double(ncol),
                    hessian  = double(ncol*ncol), PACKAGE = "coxrobust")

            hess <- matrix(res$hessian, ncol, ncol)
            hessinv <- chol2inv(chol(hess))
            break
        } else if ( i>1 && minimum_new > minimum ) {
            beta <- (beta + beta_new) / 2
        } else {
            beta <- beta_new
            minimum <- minimum_new
        }

    }

    return(list(beta = as.double(beta), hess = hess, hessinv = hessinv))

}

coxr.lambda <- function(beta, time, status, z, prev_exp_zbeta, M, nrow, f.weight) {

    cres <- .C("lambda", as.double(exp(z %*% beta)), as.double(time),
            as.integer(status), as.double(prev_exp_zbeta), as.double(M),
            as.integer(nrow), as.integer(f.weight), lmb = double(nrow))#,
            #PACKAGE = "coxrobust")

    return(cres$lmb)

}

coxr.covar <- function(beta, time, status, z, prev_exp_zbeta, M, nrow, ncol,
                  f.weight, hessinv) {

    cres <- .C("lin", as.double( exp(z %*% beta) ), as.double(time),
              as.integer(status), as.double(z), as.double(prev_exp_zbeta),
              as.double(M), as.integer(nrow), as.integer(ncol),
              as.integer(f.weight), lin = double(nrow * ncol),
              PACKAGE = "coxrobust")

    lin <- matrix(cres$lin, nrow, ncol)

    IF <- lin %*% -hessinv

    cv <- matrix(0, nrow = ncol, ncol = ncol)

    for (i in 1:nrow) {
        cv <- cv + IF[i,] %*% t(IF[i,])
    }

    return( list(covar = cv) )

}
