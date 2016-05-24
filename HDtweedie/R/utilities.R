######################################################################
## These functions are minor modifications or directly copied from the
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.
## The original comments and header are preserved.

err <- function(n, maxit, pmax) {
    if (n == 0) 
        msg <- ""
    if (n > 0) {
        #fatal error
        if (n < 7777) 
            msg <- "Memory allocation error"
        if (n >= 20000 & n < 30000) 
            msg <- paste("Predictor ", n-20000, " has zero variance")
        if (n == 10000) 
            msg <- "All penalty factors are <= 0"
        if (n == 15000)
        	msg <- "Some observation weights are negative"
        if (n >= 30000)
        	msg <- "SVD fails"
        n <- 1
        msg <- paste("in tweediegrpnet fortran code -", msg)
    }
    if (n < 0) {
        #non fatal error
        if (n > -10000) 
            msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                maxit, " iterations; solutions for larger lambdas returned", 
                sep = "")
        if (n < -10000) 
            msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", 
                pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned", 
                sep = "")
        n <- -1
        msg <- paste("from cmd fortran code -", msg)
    }
    list(n = n, msg = msg)
}



error.bars <- function(x, upper, lower, width = 0.02, ...) {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}


getoutput <- function(fit, maxit, pmax, nvars, vnames) {
    nalam <- fit$nalam
    nbeta <- fit$nbeta[seq(nalam)]
    nbetamax <- max(nbeta)
    lam <- fit$alam[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = print(errmsg$msg, 
        call. = FALSE))
    dd <- c(nvars, nalam)
    if (nbetamax > 0) {
        beta <- matrix(fit$beta[seq(nvars * nalam)], nvars, nalam, dimnames = list(vnames, 
            stepnames))
        df <- apply(abs(beta) > 0, 2, sum)
    } else {
        beta <- matrix(0, nvars, nalam, dimnames = list(vnames, stepnames))
        df <- rep(0, nalam)
    }
    b0 <- fit$b0
    if (!is.null(b0)) {
        b0 <- b0[seq(nalam)]
        names(b0) <- stepnames
    }
    list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
}



lambda.interp <- function(lambda, s) {
    ### lambda is the index sequence that is produced by the model
    ### s is the new vector at which evaluations are required.
    ### the value is a vector of left and right indices, and a
    #   vector of fractions.
    ### the new values are interpolated bewteen the two using the
    #   fraction
    ### Note: lambda decreases. you take:
    ### sfrac*left+(1-sfrac*right)
    if (length(lambda) == 1) {
        nums <- length(s)
        left <- rep(1, nums)
        right <- left
        sfrac <- rep(1, nums)
    } else {
        s[s > max(lambda)] <- max(lambda)
        s[s < min(lambda)] <- min(lambda)
        k <- length(lambda)
        sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
        lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
        coord <- approx(lambda, seq(lambda), sfrac)$y
        left <- floor(coord)
        right <- ceiling(coord)
        sfrac <- (sfrac - lambda[right])/(lambda[left] - lambda[right])
        sfrac[left == right] <- 1
    }
    list(left = left, right = right, frac = sfrac)
}


lamfix <- function(lam) {
    llam <- log(lam)
    lam[1] <- exp(2 * llam[2] - llam[3])
    lam
} 

getmin <- function(lambda, cvm, cvsd) {
    cvmin <- min(cvm, na.rm=TRUE)
    idmin <- cvm <= cvmin
    lambda.min <- max(lambda[idmin], na.rm=TRUE)
    idmin <- match(lambda.min, lambda)
    semin <- (cvm + cvsd)[idmin]
    idmin <- cvm <= semin
    lambda.1se <- max(lambda[idmin], na.rm=TRUE)
    list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

devi <- function(y, mu, rho) {
	((y^(2-rho)-y*mu^(1-rho))/(1-rho) - (y^(2-rho)-mu^(2-rho))/(2-rho))*2
}

cvcompute=function(mat,weights,foldid,nlams){
  ###Computes the weighted mean and SD within folds, and hence the se of the mean
  wisum=tapply(weights,foldid,sum)
  nfolds=max(foldid)
  outmat=matrix(NA,nfolds,ncol(mat))
  good=matrix(0,nfolds,ncol(mat))
  mat[is.infinite(mat)]=NA#just in case some infinities crept in
  for(i in seq(nfolds)) {
    mati=mat[foldid==i,]
    wi=weights[foldid==i]
    outmat[i,]=apply(mati,2,weighted.mean,w=wi,na.rm=TRUE)
    good[i,seq(nlams[i])]=1
  }
  N=apply(good,2,sum)
  list(cvraw=outmat,weights=wisum,N=N)
}

