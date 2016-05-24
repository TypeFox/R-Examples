##' This function runs generalized profiling for DDE models with time varying coefficients.
##' This function carry out the profiled optimization method for DDe models using a sum of squared errors criteria for both fit to data and the fit of the derivatives to a delay differential equation with time varying coefficients.
##' @title Profile Estimation Functions for DDE with Time Varying Coefficients.
##' @param fn fn A named list of functions giving the righthand side of a delay differential equation. The functions should have arguments
##' \describe{
##' \item{times}{he times at which the righthand side is being evaluated.}
##' \item{x}{The state values at those times.}
##' \item{p}{Parameters to be entered in the system.}
##' \item{more}{A list object containing additional inputs to \code{fn}, The distributed delay state are passed into derivative calculation as \code{more$y}.}
##' The list of functions should contain the elements:
##' \item{fn}{Function to calculate the right hand sid.}
##' \item{dfdx}{Function to calculate the derivative of each right-hand function with respect to the states.}
##' \item{dfdp}{calculates the derivative of therighthand side function with respect to parameters. }
##' \item{d2fdx2}{Function to calculate the second derivatives with respect to states.}
##' \item{d2fdxdp}{Function to calculate the cross derivatives of each right-hand function with respect to state and parameters.}
##' \item{dfdx.d}{Function to calculate the the derivative of each righthand function with respect to the delayed states.}
##' \item{d2fdx.ddp}{Function to calculate the cross derivatives of each righthand function with respect to the delayed states and parameters.}
##' \item{d2fdxdx.d}{Function to calculate the cross derivatives of each right-hand function with respect to the state and the delayed states.}
##' \item{d2fdx.d2}{Function to calculate the second derivatives of the right-hand function with respect to the delayed states.}
##' }
##' @param data Matrix of observed data values.
##' @param times Vector observation times for the data.
##' @param pars Initial values of parameters to be estimated processes.
##' @param beta Initial values of the contribution of lags for the delay.
##' @param kappa Initial values of parameters for a time varying function.
##' @param coefs Vector giving the current estimate of the coefficients in the spline.
##' @param basisvals Values of the collocation basis to be used. This should be a basis object from the fda package.
##' @param lambda Penalty value trading off fidelity to data with fidelity to dif- ferential equations.
##' @param fd.obj A functional data object; if this is non-null, coefs and basisvals is extracted from here.
##' @param more An object specifying additional arguments to fn.
##' @param weights Weights for weighted estimation.
##' @param quadrature Quadrature points, should contain two elements (if not \code{NULL})
##' \describe{
##' \item{qpts}{ sQuadrature points; defaults to midpoints between knots}
##' \item{qwts}{Quadrature weights; defaults to normalizing by the length of qpts.}
##' }
##' @param in.meth Inner optimization function currently one of \code{'nlminb'}, \code{'optim'}, or \code{'trustOptim'}.
##' @param out.meth Outer optimization function to be used, depending on the type of method.
##' \describe{
##' \item{nls}{Nonlinear least square}
##' \item{nnls.eq}{Nonlinear least square with equality or/and inequality constraints of the parameters.}
##' }
##' @param control.in Control object for inner optimization function.
##' @param control.out Control object for outer optimization function.
##' @param eps Finite differencing step size, if needed.
##' @param active Incides indicating which parameters of pars should be estimated; defaults to all of them.
##' @param posproc Should the state vector be constrained to be positive? If this is the case, the state is represented by an exponentiated basis expansion in the proc object.
##' @param discrete Is it a discrete process.
##' @param poslik Should the state be exponentiated before being compared to the data? When the state is represented on the log scale (posproc=TRUE), this is an alternative to taking the log of the data.
##' @param names The names of the state variables if not given by the column names of coefs.
##' @param sparse Should sparse matrices be used for basis values? This option can save memory when using 'trust' optimization method.
##' @param basisvals0 Values of the collocation basis to be used for the history part of the data. This should be a basis object from the fda package.
##' @param coefs0 Vector giving the  estimate of the coefficients in the spline for the history part of the data.
##' @param nbeta The number of lags for the delay.
##' @param ndelay A vector inidicating which state process has a delay term.
##' @param tau A list of delay lags.
##' @return A list with elements
##' \describe{
##' \item{data}{The matrix for the observed data.}
##' \item{res}{The inner optimization result.}
##' \item{ncoefs}{The estimated coefficients.}
##' \item{lik}{The \code{lik} object generated.}
##' \item{proc}{The \code{proc} object generated.}
##' \item{pars}{The estimated parameters.}
##' \item{beta}{The estimated contribution of lags for the delay.}
##' \item{kappa}{The estimated parameters for the time varying function.}
##' \item{times}{The times at which the data are observed.}
##' \item{fdobj.d}{The functional data object for the estimated state process.}
##' \item{fdobj0}{The functional data object for the estimated state process of the history part.}
##' \item{tau}{The lags of delays.}
##' }
##' @export
##' @author Ziqian Zhou
Profile.LS.TV.DDE <- function(fn, data, times, pars, beta, kappa, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    discrete = FALSE, poslik = FALSE, names = NULL, sparse = FALSE,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau)
{
    if (is.null(active)) {
        active = 1:length(c(pars, kappa))
    }
    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"

    profile.obj = LS.setup(pars = c(pars, kappa), coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    proc$more$more$nKappa <- length(kappa)
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times

    ## Create names for delay parameters beta
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    proc$more$betanames <- betanames

    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    proc$more$more$ndelay <- lik$more$more$ndelay <- ndelay
    proc$more$more$nbeta <- lik$more$more$nbeta <- sapply(tau, length)
    proc$more$more$tau <- lik$more$more$tau <- tau
    delay <- make.delay()
    proc$dfdc <- delay$dfdc
    proc$d2fdc2 <- delay$d2fdc2.DDE
    proc$d2fdcdp <- delay$d2fdcdp.sparse
    proc$more$delay <- delay
    proc$more$dfdtau <- dfdbeta.sparse
    proc$more$d2fdxdtau <- d2fxdbeta.sparse
    proc$more$d2fdx.ddtau <- d2fdx.ddbeta.sparse


    Ires <- inneropt.DDE(data, times, c(pars, kappa), beta, coefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ## Ires <- IresTmp
    ncoefs <- Ires$coefs
    apars = c(pars, kappa)[active]
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
     if (is.null(control.out$echo)){
         control.out$echo = TRUE
    }
    res <- nls.tv.delay(pars = pars, beta = beta, kappa = kappa, active = active, basisvals = basisvals, fdobj0 = fdobj0, times = times, data = data, coefs = ncoefs, lik = lik, proc = proc, control.out = control.out, control.in = control.in, in.meth = in.meth)
    ncoefs <- res$coefs
    return(list( data = data, res = res, ncoefs = ncoefs, lik = lik, proc = proc, pars = res$pars, beta = res$beta, kappa = res$kappa, times = times, fdobj.d = fdobj.d, fdobj0 = fdobj0,  tau = tau))
}


nls.tv.delay <- function(pars, beta, kappa, active, basisvals, fdobj0, times, data, coefs, lik, proc, start, X.index, control.out, control.in, in.meth){
    if(control.out$method == "twoStage"){
        delta <- rep(1, length(data))
    }
    pars.names <- names(pars)
    kappa.names <- names(kappa)
    f.conv <- pars.kappa.beta <- c()
    maxStep <- 10
    lambda.sparse <- control.out$lambda.sparse
    for(i in 1:control.out$maxIter){
        for(j in 1:maxStep){
            linObj <- ProfileSSE.AllPar.sparse(pars = c(pars,kappa), beta = beta, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth,control.in = control.in, basisvals = basisvals, fdobj0 = fdobj0)
            f.new <- linObj$f
            f.new <- sum(f.new^2)
            if(control.out$echo){
                print(x = c(paste("Iter:", i, f.new)))
                cat(pars, kappa, beta, "\n")
            }
            if(i == 1){
                break
            }else{
                if(f.conv[i - 1] - f.new > 0 & f.conv[i - 1] - f.new < control.out$tol){
                    return(list(pars=pars.old, kappa = kappa.old, beta = beta.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa.beta=pars.kappa.beta, conv.message = "Converged.")))
                }
                if(f.conv[i - 1] - f.new > 0){
                    break
                }
                if(f.conv[i - 1] - f.new < 0 & j == maxStep){
                    return(list(pars=pars.old, kappa = kappa.old, beta = beta.old, coefs = coefs, f = f.new,  y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa.beta=pars.kappa.beta, conv.message = "Non-dereasing objective.")))
                }
                pars <- 0.5*(pars - pars.old) + pars.old
                kappa <- 0.5*(kappa - kappa.old) + kappa.old
                beta <- 0.5*(beta - beta.old) + beta.old
            }
        }
        pars.old <- pars
        kappa.old <- kappa
        beta.old <- beta
        f.conv <- c(f.conv, f.new)
        pars.kappa.beta <- rbind(pars.kappa.beta, c(pars, kappa, beta))
        Xdf <- - linObj$df[, 1:length(pars), drop = FALSE]
        Zdf <- - linObj$df[, (length(pars) +1): dim(linObj$df)[2]]
        y <- - linObj$df %*% c(pars, kappa, beta) + linObj$f
        coefs <- linObj$coefs
        if(control.out$method == "nnls.eq"){
            E <- t(c(rep(0, length(pars)) , rep(0, length(kappa)), rep(1, length(beta))))
            F <- 1
            G <- diag(length(c(pars, kappa, beta)))
            H <- rep(0, length(c(pars, kappa, beta)))
            res <- limSolve::lsei(A= cbind(Xdf, Zdf), B = y, E = E, F=F, G = G, H = H)
            ## res <- nnls(A = cbind(Xdf, Zdf), b= y)
            pars <- res$X[1:length(pars)]
            kappa <- res$X[(length(pars) + 1) : (length(pars) + length(kappa))]
            beta <- res$X[(length(pars) + length(kappa) + 1) : (length(pars) + length(kappa) + length(beta))]
        }
        names(pars) <- pars.names
        names(kappa) <- kappa.names
    }
    return(list(pars=pars.old, kappa = kappa.old, beta = beta.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa.beta=pars.kappa.beta,  conv.message = "Maximum iterations reached.")))
}

##' Sparsity selection for the lags of delay and time varying coefficients
##' This function carry out one step sparsity selection for the lags of delay given the profiled optimization result.
##' @title Sparsity selection for the lags of delay and time varying coefficients
##' @param fn A named list of functions giving the righthand side of a delay differential equation. The functions should have arguments
##' \describe{
##' \item{times}{he times at which the righthand side is being evaluated.}
##' \item{x}{The state values at those times.}
##' \item{p}{Parameters to be entered in the system.}
##' \item{more}{A list object containing additional inputs to \code{fn}, The distributed delay state are passed into derivative calculation as \code{more$y}.}
##' The list of functions should contain the elements:
##' \item{fn}{Function to calculate the right hand sid.}
##' \item{dfdx}{Function to calculate the derivative of each right-hand function with respect to the states.}
##' \item{dfdp}{calculates the derivative of therighthand side function with respect to parameters. }
##' \item{d2fdx2}{Function to calculate the second derivatives with respect to states.}
##' \item{d2fdxdp}{Function to calculate the cross derivatives of each right-hand function with respect to state and parameters.}
##' \item{dfdx.d}{Function to calculate the the derivative of each righthand function with respect to the delayed states.}
##' \item{d2fdx.ddp}{Function to calculate the cross derivatives of each righthand function with respect to the delayed states and parameters.}
##' \item{d2fdxdx.d}{Function to calculate the cross derivatives of each right-hand function with respect to the state and the delayed states.}
##' \item{d2fdx.d2}{Function to calculate the second derivatives of the right-hand function with respect to the delayed states.}
##' }
##' @param data Matrix of observed data values.
##' @param times Vector observation times for the data.
##' @param basisvals  Values of the collocation basis to be used. This should be a basis object from the fda package.
##' @param lambda Penalty value trading off fidelity to data with fidelity to dif- ferential equations.
##' @param fd.obj A functional data object; if this is non-null, coefs and basisvals is extracted from here.
##' @param more An object specifying additional arguments to fn.
##' @param weights Weights for weighted estimation.
##' @param quadrature Quadrature points, should contain two elements (if not \code{NULL})
##' \describe{
##' \item{qpts}{ sQuadrature points; defaults to midpoints between knots}
##' \item{qwts}{Quadrature weights; defaults to normalizing by the length of qpts.}
##' }
##' @param in.meth Inner optimization function currently one of \code{'nlminb'}, \code{'optim'}, or \code{'trustOptim'}.
##' @param out.meth Outer optimization selection function to be used, depending on the type of method.
##' \describe{
##' \item{"penalized"}{Uses LASSO method from \code{penalized} package.}
##' }
##' @param control.in Control object for inner optimization function.
##' @param control.out Control object for outer optimization function.
##' @param eps Finite differencing step size, if needed.
##' @param active Incides indicating which parameters of pars should be estimated; defaults to all of them.
##' @param posproc Should the state vector be constrained to be positive? If this is the case, the state is represented by an exponentiated basis expansion in the proc object.
##' @param poslik Should the state be exponentiated before being compared to the data? When the state is represented on the log scale (posproc=TRUE), this is an alternative to taking the log of the data.
##' @param discrete Is it a discrete process.
##' @param names The names of the state variables if not given by the column names of coefs.
##' @param sparse Should sparse matrices be used for basis values? This option can save memory when using 'trust' optimization method.
##' @param basisvals0 Values of the collocation basis to be used for the history part of the data. This should be a basis object from the fda package.
##' @param coefs0 Vector giving the  estimate of the coefficients in the spline for the history part of the data.
##' @param nbeta The number of lags for the delay.
##' @param ndelay A vector inidicating which state process has a delay term.
##' @param tau A list of delay lags.
##' @param nnls.res nnls.res \code{res} item returned from \code{\link{Profile.LS.DDE}}
##' @return  A list with elements
##' \describe{
##' \item{data}{The matrix for the observed data.}
##' \item{res}{The inner optimization result.}
##' \item{select}{A list containing the result after selection, the parameter, delay contribution and coefficients after the selection.}
##' }
##' @seealso \code{\link{Profile.LS.TV.DDE}}
##' @export
##' @author Ziqian Zhou
sparse.TV.DDE <- function(fn, data, times, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    basisvals0 = NULL, coefs0 = NULL,  nbeta, ndelay, tau, nnls.res)
{
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    pars <- nnls.res$pars
    beta <- nnls.res$beta
    kappa <- nnls.res$kappa
    coefs <- nnls.res$coefs
    if (is.null(active)) {
        active = 1:length(c(pars,beta, kappa))
    }
    apars <- pars[active]

    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"

    profile.obj <- LS.setup(pars = c(pars, kappa), coefs = coefs, fn = fn,
                            basisvals, lambda = lambda, fd.obj, more, data, weights,
                            times, quadrature, eps = 1e-06, posproc, poslik, discrete,
                            names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    proc$more$more$nKappa <- length(kappa)
    coefs <- profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times
    proc$more$betanames <- betanames
    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    proc$more$more$ndelay <- lik$more$more$ndelay <- ndelay
    proc$more$more$nbeta <- lik$more$more$nbeta <- sapply(tau, length)
    proc$more$more$tau <- lik$more$more$tau <- tau
    delay <- make.delay()
    proc$dfdc <- delay$dfdc
    proc$d2fdc2 <- delay$d2fdc2.DDE
    proc$d2fdcdp <- delay$d2fdcdp.sparse
    proc$more$delay <- delay
    proc$more$dfdtau <- dfdbeta.sparse
    proc$more$d2fdxdtau <- d2fxdbeta.sparse
    proc$more$d2fdx.ddtau <- d2fdx.ddbeta.sparse

    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
    if(is.null(control.out$selection.method)){
        control.out$selection.method <- out.meth
    }
    res <- nnls.res
    if(is.null(control.out$pars.c))
        control.out$pars.c <- 100
    if(control.out$lambda.sparse == "penalized"){
        if(is.null(control.out$maxInnerIter)) maxInnerIter <- 50
        if(is.null(control.out$lambda.len)) nlambda <- 10

        Zdf1 <- res$Zdf[, 1: length(kappa), drop = FALSE]
        Zdf2 <- res$Zdf[,(1+length(kappa)):(length(kappa)+length(beta)),drop = FALSE]
        Xdf <- res$Xdf
        y <- res$y
        Z1 <- rowSums(Zdf1)
        Z2 <- rowSums(Zdf2)

        beta0 <- sum((y- Xdf %*% pars) * Z1) / sum(Z1^2)
        lambda10 <- max(abs(as.vector(t(y -  Xdf %*% pars - Z1 * beta0) %*% sweep(Zdf1, 1, mean(Zdf1)))))
        lambda1 = exp(seq(log(lambda10 * 0.1), log(lambda10 * 0.001), len = nlambda))
        lambda20 <- max(abs(as.vector(t(y) %*% Zdf2)))
        lambda2 = exp(seq(log(lambda20), log(lambda20 * 0.001), len = nlambda))

        pars.pen <- kappa.pen <- beta.pen <- coefs.pen <- list()
        bic <- rep(NA, length(lambda1) * length(lambda2))
        for(i in 1:nlambda){
            for(j in 1:nlambda){
                ij <- (i - 1) * nlambda + j
                lambda.i1 <- lambda1[i]
                lambda.j1 <- lambda2[j]
                y1 <- y - Zdf2 %*% beta
                y2 <- y - Zdf1 %*% kappa
                kappa.pen[[ij]] <- kappa
                pars.pen[[ij]] <- pars
                beta.pen[[ij]] <- beta
                for(k in 1:maxInnerIter){
                    res.sparse1 <- penalized::penalized(response = y1, penalized = Zdf1, unpenalized = Xdf, lambda2 = lambda1[i], fusedl = TRUE, positive = TRUE, maxiter = 50, trace = FALSE)
                    res.sparse2 <- penalized::penalized(response = y2, penalized = Zdf2, unpenalized = Xdf, lambda1 = lambda2[j], positive = TRUE, maxiter = 50, trace = FALSE)
                    ifelse(sum(res.sparse2@penalized) == 0, sumbeta <- 1, sumbeta <- sum(res.sparse2@penalized))
                    if(sum((kappa.pen[[ij]] - res.sparse1@penalized)^2, (pars.pen[[ij]] - (res.sparse1@unpenalized + res.sparse2@unpenalized)/2)^2, (res.sparse1@unpenalized  - res.sparse2@unpenalized)^2, (beta.pen[[ij]] - res.sparse2@penalized / sumbeta)^2) < eps){
                        kappa.pen[[ij]] <- res.sparse1@penalized
                        pars.pen[[ij]] <- (res.sparse2@unpenalized + res.sparse1@unpenalized) / 2
                        beta.pen[[ij]] <- res.sparse2@penalized / sumbeta
                        break
                    }
                    kappa.pen[[ij]] <- res.sparse1@penalized
                    pars.pen[[ij]] <- (res.sparse2@unpenalized + res.sparse1@unpenalized) / 2
                    beta.pen[[ij]] <- res.sparse2@penalized / sumbeta
                    y1 <- y - Zdf2 %*% beta.pen[[ij]]
                    y2 <- y - Zdf1 %*% kappa.pen[[ij]]
                }
                Ires <- inneropt.DDE(data, times, c(pars.pen[[ij]],kappa.pen[[ij]]), beta.pen[[ij]], coefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
                devals <- as.matrix(lik$bvals%*%Ires$coefs)
                f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
                coefs.pen[[ij]] <- Ires$coefs
                sd.pen <- sd(f)
                ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
                bic[(i-1) * length(lambda1) + j] <- -2 * ll.pen + (length(unique(kappa.pen[[ij]])) + length(pars.pen[[ij]]) + sum(beta.pen[[ij]] != 0)) * log(length(data))
            }
        }
        ij.select <- which(bic == min(bic))
        sel.res <- list(pars.pen = pars.pen[[ij.select]], kappa.pen = kappa.pen[[ij.select]], beta.pen = beta.pen[[ij.select]], bic = bic[ij.select], coefs.pen = coefs.pen[[ij.select]], lambda = c(lambda1[ceiling(ij.select / nlambda)], lambda2[ifelse(ij.select %% nlambda == 0, nlambda, ij.select %% nlambda)]))
    }
    return(list(data = data, res = res, sel.res))
}

## Initialize with unobserved data
init.unob.LS.tv.delay <- function(fn, data, times, pars, beta, kappa, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, delay = NULL, tauMax = NULL,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau, unob = 1)
{
    if (is.null(active)) {
        active = 1:length(c(pars, kappa))
    }

    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"

    profile.obj = LS.setup(pars = c(pars, kappa), coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    proc$more$more$nKappa <- length(kappa)
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times

    ## Create names for delay parameters beta
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    proc$more$betanames <- betanames

    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    proc$more$more$ndelay <- lik$more$more$ndelay <- ndelay
    proc$more$more$nbeta <- lik$more$more$nbeta <- sapply(tau, length)
    proc$more$more$tau <- lik$more$more$tau <- tau
    delay <- make.delay()
    proc$dfdc <- delay$dfdc
    proc$d2fdc2 <- delay$d2fdc2.DDE
    proc$d2fdcdp <- delay$d2fdcdp.sparse
    proc$more$delay <- delay
    proc$more$dfdtau <- dfdbeta.sparse
    proc$more$d2fdxdtau <- d2fxdbeta.sparse
    proc$more$d2fdx.ddtau <- d2fdx.ddbeta.sparse
    ## Ires <- inneropt.DDE.unob(data, times, c(pars, kappa), beta, coefs[,1:unob], lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0, coefs.fix = coefs[, (unob+1):dim(coefs)[2]])
    ncoefs <- inneropt.LS.unob(data, pars, kappa, beta, coefs, lik, proc)
    return(ncoefs)
}

inneropt.LS.unob <- function(data, pars, kappa, beta, coefs, lik, proc){
    coefsS <-  coefs[,1, drop = FALSE]
    coefsI <- coefs[,2, drop = FALSE]
    I <- as.matrix(proc$bvals$bvals %*% coefsI)
    I.d <- proc$more$more$y.d
    dI <- as.matrix(proc$bvals$dbvals %*% coefsI)
    X1 <- sweep(proc$bvals$bvals, 1, tvtrans(proc$more$qpts, kappa) * I.d ,"*" ) + proc$bvals$dbvals
    X2 <- sweep(proc$bvals$bvals, 1, tvtrans(proc$more$qpts, kappa) * I.d, "*")
    X <- rbind(X1, X2)
    y <- c(proc$more$more$b, as.vector(dI + pars["gamma"] * I))
    coefs.fit <- lm.fit(x = X, y=y)
    return(coefs.fit)
}
