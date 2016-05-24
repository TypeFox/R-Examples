##' @import fda
##' @import CollocInfer
##' @importFrom graphics lines par plot points title
##' @importFrom methods as
##' @importFrom stats arima lm.fit nlminb optim sd tsdiag

nls.sparse <- function(pars, beta, active, basisvals, fdobj0, times, data, coefs, lik, proc, start, X.index, control.out, control.in, in.meth){
    if(control.out$method == "twoStage"){
        delta <- rep(1, length(data))
    }
    pars.names <- names(pars)
    f.conv <- pars.beta <- c()
    maxStep <- 8
    lambda.sparse <- control.out$lambda.sparse
    for(i in 1:control.out$maxIter){
        for(j in 1:maxStep){
            linObj <- ProfileSSE.AllPar.sparse(pars = pars, beta = beta, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth,control.in = control.in, basisvals = basisvals, fdobj0 = fdobj0)
            f.new <- sum(linObj$f^2)
            if(control.out$echo == TRUE){
                print(x = c(paste("Iter:", i, f.new)))
                cat(pars, beta, "\n")
            }
            if(i == 1){
                break
            }else{
                if(f.conv[i - 1] - f.new >= 0 & f.conv[i - 1] - f.new < control.out$tol){
                    return(list(pars=pars, beta = beta, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.beta=pars.beta, conv.message = "Converged.")))
                }
                if(f.conv[i - 1] - f.new > 0){
                    break
                }
                if(f.conv[i - 1] - f.new < 0 & j == maxStep){
                    return(list(pars=pars.old, beta = beta.old, coefs = coefs, f = f.new,  y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.beta=pars.beta, conv.message = "Non-dereasing objective.")))
                }
                pars <- 0.5*(pars - pars.old) + pars.old
                beta <- 0.5*(beta - beta.old) + beta.old
            }
        }
        pars.old <- pars
        beta.old <- beta
        f.conv <- c(f.conv, f.new)
        pars.beta <- rbind(pars.beta, c(pars, beta))
        Xdf <- - linObj$df[, 1:length(pars), drop = FALSE]
        Zdf <- - linObj$df[, (length(pars) + 1): dim(linObj$df)[2]]
        y <- - linObj$df %*% c(pars, beta) + linObj$f
        coefs <- linObj$coefs
        if(control.out$method == "penalized"){
            res <- penalized::penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda1 = lambda.sparse, positive = TRUE, trace = FALSE)
            pars <- res@unpenalized
            beta <- res@penalized
        }
        if(control.out$method == "enet"){
            res <- penalized::penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda1 = lambda.sparse, lambda2 = 0.00001, positive = TRUE, trace = FALSE)
            pars <- res@unpenalized
            beta <- res@penalized
        }
        if(control.out$method == "ols"){
            res <- lm.fit(x = cbind(Xdf, Zdf), y= y)
            pars <- res$coefficients[1:length(pars)]
            beta <- res$coefficients[(length(pars) + 1) : length(res$coefficients)]
        }
        if(control.out$method == "nnls.eq"){
            if(is.null(control.out$E)){
                E <- t(c(rep(0, length(pars)), rep(1, length(beta))))
            }
            if(is.null(control.out$F)){
                F <- 1
            }
            if(is.null(control.out$G)){
                G <- diag(length(c(pars, beta)))
            }
            if(is.null(control.out$H)){
                H <- rep(0, length(c(pars, beta)))
            }
            res <- limSolve::lsei(A= cbind(Xdf, Zdf), B = y, E = E, F=F, G = G, H = H)
            pars <- res$X[1:length(pars)]
            beta <- res$X[(length(pars) + 1) : (length(pars) + length(beta))]
        }
        if(control.out$method == "nnls"){
            res <- limSolve::nnls(A = cbind(Xdf, Zdf), B= y)
            pars <- res$x[1:length(pars)]
            beta <- res$x[(length(pars) + 1) : length(res$x)]
        }
        if(control.out$method == "nnls.old"){
            res <- nnls::nnls(A = cbind(Xdf, Zdf), b= y)
            pars <- res$x[1:length(pars)]
            beta <- res$x[(length(pars) + 1) : length(res$x)]
        }
        names(pars) <- pars.names
    }
    return(list(pars=pars.old, beta = beta.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.beta=pars.beta, conv.message = "Maximum iterations reached.")))
}


delay.fit.sparse <- function(fd0, fd.d, times, tau, beta, ndelay, basis = NULL, lik = FALSE, in.meth = "nlminb"){
    basisvals0 <- fd0$basis
    basisvals.d <- fd.d$basis
    start.d <- fd.d$basis$rangeval[1]
    y.d.beta <- matrix(0, length(times), length(ndelay))
    j <- i <- 0
    y.d.list <- list()
    if(lik == TRUE){
        for(idelay in ndelay){
            j <- j + 1
            for(itau in tau[[j]]){
                i <- i + 1
                times.d <- times - itau
                y.d <- eval.fd(times.d[times.d >= start.d], fd.d)
                if(sum(times.d < start.d)){
                    y.d0 <- eval.fd(times.d[times.d < start.d], fd0)
                    y.d <- rbind(y.d0, y.d)
                }
                y.d.beta[,j] <- y.d.beta[,j] + beta[i] * y.d[,idelay]
            }
        }
        return(list(y.d = y.d.beta))
    }
    bvals.d <- bvals.d.list <- list()
    if(is.null(basis)){
        for(idelay in ndelay){
            j <- j + 1
            for(itau in tau[[j]]){
                i <- i + 1
                times.d <- times - itau
                y.d <- eval.fd(times.d[times.d >= start.d], fd.d)
                bvals <- eval.basis(times.d[times.d >= start.d], basisvals.d, 0)
                if(sum(times.d < start.d)){
                    y.d0 <- eval.fd(times.d[times.d < start.d], fd0)
                    y.d <- rbind(y.d0, y.d)
                    bvals <- rbind(matrix(0, nrow = length(times) - dim(bvals)[1], ncol = dim(bvals)[2]), bvals)
                }
                y.d.list[[i]] <- y.d[,idelay]
                if(in.meth == "trustOptim"){
                    bvals.d.list[[i]] <- as(bvals, "CsparseMatrix")
                }
                else{
                    bvals.d.list[[i]] <- bvals
                }
                y.d.beta[,j] <- y.d.beta[,j] + beta[i] * y.d[,idelay]
                if(i == 1){
                    bvals.d[[j]] <- bvals.d.list[[i]] * beta[i]
                }
                else{
                    bvals.d[[j]] <- bvals.d[[j]] + bvals.d.list[[i]] * beta[i]
                }
            }
        }
    }
    else{
        bvals.d.list <- basis
        for(idelay in ndelay){
            j <- j + 1
            for(itau in tau[[j]]){
                i <- i + 1
                times.d <- times - itau
                y.d <- eval.fd(times.d[times.d >= start.d], fd.d)
                if(sum(times.d < start.d)){
                    y.d0 <- eval.fd(times.d[times.d < start.d], fd0)
                    y.d <- rbind(y.d0, y.d)
                }
                y.d.list[[i]] <- y.d[,idelay]
                y.d.beta[,j] <- y.d.beta[,j] + beta[i] * y.d[,idelay]
                if(i == 1){
                    bvals.d[[j]] <- bvals.d.list[[i]] * beta[i]
                }
                else{
                    bvals.d[[j]] <- bvals.d[[j]] + bvals.d.list[[i]] * beta[i]
                }
            }
        }
    }
    ## Names ??
    ## Returning only bvals.d[[1]] for now. Need to be fixed.
    return(list(y.d = y.d.beta, bvals.d = bvals.d[[1]], y.d.list = y.d.list, bvals.d.list = bvals.d.list))
}

##' Estmates spline coefficients given parameters for DDE models.
##'
##' This minimizes the objective function for DDE models defined by the addition of the \code{lik} and \code{proc} objectives with respect to the coefficients. A number of generic optimization routines can be used and some experimentation is recommended.
##' @title Inner optimization for estmating coefficients given parameters.
##' @param data Matrix of observed data values.
##' @param times Vector observation times for the data.
##' @param pars Initial values of parameters to be estimated processes.
##' @param beta Initial values of the contribution of lags for the delay.
##' @param coefs Vector giving the current estimate of the coefficients in the spline.
##' @param lik \code{lik} object defining the observation process.
##' @param proc \code{proc} object defining the state process.
##' @param in.meth Inner optimization function currently one of \code{'nlminb'}, \code{'optim'}, or \code{'trustOptim'}.
##' @param control.in Control object for inner optimization function.
##' @param basisvals Values of the collocation basis to be used. This should be a basis object from the fda package.
##' @param fdobj0  A functional data object fitted with the history part of the data.
##' @return A list with elements
##' \describe{
##'     \item{coefs}{A matrix giving he optimized coefficients.}
##'     \item{res}{The results of the inner optimization function.}
##' }
##' @export
##' @author Ziqian Zhou
inneropt.DDE <- function(data, times, pars, beta, coefs, lik, proc,
                         in.meth = "nlminb", control.in = list(),
                         basisvals, fdobj0)
{
    check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    if (in.meth == "optim") {

        if(is.null(control.in$trace)){
            control.in$trace = 0
        }

        if (is.null(control.in$maxit)) {
            control.in$maxit = 1000
        }
        if (is.null(control.in$reltol)){
            control.in$reltol = 1e-12
        }
        if (is.null(control.in$meth)) {
            control.in$meth = "BFGS"
        }
        if(is.null(control.in$reportHessian)){
            control.in$reportHessian = TRUE
        }
        imeth = control.in$meth
        control.in$meth = NULL

        res <- optim(coefs, SplineCoefsErr.DDE, gr = SplineCoefsDC.DDE,
        hessian = control.in$reportHessian, control = control.in,
        times = times, data = data, lik = lik, proc = proc, pars = pars,
        beta = beta, method = imeth, basisvals = basisvals, fdobj0 = fdobj0)

        ncoefs = matrix(res$par, ncol(lik$bvals), length(res$par)/ncol(lik$bvals))

    }
    else if (in.meth == "nlminb") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$eval.max)) {
            control.in$eval.max = 2000
        }
        if (is.null(control.in$iter.max)) {
            control.in$iter.max = 1000
        }
        if (is.null(control.in$rel.tol)) {
            control.in$rel.tol = 1e-12
        }
        if (is.null(control.in$useHessian)) {
            Hessian = SplineCoefsDC2.DDE
        }
        else {
            Hessian = NULL
        }
        ## SplineCoefsErr do not need to be changed.
        ##

        res <- nlminb(coefs, SplineCoefsErr.DDE, gradient = SplineCoefsDC.DDE,
                      hessian = Hessian, control = control.in, times = times,
                      data = data, lik = lik, proc = proc, pars = pars,
                      basisvals = basisvals, fdobj0 = fdobj0, beta = beta)

        ncoefs = matrix(res$par, ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
    }
    else if(in.meth == "trustOptim"){
        if(is.null(control.in$useHessian)){
            Hessian <- SplineCoefsDC2.DDE.sparse
            methodTO <- "Sparse"
        }
        if(is.null(control.in$report.level)){
            control.in$report.level <- -1L
        }

        res <- trustOptim::trust.optim(x = coefs, fn = SplineCoefsErr.DDE, gr = SplineCoefsDC.DDE,
                           hs = Hessian, method = methodTO, control = control.in, times = times,
                           data = data, lik = lik, proc = proc, pars = pars,
                           basisvals = basisvals, fdobj0 = fdobj0, beta = beta)
        ncoefs <- matrix(res$solution, ncol(lik$bvals), length(res$solution)/ncol(lik$bvals))
    }
    else {
        stop("Unknown optimizer specified")
    }

    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
    return(list(coefs = ncoefs, res = res))
}

##' This function runs generalized profiling for DDE models.
##' This function carry out the profiled optimization method for DDe models using a sum of squared errors criteria for both fit to data and the fit of the derivatives to a delay differential equation.
##' @title Profile Estimation Functions for DDE
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
##' @param pars Initial values of parameters to be estimated processes.
##' @param beta Initial values of the contribution of lags for the delay.
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
##' @param poslik Should the state be exponentiated before being compared to the data? When the state is represented on the log scale (posproc=TRUE), this is an alternative to taking the log of the data.
##' @param discrete Is it a discrete process?
##' @param names The names of the state variables if not given by the column names of coefs.
##' @param sparse Should sparse matrices be used for basis values? This option can save memory when using 'trust' optimization method.
##' @param basisvals0 Values of the collocation basis to be used for the history part of the data. This should be a basis object from the fda package.
##' @param coefs0  Vector giving the  estimate of the coefficients in the spline for the history part of the data.
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
##' \item{times}{The times at which the data are observed.}
##' \item{fdobj.d}{The functional data object for the estimated state process.}
##' \item{fdobj0}{The functional data object for the estimated state process of the history part.}
##' \item{tau}{The lags of delays.}
##' }
##' @examples
##' yout <- DSIRdata
##' times <- seq(-0.5, 30, by = 0.1)
##' yout0 <- yout[times >= 0, ]
##' yout.d <- yout[times >= 5, ]
##' colnames(yout.d) <-  c("S","I")
##' times0 <- times[times>=0]
##' times.d <- times[times>=5]
##' norder = 3
##' nbasis.d = length(times.d) + norder - 2
##' nbasis0 <- length(times0) + norder - 2
##' basis0 <- create.bspline.basis(range=range(times0),
##'     nbasis=nbasis0, norder=norder, breaks=times0)
##' basis.d <- create.bspline.basis(range=range(times.d),
##'     nbasis=nbasis.d, norder=norder, breaks=times.d)
##' fdnames=list(NULL,c('S', 'I'),NULL)
##' bfdPar0 = fdPar(basis0,lambda=1,int2Lfd(1))
##' bfdPar.d <- fdPar(basis.d,lambda=1,int2Lfd(1))
##' DEfd0 <- smooth.basis(times0, yout0, bfdPar0,fdnames=fdnames)$fd
##' coefs0 <- DEfd0$coefs
##' colnames(coefs0) = c("S","I")
##' initPars <- c(5, 0.0012)
##' names(initPars) <- c("gamma", "beta")
##' initBeta <- rep(0, 11)
##' initBeta[c(4,5,11)] <- c(0.611, 0.362, 0.026)
##' tau <- list(seq(0,1, length.out = 11))
##' lambda = 1000
##' DSIRfn <- DSIRfn.make()
##' \dontrun{
##' dde.fit <- Profile.LS.DDE(DSIRfn, yout.d, times.d, pars = initPars,
##'     beta = initBeta, coefs = DSIRInitCoefs, basisvals = basis.d,
##'     lambda = 1000,
##'     in.meth='nlminb', basisvals0 = basis0, coefs0 = coefs0,
##'     nbeta = length(initBeta), ndelay = 2, tau = tau,
##'     control.out = list(method = "nnls.eq", maxIter = 2, echo = TRUE))
##'     }
##' @export
##' @author Ziqian Zhou
Profile.LS.DDE <- function(fn, data, times, pars, beta, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau)
{
    if (is.null(active)) {
        active = 1:length(pars)
    }
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"
    profile.obj = LS.setup(pars = pars, coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times
    proc$more$betanames <- betanames
    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay, in.meth = in.meth)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay, in.meth = in.meth)
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
    Ires <- inneropt.DDE(data, times, pars, beta = beta,  coefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ## Ires <- IresTmp
    ncoefs <- Ires$coefs
    apars = pars[active]
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
     if (is.null(control.out$method)){
        control.out$method <- out.meth
    }
    res <- nls.sparse(pars = pars, beta = beta, active = active, basisvals = basisvals, fdobj0 = fdobj0, times = times, data = data, coefs = ncoefs, lik = lik, proc = proc, control.out = control.out, control.in = control.in, in.meth = in.meth)
    ncoefs <- res$coefs
    fdobj.d <- list(coefs = ncoefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    return(list( data = data, res = res, ncoefs = ncoefs, lik = lik, proc = proc, pars = res$pars, beta = res$beta, times = times, fdobj.d = fdobj.d, fdobj0 = fdobj0,  tau = tau))
}


ProfileSSE.AllPar.sparse <- function(pars, beta, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1,
                             basisvals, fdobj0)
{
    ## Squared Error outer criterion
    Ires = inneropt.DDE(data,times,pars, beta, coefs,lik,proc, in.meth,control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ncoefs = Ires$coefs
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj.d <- list(coefs = ncoefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta= beta, ndelay = proc$more$more$ndelay, in.meth = in.meth)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, in.meth = in.meth)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names
    ## Squared errors: No need to change for DDE
    weights <- checkweights(lik$more$weights,lik$more$whichobs,data)
    weights <- mat(weights)
    f = as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more))*sqrt(weights))
    isnaf = is.na(f)
    f[isnaf] = 0
    dlikdp = lik$more$dfdp(times,devals,pars,lik$more$more)
    dlikdp = matrix(dlikdp,dim(dlikdp)[1]*dim(dlikdp)[2],dim(dlikdp)[3])
    ## dlikdp will be zero if the likelihood doesn't directly have ode parameters,
    ## which is true for least square case.
    dlikdx = lik$more$dfdx(times,devals,pars,lik$more$more)
    ## dlikdx[i,,] is an identity matrix for every i.
    dlikdc = c()
    for(i in 1:dim(dlikdx)[2]){
        tH = c()
        for(j in 1:dim(dlikdx)[3]){
            ## Only dlikdx[,i,i] are non-zero (all 1's)
            tH = cbind(tH,as.matrix(diag(dlikdx[,i,j])%*%lik$bvals))
        }
        dlikdc = rbind(dlikdc,tH)
    }
    ## ??dlikdc: why 0.5 and 0.5 ??
    d2Hdc2  = SplineCoefsDC2sparse(ncoefs,times,data,lik,proc,pars)
    ## need not be changed?
    d2Hdcdp = SplineCoefsDCDP.sparse(ncoefs, times, data, lik, proc, pars)
    ## Got warning message:
    ## In dim(weights[whichrows, ]) == dim(diffs):
    ## longer object length is not a multiple of shorter object length?

    ## Use Implicit function theorem:
    if(is.matrix(d2Hdc2)){
        ## When will it not be a matrix? How to use solve in that case?
        dcdp = MASS::ginv(d2Hdc2) %*% d2Hdcdp
    } else {
        dcdp = as.matrix(solve(d2Hdc2,d2Hdcdp))
    }
    ## Chain rule:
    df = dlikdc%*%dcdp ## + dlikdp
    df[isnaf,] = 0
    colnames(df) = c(proc$more$parnames, proc$more$betanames)
    if(!is.null(lik$report)){ print(f) }
    f = sgn*f
    df = sgn*df
    return(list(f = f, df = df, coefs = ncoefs))
}



SplineCoefsDC.DDE <- function(coefs, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay, basis = proc$more$more$bvals.d.list)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, basis = lik$more$more$bvals.d.list)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals,
        pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars,
        proc$more)
    g = as.vector(g)
    return(sgn * g)
}

SplineCoefsDC2.DDE <- function(coefs, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay, basis = proc$more$more$bvals.d.list)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, basis = lik$more$more$bvals.d.list)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    result = as.matrix(SplineCoefsDC2sparse(coefs, times, data, lik, proc, pars, sgn))
    return(result)
}

SplineCoefsDC2.DDE.sparse <- function(coefs, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay, basis = proc$more$more$bvals.d.list)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, basis = lik$more$more$bvals.d.list)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    result <- as(SplineCoefsDC2sparse(coefs, times, data, lik, proc, pars, sgn), "dgCMatrix")
    return(result)
}

SplineCoefsErr.DDE <- function(coefs, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay, basis = proc$more$more$bvals.d.list)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, basis = lik$more$more$bvals.d.list, lik = TRUE)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    ## lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2,
        proc$bvals, pars, proc$more)
    if (!is.null(proc$report)) {
        print(f)
    }
    return(sgn * f)
}

SplineCoefsDCDP.sparse <- function (coefs, times, data, lik, proc, pars, sgn = 1)
{
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    ## Till now, H has all 0's
    H <- proc$d2fdcdp(coefs2, proc$bvals, pars, proc$more)
    return(as.matrix(sgn * H))
}


##'  This function carry out one step sparsity selection for the lags of delay given the profiled optimization result.
##' @title  Sparsity selection for the lags of delay.
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
##' @param out.meth Outer optimization selection function to be used, depending on the type of method.
##' \describe{
##' \item{"penalized"}{Uses LASSO method from \code{penalized} package.}
##' \item{"addaptive"}{Positive addaptive lasso using lars algorithm.}
##' \item{"lars"}{Positive lasso using lars algorithm.}
##' }
##' @param control.in Control object for inner optimization function.
##' @param control.out Control object for outer optimization function.
##' @param eps Finite differencing step size, if needed.
##' @param active Incides indicating which parameters of pars should be estimated; defaults to all of them.
##' @param posproc Should the state vector be constrained to be positive? If this is the case, the state is represented by an exponentiated basis expansion in the proc object.
##' @param poslik Should the state be exponentiated before being compared to the data? When the state is represented on the log scale (posproc=TRUE), this is an alternative to taking the log of the data.
##' @param names The names of the state variables if not given by the column names of coefs.
##' @param sparse Should sparse matrices be used for basis values? This option can save memory when using 'trust' optimization method.
##' @param discrete Is it a discrete process?
##' @param basisvals0 Values of the collocation basis to be used for the history part of the data. This should be a basis object from the fda package.
##' @param coefs0 Vector giving the  estimate of the coefficients in the spline for the history part of the data.
##' @param nbeta The number of lags for the delay.
##' @param ndelay A vector inidicating which state process has a delay term.
##' @param tau A list of delay lags.
##' @param nnls.res \code{res} item returned from \code{\link{Profile.LS.DDE}}
##' @return A list with elements
##' \describe{
##' \item{data}{The matrix for the observed data.}
##' \item{res}{The inner optimization result.}
##' \item{select}{A list containing the result after selection, the parameter, delay contribution and coefficients after the selection.}
##' }
##' @seealso \code{\link{Profile.LS.DDE}}
##' @export
##' @author Ziqian Zhou
sparse.DDE <- function(fn, data, times, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE,  names = NULL, sparse = FALSE, discrete = FALSE,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau, nnls.res)
{
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    pars <- nnls.res$pars
    ncoefs <- coefs <- nnls.res$coefs
    beta <- nnls.res$beta
    if (is.null(active)) {
        active = 1:length(pars)
    }
    apars <- pars[active]
    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"
    profile.obj = LS.setup(pars = pars, coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
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
    if(control.out$selection.method == "penalized"){
        Zdf <- res$Zdf
        Xdf <- res$Xdf
        y <- res$y
        lambda0 = max(abs(as.vector(t(y) %*% Zdf)))
        lambda = exp(seq(log(lambda0), log(lambda0 * 0.001), len = 20))
        pars.pen <- beta.pen <- coefs.pen <- list()
        bic <- f <- rep(NA, length(lambda))
        for(i in 1:length(lambda)){
            selection.method <- lambda[i]
            res.sparse <- penalized::penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda1 = lambda[i], positive = TRUE, trace = FALSE)
            pars.pen[[i]] <- res.sparse@unpenalized
            beta.pen[[i]] <- res.sparse@penalized
            Ires <- inneropt.DDE(data, times, pars = pars.pen[[i]], beta = beta.pen[[i]],  ncoefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (sum(beta.pen[[i]] > .Machine$double.eps) + length(pars.pen[[i]])) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        sel.res <- list(pars.pen = pars.pen[[i.select]], beta.pen = beta.pen[[i.select]], bic = bic[i.select], coefs.pen = coefs.pen[[i.select]], lambda = lambda[i.select])
    }
    ## Positive addaptive lasso using lars: Partial eliminatioin:
    if(control.out$selection.method == "addaptive"){
        y <- res$y - res$Xdf %*% res$pars
        w.beta <- abs(res$beta)                      # weights for adaptive lasso
        w.beta[w.beta == 0] <- min(w.beta[w.beta > 0]) / 2
        x <- res$Zdf
        n <- nrow(x)
        one <- rep(1, n)
        meanx <- drop(one %*% x)/n
        xc <- scale(x, meanx, FALSE)         # first subtracts mean
        normx <- sqrt(drop(one %*% (xc^2)))
        names(normx) <- NULL
        xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
        xs <- scale(xs, center=FALSE, scale=1/w.beta)  # xs times the weights
        object <- lars.pos(xs, y, type="lasso",normalize=FALSE, positive = TRUE)
        object$beta <- sweep(object$beta, 2, w.beta / normx, "*")
        bic <- rep(0, dim(object$beta)[1])
        coefs.pen <- list()
        for(i in 1:dim(object$beta)[1]){
            beta.al <- object$beta[i,]
            Ires <- inneropt.DDE(data, times, pars = res$pars, beta = beta.al,  ncoefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, res$pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (sum(beta.al > .Machine$double.eps) + length(res$pars)) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        beta.al <- object$beta[i.select,]
        sel.res <- list(pars = res$pars, beta = beta.al, bic = bic[i.select], coefs = coefs.pen[[i.select]])
    }
    ## Positive Lars:
    if(control.out$selection.method == "lars"){
        Zdf <- res$Zdf
        y <- res$y - res$Xdf %*% res$pars
        object <- lars.pos(Zdf,y, positive = TRUE)
        bic <- rep(0, dim(object$beta)[1])
        coefs.pen <- list()
        for(i in 1:dim(object$beta)[1]){
            beta.pl <- object$beta[i,]
            Ires <- inneropt.DDE(data, times, pars = res$pars, beta = beta.pl,  ncoefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (sum(beta.pl > .Machine$double.eps) +length(res$pars)) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        beta.al <- object$beta[i.select ,]
        sel.res <- list(pars = res$pars, beta = beta.al, bic = bic[i.select], coefs = coefs.pen[[i.select]])
    }
    return(list( data = data,res = res, select = sel.res))
}


SplineCoefsDC.DDE.unob <- function(coefs.unob, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0, coefs.fix){
    coefs <- cbind(coefs.unob, coefs.fix)
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals,
        pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars,
        proc$more)
    g = as.vector(g)
    g <- g[1:length(coefs.unob)]
    return(sgn * g)
}

## ordering ##
SplineCoefsErr.DDE.unob <- function(coefs.unob, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0, coefs.fix){
    fdnames <- list(NULL, NULL, NULL)
    coefs <- cbind(coefs.unob, coefs.fix)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, lik = TRUE)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    ## lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = proc$fn(coefs2, proc$bvals, pars, proc$more)
    if (!is.null(proc$report)) {
        print(f)
    }
    return(sgn * f)
}


SplineCoefsDC2.DDE.unob <- function(coefs.unob, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0, coefs.fix){
    coefs <- cbind(coefs.unob, coefs.fix)
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    result = as.matrix(SplineCoefsDC2sparse(coefs, times, data, lik, proc, pars, sgn))
    result <- result[1:length(coefs.unob), 1:length(coefs.unob)]
    return(result)
}

inneropt.DDE.unob <- function(data, times, pars, beta, coefs, lik, proc,
                         in.meth = "nlminb", control.in = list(),
                         basisvals, fdobj0, coefs.fix)
{
    check.lik.proc.data.coefs(lik, proc, data, times, cbind(coefs, coefs.fix))
    if (in.meth == "optim") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$maxit)) {
            control.in$maxit = 1000
        }
        if (is.null(control.in$reltol)){
            control.in$reltol = 1e-12
        }
        if (is.null(control.in$meth)){
            control.in$meth = "BFGS"
        }
        if (is.null(control.in$reportHessian)){
            control.in$reportHessian = TRUE
        }
        imeth = control.in$meth
        control.in$meth = NULL
        res = optim(coefs, SplineCoefsErr.DDE.unob, gr = SplineCoefsDC.DDE.unob,
        hessian = control.in$reportHessian, control = control.in,
        times = times, data = data, lik = lik, proc = proc, pars = pars,
        beta = beta, method = imeth, basisvals = basisvals, fdobj0 = fdobj0, coefs.fix = coefs.fix)
        ncoefs = matrix(c(res$par, coefs.fix), ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
    }
    else if (in.meth == "nlminb") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$eval.max)) {
            control.in$eval.max = 2000
        }
        if (is.null(control.in$iter.max)) {
            control.in$iter.max = 1000
        }
        if (is.null(control.in$rel.tol)) {
            control.in$rel.tol = 1e-12
        }
        if (is.null(control.in$useHessian)) {
            Hessian = SplineCoefsDC2.DDE.unob
        }
        else {
            Hessian = NULL
        }
        ## SplineCoefsErr do not need to be changed.
        ##
        res <- nlminb(coefs, SplineCoefsErr.DDE.unob, gradient = SplineCoefsDC.DDE.unob,
                      hessian = Hessian, control = control.in, times = times,
                      data = data, lik = lik, proc = proc, pars = pars,
                      basisvals = basisvals, fdobj0 = fdobj0, beta = beta, coefs.fix = coefs.fix)
        coefs <- c(res$par, coefs.fix)
        ncoefs = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    }
    else {
        stop("Unknown optimizer specified")
    }
    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
    return(list(coefs = ncoefs, res = res))
}


ProfileDP.sparse <- function(pars, beta, fn, data, times, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, delay = NULL, tauMax = NULL,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau)
{
    allpars <- pars
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"
    profile.obj = LS.setup(pars = pars, coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times
    proc$more$betanames <- betanames
    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay, in.meth = in.meth)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay, in.meth = in.meth)
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
    ## Ires <- IresTmp

    apars = pars[active]
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
    control.in$iter.max = 1

    res <- ProfileDP.AllPar.sparse(pars = pars, beta = beta, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth, control.in = control.in, basisvals = basisvals, fdobj0 = fdobj0)
    Xdf <- -res$Xdf
    Zdf <- -res$Zdf
    return(list(Xdf = Xdf, Zdf = Zdf))
}
