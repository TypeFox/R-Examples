nls.tv <- function(pars, kappa, active, basisvals, times, data, coefs, lik, proc, start, X.index, control.out, control.in, in.meth){
    if(control.out$method == "twoStage"){
        delta <- rep(1, length(data))
    }
    pars.names <- names(pars)
    f.conv <- pars.kappa <- c()
    maxStep <- 10
    lambda.sparse <- control.out$lambda.sparse
    for(i in 1:control.out$maxIter){
        for(j in 1:maxStep){
            linObj <- ProfileSSE.AllPar.tv(pars = c(pars,kappa), times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth,control.in = control.in)
            f.new <- linObj$f
            f.new <- sum(f.new^2)
            if(control.out$method == "penalized"){
                f.new <- f.new + lambda.sparse * sum(abs(kappa[-1] - kappa[-length(kappa)]))
            }
            if(control.out$echo){
                print(x = c(paste("Iter:", i, f.new)))
                cat(pars, kappa, "\n")
            }
            if(i == 1){
                break
            }else{
                if(f.conv[i - 1] - f.new > 0 & f.conv[i - 1] - f.new < control.out$tol){
                    return(list(pars=pars.old, kappa = kappa.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa=pars.kappa, conv.message = "Converged.")))
                }
                if(f.conv[i - 1] - f.new > 0){
                    break
                }
                if(f.conv[i - 1] - f.new < 0 & j == maxStep){
                    return(list(pars=pars.old, kappa = kappa.old, coefs = coefs, f = f.new,  y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa=pars.kappa, conv.message = "Non-dereasing objective.")))
                }
                pars <- 0.5*(pars - pars.old) + pars.old
                kappa <- 0.5*(kappa - kappa.old) + kappa.old
            }
        }
        pars.old <- pars
        kappa.old <- kappa
        f.conv <- c(f.conv, f.new)
        pars.kappa <- rbind(pars.kappa, c(pars, kappa))
        Xdf <- - linObj$df[, 1:length(pars), drop = FALSE]
        Zdf <- - linObj$df[, (length(pars) +1): dim(linObj$df)[2]]
        y <- - linObj$df %*% c(pars, kappa) + linObj$f
        coefs <- linObj$coefs
        if(control.out$method == "penalized"){
            res <- penalized::penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda2 = lambda.sparse, positive = TRUE, fusedl = TRUE, trace = FALSE)
            pars <- res@unpenalized
            kappa <- res@penalized
        }
        if(control.out$method == "nnls"){
            res <- nnls::nnls(A = cbind(Xdf, Zdf), b= y)
            pars <- res$x[1:length(pars)]
            kappa <- res$x[(length(pars) + 1) : length(res$x)]
        }
        names(pars) <- pars.names
    }
    return(list(pars=pars.old, kappa = kappa.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa=pars.kappa, conv.message = "Maximum iterations reached.")))
}

Profile.LS.tv <- function(fn, data, times, pars, kappa, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL)
{
    if (is.null(active)) {
        active = 1:length(c(pars, kappa))
    }
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
    Ires <- inneropt(data, times, c(pars, kappa), coefs, lik, proc, in.meth, control.in)
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
    res <- nls.tv(pars = pars, kappa = kappa, active = active, basisvals = basisvals, times = times, data = data, coefs = ncoefs, lik = lik, proc = proc, control.out = control.out, control.in = control.in, in.meth = in.meth)
    ncoefs <- res$coefs
    return(list( data = data,res = res, ncoefs = ncoefs))
}

LS.tv <- function(fn, data, times, pars, kappa, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, nnls.res)
{
    if (is.null(active)) {
        active = 1:length(c(pars, kappa))
    }
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
    ## Ires <- IresTmp
    apars = c(pars, kappa)[active]
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
    res <- nnls.res
    ncoefs <- res$coefs
    if(is.null(control.out$pars.c))
        control.out$pars.c <- 100
    if(control.out$lambda.sparse == -1){
        Zdf <- res$Zdf
        Xdf <- res$Xdf
        y <- res$y
        pars <- res$pars
        Z0 <- rowSums(Zdf)
        beta0 <- sum((y - Xdf %*% pars ) * Z0) / sum(Z0^2)
        lambda0 <- max(abs(as.vector(t(y -Xdf %*% pars - Z0 * beta0) %*% sweep(Zdf, 1, mean(Zdf)))))
        lambda = exp(seq(log(lambda0), log(lambda0 * 0.01), len = 50))
        pars.pen <- kappa.pen <- coefs.pen <- list()
        bic <- f <- rep(NA, length(lambda))
        for(i in 1:length(lambda)){
            lambda.sparse <- lambda[i]
            res.sparse <- penalized::penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda2 = lambda.sparse, fusedl = TRUE, positive = TRUE, trace = FALSE)
            pars.pen[[i]] <- res.sparse@unpenalized
            kappa.pen[[i]] <- res.sparse@penalized
            Ires <- inneropt(data, times, pars = c(pars.pen[[i]],kappa.pen[[i]]),  ncoefs, lik, proc, in.meth, control.in)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (length(unique(kappa.pen[[i]])) + length(pars.pen[[i]])) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        sel.res <- list(pars.pen = pars.pen[[i.select]], kappa.pen = kappa.pen[[i.select]], bic = bic[i.select], coefs.pen = coefs.pen[[i.select]], lambda = lambda[i.select])
    }
    return(list(select = sel.res))
}


ProfileSSE.AllPar.tv <- function(pars, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1)
{
    ## Squared Error outer criterion
    ##    coefs = as.vector(coefs)  # First run the inner optimization
    ##
    ## f1 = SplineCoefsErr.DDE(coefs,times,data,lik,proc,pars, kappa, basisvals = basisvals, fdobj0 = fdobj0)
    ## The coefficients for delayed times:
    ##################################################
    ## Not Sure:
    ## if(!is.null(dcdp)){
    ##     tcoefs = as.vector(coefs) + dcdp%*%(pars-oldpars);
    ##     f2 = SplineCoefsErr(tcoefs,times,data,lik,proc,pars)
    ##     if(f2 < f1){
    ##         coefs = tcoefs
    ##         f1 = f2
    ##     }
    ## }
    ###################################################
    ## Inner optimization need to be changed as well.
    Ires = inneropt(data,times,pars, coefs,lik,proc, in.meth,control.in)
    ncoefs = Ires$coefs
    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names
    ## Squared errors: No need to change for DDE
    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
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
    d2Hdcdp = SplineCoefsDCDP(ncoefs, times, data, lik, proc, pars)
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
    df = dlikdc%*%dcdp + dlikdp
    df[isnaf,] = 0
    colnames(df) = c(proc$more$parnames)
    if(!is.null(lik$report)){ print(f) }
    f = sgn*f
    df = sgn*df
    return(list(f = f, df = df, coefs = ncoefs))
}



tvtrans <- function(t,k){
    month <- t %% 1
    r <- rep(0, length(t))
    for(i in 1:length(k)){
        r[month  >= (i-1)/length(k) & month < i /length(k)] <- k[i]
    }
    return(r)
}

tvDSIRfn <- list()

tvDSIRfn$fn <- function (t, y, p, more)
{
    r = y
    nKappa <- more$nKappa
    pk <- p[(length(p) - (nKappa - 1)):length(p)]
    r[, "S"] =  - tvtrans(t, pk) * y[,"I"] * y[, "S"] + 8000 * (sin(t / pi) / 2 + 2)
    r[, "I"] =  tvtrans(t, pk) * y[,"I"] * y[, "S"] - p["gamma"] * y[, "I"]
    return(r)
}

tvDSIRfn$dfdx <- function (t, y, p, more)
{
    nKappa <- more$nKappa
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    pk <- p[(length(p) - (nKappa - 1)):length(p)]
    r[, "S", "S"] = - tvtrans(t, pk) * y[,"I"]
    r[, "I", "S"] =  tvtrans(t, pk) * y[,"I"]
    r[, "I", "I"] = -p["gamma"]
    return(r)
}

tvDSIRfn$dfdp <- function (t, y, p, more)
{
    nKappa <- more$nKappa
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    r[, "I", "gamma"] = - y[, "I"]
    month <- t %% 1
    for(i in 1:nKappa){
        r[ , "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa]
    }
    return(r)
}

tvDSIRfn$d2fdx2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

tvDSIRfn$d2fdxdp <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    nKappa <- more$nKappa
    r[, "I", "I", "gamma"] = -1
    month <- t %% 1
    for(i in 1:nKappa){
        r[ , "S", "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - y[,"I"][month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "I", "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = y[,"I"][month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "S", "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - y[,"S"][month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "I", "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = y[,"S"][month  >= (i-1)/nKappa & month < i /nKappa]
    }
    return(r)
}

