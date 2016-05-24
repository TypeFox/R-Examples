##' Newey-West estimate of covariance of parameter estimates from profiling for DDE models.
##' Currently assumes a lag-5 auto-correlation among observation vectors.
##' @title ProfileSSE.covariance.DDE
##' @param pars The estimated parameters.
##' @param beta The estimated parameters.
##' @param active Incides indicating which parameters of pars should be estimated; defaults to all of them.
##' @param eps Step-size for finite difference estimate of second derivatives.
##' @param ... Additional arguments used for profiling estimation
##' @return Returns a Newey-West estimate of the covariance matrix of the parameter estimates.
##' @seealso \code{\link{Profile.LS.DDE}}
##' @export
##' @author Ziqian Zhou
ProfileSSE.covariance.DDE <- function(pars, beta, active = NULL, eps = 1e-06, ...)
{
    if (is.null(active)) {
        active = 1:length(pars)
    }
    apars <- pars[active]
    H <- matrix(0, length(apars) + length(beta), length(apars) + length(beta))
    g <- ProfileDP.sparse(pars = pars, beta = beta, active = active, ...)
    gg <- c(colSums(g$Xdf), colSums(g$Zdf))
    for(i in 1:(length(apars) + length(beta))){
        if(i <= length(apars)){
            tpars <- pars
            tpars[active][i] <-tpars[active][i] + eps
            tg <- ProfileDP.sparse(tpars, beta, active = active, ...)
            tg <- c(colSums(tg$Xdf), colSums(tg$Zdf))
        } else {
            tbeta <- beta
            tbeta[i - length(apars)] <- beta[i - length(apars)] + eps
            tbeta <- tbeta / sum(tbeta)
            tg <- ProfileDP.sparse(pars, tbeta, active = active, ...)
            tg <- c(colSums(tg$Xdf), colSums(tg$Zdf))
        }
        H[,i] <- (tg - gg)/eps
    }
    Covar <- NeweyWest.Var( 0.5*(t(H)+H) ,cbind(g$Xdf, g$Zdf) ,5)
    return(Covar)
}



ProfileDP.AllPar.sparse <- function(pars, beta, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1,
                             basisvals, fdobj0)
{
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta= beta, ndelay = proc$more$more$ndelay, in.meth = in.meth)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, in.meth = in.meth)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names
    ## Squared errors: No need to change for DDE
    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
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
    d2Hdc2  = SplineCoefsDC2sparse(coefs,times,data,lik,proc,pars)
    ## need not be changed?
    d2Hdcdp = SplineCoefsDCDP.sparse(coefs, times, data, lik, proc, pars)
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
    Xdf <- df[, 1:length(pars), drop = FALSE]
    Zdf <- df[, (length(pars) + 1): dim(df)[2]]
    return(list(Xdf = Xdf, Zdf = Zdf))
}
