CovNAOgk <- function(x,
                     niter = 2,
                     beta = 0.9,
                     impMeth = c("norm" , "seq", "rseq"),
                     control)
{
    ## Analize and validate the input parameters ...

    ## If a control object was supplied, take the option parameters from it,
    ##  but if single parameters were passed (not defaults) they will override the
    ##  control object.
    ## The functions 'mrob()' and 'vrob()' can be supplied only via the control
    ##  object. If no control object is passed these function will be taken
    ##  from the default one

    defcontrol <- CovControlOgk()          # default control
    mrob <- defcontrol@mrob
    vrob <- defcontrol@vrob
    if(!missing(control)){                  # a control object was supplied
        if(niter == defcontrol@niter)       niter <- control@niter
        if(beta == defcontrol@beta)         beta <- control@beta
        mrob <- control@mrob
        vrob <- control@vrob
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    call <- match.call()


    ## drop all rows which contain only missings
    na.x <- rowSums(ifelse(is.na(x),1,0)) == ncol(x)
    ok <- !na.x
    x <- x[ok, , drop = FALSE]

    dimn <- dimnames(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    if(p < 2)
        stop("Need at least 2 columns ")

##    s <- prelim.norm(x)                     # do preliminary manipulations
##    thetahat <- em.norm(s, showits=FALSE)   # find the mle estimates
##    rngseed(1234567)                        # set random number generator seed
##    ximp <- imp.norm(s, thetahat, x)        # impute missing data under the MLE
##    xx<-imp.norm(s, thetahat, x)

    ximp <- .imputation(x, impMeth = impMeth)
    ogk <- CovOgk(ximp, niter = niter, beta = beta, control=control)
    ogk.rew <- .cov.na.rew(x, center=ogk@raw.center, cov=ogk@raw.cov, cutoff=beta, method="ML")

    raw.mah <- .mah.na(x, ogk@raw.center, ogk@raw.cov)
    raw.wt <- .dflag(raw.mah, pr=0.975)

#    dist2 <- opw$distances
#    quantiel <- qchisq(beta, p)
#    qq <- (quantiel * median(dist2))/qchisq(0.5, p)
#    wt <- ifelse(dist2 < qq, 1, 0)
#    sum.wt <- sum(wt)

##  compute the reweighted estimates:  OGK2(0.9)
#    wcenter <- colSums(x*wt)/sum.wt
#    X <- sqrt(wt) * sweep(x, 2, wcenter)
#    wcov <- (t(X) %*% X)/sum.wt

##  Compute consistency correction factor for the reweighted  cov
##    qdelta.rew <- qchisq(sum(wt)/n, p)
##    cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum(wt)/n)
##    cnp2 <- 1/cdeltainvers.rew


    method="Orthogonalized Gnanadesikan-Kettenring Estimator for incomplete data"
    ans <- new("CovNAOgk",
               call = call,
               iter=niter,
               crit=1,
               cov=ogk.rew$cov,
               center=ogk.rew$center,
               n.obs=ogk@n.obs,
               raw.cov=ogk@raw.cov,
               raw.center=ogk@raw.center,
               raw.mah = raw.mah,
               raw.wt = raw.wt,
               X = x,
               method=method)
    ans
}
