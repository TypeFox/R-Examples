### R code from vignette source 'variance.Rnw'

###################################################
### code chunk number 1: variance.Rnw:19-20
###################################################
options(continue=' ', width=60)


###################################################
### code chunk number 2: variance.Rnw:56-58 (eval = FALSE)
###################################################
## fit <- coxme(Surv(time, status) ~ stage + histology + (z|1),
##              data=mydata, varlist=gexchange)


###################################################
### code chunk number 3: variance.Rnw:65-72 (eval = FALSE)
###################################################
## gexchange <- function() {
##     out <- list(initialize = geinit,
##                 generate  = gegenerate,
##                 wrapup    = gewrapup)
##     class(out) <- "coxmevar"
##     out
## }


###################################################
### code chunk number 4: variance.Rnw:114-161
###################################################
geinit <- function(vinit, vfixed, intercept, G, X, sparse) {
    if (length(G) >0)  
        stop("This function does not handle groupings")
    if (length(X)==0)
        stop("No covariates given!")
    if (ncol(X) <2)
        stop("Only one random slope, correlation is irrelevant")
    
    if (length(vinit) >0) {
        if (length(vinit) !=2) stop("Wrong length for vinit")
        indx <- match(names(vinit), c("sigma", "rho"))
        if (any(is.na(indx))) 
            stop("Unrecognized parameter name in vinit values")
                      
        theta <- list(vinit[[indx[1]]], vinit[[indx[2]]])
    }
    else theta <- list(c(.05, .2, .6)^2, c(.01, .5, .9))
        
    which.fixed <- c(FALSE, FALSE)    
    if (length(vfixed) >0) {
        indx <- match(names(vfixed), c("sigma", "rho"))
        if (any(is.na(indx))) 
            stop("Unrecognized parameter name in vfixed values")
        if (is.list(vfixed)) {
            if (any(sapply(vfixed, length) !=1))
                stop("Any fixed parameter must have a single value")
            vfixed <- unlist(vfixed)
        }
        which.fixed[indx] <- TRUE
        theta[which.fixed] <- vfixed[which.fixed]
    }

    if (any(theta[[1]] <=0)) stop("Variance must be >0")
    if (any(theta[[2]] <0) || any(theta[[2]] >=1))
        stop("Correlation must be between 0 and 1")

    theta[[2]] <- theta[[2]]/(1-theta[[2]])
    
    # In the shape of X, first column=1, second =2, etc
    xmap <- matrix(rep(1:ncol(X), each=nrow(X)), nrow(X))

    list(theta=lapply(theta, log)[!which.fixed], imap=NULL, 
         X=X, xmap=xmap, penalty=FALSE,
         parms=list(theta=sapply(theta, function(x) x[1]), 
                    fixed=which.fixed,
                    xname=dimnames(X)[[2]], nvar=ncol(X)))
}


###################################################
### code chunk number 5: generate
###################################################
gegenerate <- function(newtheta, parms) {
    safe.exp <- function(x, emax=20) {
        exp(pmax(-emax, pmin(emax, x)))
    }

    theta <- parms$theta
    if (!all(parms$fixed))
        theta[!parms$fixed] <- safe.exp(newtheta)

    correlation <- theta[2]/(1+theta[2])
    variance    <- min(theta[1],20)   #keep it out of trouble
    varmat <- matrix(variance*correlation, nrow=parms$nvar, 
                     ncol=parms$nvar)    
    diag(varmat) <- variance
    varmat
}


###################################################
### code chunk number 6: gewrapup
###################################################
gewrapup <- function(newtheta, b, parms) {
    theta <- parms$theta
    theta[!parms$fixed] <- exp(newtheta)
    correlation <- theta[2]/(1+theta[2])

    rtheta <- list('(Shrinkage)' = c(variance=theta[1], 
                                     correlation=correlation))
    names(b) <- parms$xname
    list(theta=rtheta, b=list(X=list(b)))
}


###################################################
### code chunk number 7: example
###################################################
gexchange <- function() {
    out <- list(initialize = geinit,
                generate  = gegenerate,
                wrapup    = gewrapup)
    class(out) <- "coxmevar"
    out
}
require(coxme)
set.seed(1960)
n <- nrow(lung)
dgene <- matrix(rnorm(n*12, mean=8, sd=2), ncol=12)
fit <- coxme(Surv(time, status) ~ age + ph.ecog + (dgene |1),
             data=lung, varlist=gexchange)
print(fit)

coxme(Surv(time, status) ~ age + ph.ecog + (dgene |1),
             data=lung, varlist=gexchange,
             vfixed=c(sigma=.01))


###################################################
### code chunk number 8: variance.Rnw:338-340 (eval = FALSE)
###################################################
## coxme(Surv(time, status) ~ stage + trt + (1+trt | site), data=mydata,
##       varlist=myvar(c(.1, .2, .4)))


###################################################
### code chunk number 9: myvar
###################################################
myvar <- function(var) {
    initialize <- function(vinit, fixed, intercept, G, X, sparse) {
        imap <- as.matrix(as.numeric(G[[1]]))
        ngroup <- max(imap)
        v2 <- var
        v2[3] <- v2[3]*sqrt(v2[1]*v2[2])  #covariance
        list(theta=NULL, imap=imap, X=X, xmap=imap+ngroup,
             parms=list(v2=v2, theta=var, ngroup=ngroup))
    }
    generate <- function(newtheta, parms) {
        theta <- parms$v2
        varmat <- diag(rep(theta[1:2], each=parms$ngroup))
        for (i in 1:4) varmat[i,i+ngroup] <- varmat[i+ngroup,i] <- theta[3]
        varmat
    }
    
    wrapup <- function(newtheta, b, parms) {
        theta <- parms$theta
        rtheta <- list(site=c(var1=theta[1], var2=theta[2], cor=theta[3]))
        b <- matrix(b, ncol=2)
        dimnames(b) <- list(NULL, c("Intercept", "Slope"))
        list(theta=rtheta, b=list(site=b))
    }
    out <- list(initialize=initialize, generate=generate, wrapup=wrapup)
    class(out) <- 'coxmevar'
    out
    }


