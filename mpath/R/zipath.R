zipath <- function(formula, data, weights, subset, na.action, offset, standardize=TRUE, 
                   family = c("poisson", "negbin", "geometric"),
                   link = c("logit", "probit", "cloglog", "cauchit", "log"),
                   penalty = c("enet", "mnet", "snet"),
                   start = NULL, model = TRUE, y = TRUE, x = FALSE, nlambda=100, lambda.count=NULL, lambda.zero=NULL, penalty.factor.count=NULL, penalty.factor.zero=NULL, lambda.count.min.ratio=.0001, lambda.zero.min.ratio=.1, alpha.count=1, alpha.zero=alpha.count, gamma.count=3, gamma.zero=gamma.count, rescale=FALSE, init.theta, theta.fixed=FALSE, EM=TRUE, maxit.em=200, convtype=c("count", "both"), maxit= 1000, maxit.theta =1, reltol = 1e-5, eps.bino=1e-5, shortlist=FALSE, trace=FALSE, ...)
{
    if(!attr(terms(formula), "intercept"))
        stop("non-intercept model is not implemented")
    if(length(gamma.count) > 1 || length(gamma.zero) > 1)
        stop("gamma.count or gamma.zero must be a scalar")
    family <- match.arg(family)
    penalty <- match.arg(penalty)
    link <- match.arg(link)
    lambda.min.ratio <- lambda.count.min.ratio
    ## set up likelihood
    ziPoisson <- function(parms) {
        ## count mean
        mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
        ## binary mean
        phi <- as.vector(linkinv(Z %*% parms[(kx+1):(kx+kz)] + offsetz))
        
        ## log-likelihood for y = 0 and y >= 1
        loglik0 <- log( phi + exp( log(1-phi) - mu ) ) ## -mu = dpois(0, lambda = mu, log = TRUE)
        loglik1 <- log(1-phi) + dpois(Y, lambda = mu, log = TRUE)
        
        ## collect and return
        loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
        loglik
    }

    ziNegBin <- function(parms) {
        ## count mean
        mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
        ## binary mean
        phi <- as.vector(linkinv(Z %*% parms[(kx+1):(kx+kz)] + offsetz))
        ## negbin size
        theta <- exp(parms[(kx+kz)+1])
        ## log-likelihood for y = 0 and y >= 1
        loglik0 <- log( phi + exp( log(1-phi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE)) ) )
        loglik1 <- log(1-phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

        ## collect and return
        loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
        loglik
    }

                                        #log negative binomial density
    dnbinom1 <- function(x, size, mu){
        n <- size
        p <- size/(size+mu)
        lgamma(x+n) - lgamma(n) - lgamma(x+1) + n*log(p) + x*log(1-p)   
    }        
    ziNegBin1 <- function(parms) {
        ## count mean
        mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
        ## binary mean
        phi <- as.vector(linkinv(Z %*% parms[(kx+1):(kx+kz)] + offsetz))
        ## negbin size
        theta <- exp(parms[(kx+kz)+1])
        ## log-likelihood for y = 0 and y >= 1
        loglik0 <- log( phi + exp( log(1-phi) + dnbinom1(0, size = theta, mu = mu))  )
        loglik1 <- log(1-phi) + dnbinom1(Y, size = theta, mu = mu)
        ## collect and return
        loglik <- rep(NA, n)
        for(i in 1:n){
            if(Y[i]==0) loglik[i] <- weights[i]*loglik0[i]
            else 
                loglik[i] <- weights[i] * loglik1[i]
        }
        loglik
    }
    ziGeom1 <- function(parms) ziNegBin1(c(parms, 0))
    
    ziGeom <- function(parms) ziNegBin(c(parms, 0))
    ziNegBin2 <- function(parms) ziNegBin(c(parms, log(init.theta)))
    pen_zipath <- function(start){
        countcoef1 <- start$count[-1]
        zerocoef1 <- start$zero[-1]
        lc <- length(countcoef1) ### excluding intercept
        if(lc > 0){
            pen1 <- rep(0, lc)
            for(j in 1:lc){
                pen1[j] <- pen_eval(countcoef1[j], lone=lambda.count[k]*alpha.count, ltwo=lambda.count[k]*(1-alpha.count), gamma=gamma.count, penalty=penalty)
            }
        }
        else pen1 <- 0 
        lz <- length(zerocoef1)  ### excluding intercept
        if(lz > 0){
            pen2 <- rep(0, lz)
            for(j in 1:lz)
                pen2[j] <- pen_eval(zerocoef1[j], lone=lambda.zero[k]*alpha.zero, ltwo=lambda.zero[k]*(1-alpha.zero), gamma=gamma.zero, penalty=penalty)
        }
        else pen2 <- 0
        return(n*(sum(pen1) + sum(pen2)))
    }
    pentype <- switch(penalty,
                      "enet"=1,
                      "mnet"=2,
                      "snet"=3)
    pen_zipath1 <- function(start){
        m=length(start$count[-1])
        pencount <- 0
        if(m > 0)
            pencount <- .Fortran("penGLM",
                                 start=as.double(start$count[-1]),
                                 n=as.integer(n),
                                 m=as.integer(length(start$count[-1])),
                                 lambda=as.double(rep(lambda.count[k], m)),
                                 alpha=as.double(alpha.count),
                                 gam=as.double(gamma.count),
                                 penalty=as.integer(pentype),
                                 pen=as.double(0),
                                 package="mpath")$pen
        m=length(start$zero[-1])
        penzero <- 0
        if(m > 0)
            penzero <- .Fortran("penGLM",
                                start=as.double(start$zero[-1]),
                                n=as.integer(n),
                                m=as.integer(length(start$zero[-1])),
                                lambda=as.double(rep(lambda.zero[k], m)),
                                alpha=as.double(alpha.zero),
                                gam=as.double(gamma.zero),
                                penalty=as.integer(pentype),
                                pen=as.double(0),
                                package="mpath")$pen
        return(n*(pencount + penzero))
    } 

    gradPoisson <- function(parms, total=TRUE) {
        ## count mean
        eta <- as.vector(X %*% parms[1:kx] + offsetx)
        mu <- exp(eta)
        ## binary mean
        etaz <- as.vector(Z %*% parms[(kx+1):(kx+kz)] + offsetz)
        muz <- linkinv(etaz)
        
        ## densities at 0
        clogdens0 <- -mu
        dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)

        ## working residuals  
        wres_count <- ifelse(Y1, Y - mu, -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(mu)))
        wres_zero <- ifelse(Y1, -1/(1-muz) * linkobj$mu.eta(etaz),
                            (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
        
        if(total)
            colSums(cbind(wres_count * weights * X, wres_zero * weights * Z))
        else cbind(wres_count * weights * X, wres_zero * weights * Z)
    }
    
    gradGeom <- function(parms, total=TRUE) {
        ## count mean
        eta <- as.vector(X %*% parms[1:kx] + offsetx)
        mu <- exp(eta)
        ## binary mean
        etaz <- as.vector(Z %*% parms[(kx+1):(kx+kz)] + offsetz)
        muz <- linkinv(etaz)

        ## densities at 0
        clogdens0 <- dnbinom(0, size = 1, mu = mu, log = TRUE)
        dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)

        ## working residuals  
        wres_count <- ifelse(Y1, Y - mu * (Y + 1)/(mu + 1), -exp(-log(dens0) +
                                                                 log(1 - muz) + clogdens0 - log(mu + 1) + log(mu)))
        wres_zero <- ifelse(Y1, -1/(1-muz) * linkobj$mu.eta(etaz),
                            (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
        
        if(total)
            colSums(cbind(wres_count * weights * X, wres_zero * weights * Z))
        else cbind(wres_count * weights * X, wres_zero * weights * Z)
    }

    gradNegBin <- function(parms, total=TRUE) {
        ## count mean
        eta <- as.vector(X %*% parms[1:kx] + offsetx)
        mu <- exp(eta)
        ## binary mean
        etaz <- as.vector(Z %*% parms[(kx+1):(kx+kz)] + offsetz)
        muz <- linkinv(etaz)
        ## negbin size
        theta <- exp(parms[(kx+kz)+1])

        ## densities at 0
        clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
        dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)

        ## working residuals  
        wres_count <- ifelse(Y1, Y - mu * (Y + theta)/(mu + theta), -exp(-log(dens0) +
                                                                         log(1 - muz) + clogdens0 + log(theta) - log(mu + theta) + log(mu)))
        wres_zero <- ifelse(Y1, -1/(1-muz) * linkobj$mu.eta(etaz),
                            (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
        wres_theta <- theta * ifelse(Y1, digamma(Y + theta) - digamma(theta) +
                                     log(theta) - log(mu + theta) + 1 - (Y + theta)/(mu + theta),
                                     exp(-log(dens0) + log(1 - muz) + clogdens0) *
                                     (log(theta) - log(mu + theta) + 1 - theta/(mu + theta)))

        if(total)
            colSums(cbind(wres_count * weights * X, wres_zero * weights * Z, wres_theta))
        else cbind(wres_count * weights * X, wres_zero * weights * Z, wres_theta)
    }
    checkconv <- function(start, start_old, type=c("count","both"), eps=.Machine$double.eps, thresh=.001){
        type <- match.arg(type)
        checkcount <- .Fortran("checkConvergence",
                               m = as.integer(length(start$count)),
                               beta = as.double(start$count),
                               beta_old = as.double(start_old$count),
                               eps = as.double(eps),
                               thresh = as.double(thresh),
                               converged = as.integer(0),
                               package = "mpath")$converged
        if(type=="count") return(checkcount)
        else
            if(checkcount==0) return(FALSE)
            else {
                checkzero <- .Fortran("checkConvergence",
                                      m = as.integer(length(start$zero)),
                                      beta = as.double(start$zero),
                                      beta_old = as.double(start_old$zero),
                                      eps = as.double(eps),
                                      thresh = as.double(thresh),
                                      converged = as.integer(0),
                                      package = "mpath")$converged
                return(checkzero)
            }
    }
    family <- match.arg(family)
    convtype <- match.arg(convtype)
    penalty <- match.arg(penalty)
    loglikfun <- switch(family,
                        "poisson" = ziPoisson,
                        "geometric" = ziGeom,
                        "negbin" = ziNegBin)
    gradfun <- switch(family,
                      "poisson" = gradPoisson,
                      "geometric" = gradGeom,
                      "negbin" = gradNegBin)
    ## binary link processing
    linkstr <- match.arg(link)
    linkobj <- make.link(linkstr)
    linkinv <- linkobj$linkinv

    if(trace) cat("Zero-inflated Count Model\n",
                  paste("count model:", family, "with log link\n"),
                  paste("zero-inflation model: binomial with", linkstr, "link\n"), sep = "")
    
    
    ## call and formula
    cl <- match.call()
    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    
    ## extended formula processing
    if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
        {
            ff <- formula
            formula[[3]][1] <- call("+")
            mf$formula <- formula
            ffc <- . ~ .
            ffz <- ~ .
            ffc[[2]] <- ff[[2]]
            ffc[[3]] <- ff[[3]][[2]]
            ffz[[3]] <- ff[[3]][[3]]
            ffz[[2]] <- NULL
        } else {
            ffz <- ffc <- ff <- formula
            ffz[[2]] <- NULL
        }
    if(inherits(try(terms(ffz), silent = TRUE), "try-error")) {
        ffz <- eval(parse(text = sprintf( paste("%s -", deparse(ffc[[2]])), deparse(ffz) )))
    }

    ## call model.frame()
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    ## extract terms, model matrices, response
    mt <- attr(mf, "terms")
    mtX <- terms(ffc, data = data)
    X <- model.matrix(mtX, mf)
    Xold <- X
### This means X possibly has intercept (1st) column depending on the specifications 
    Xnew <- as.matrix(X[,-1])
    mtZ <- terms(ffz, data = data)
    mtZ <- terms(update(mtZ, ~ .), data = data)
    Z <- model.matrix(mtZ, mf)
    Zold <- Z
    Znew <- as.matrix(Z[,-1])
    Y <- model.response(mf, "numeric")

    ## sanity checks
    if(length(Y) < 1) stop("empty model")
    if(all(Y > 0)) stop("invalid dependent variable, minimum count is not zero")  
    if(!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
        stop("invalid dependent variable, non-integer values")
    Y <- as.integer(round(Y + 0.001))
    if(any(Y < 0)) stop("invalid dependent variable, negative counts")

    if(trace) {
        cat("dependent variable:\n")
        tab <- table(factor(Y, levels = 0:max(Y)), exclude = NULL)
        names(dimnames(tab)) <- NULL
        print(tab)
    }
    
    ## convenience variables
    n <- length(Y)
    kx <- NCOL(X)
    kz <- NCOL(Z)
    Y0 <- Y <= 0
    Y1 <- Y > 0

    ## weights and offset
    weights <- model.weights(mf)
    if(is.null(weights)) weights <- 1
    if(length(weights) == 1) weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)
### Todo: offset  
    offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
    ##  offsetx <- NULL
    if(is.null(offsetx)) offsetx <- 0
    if(length(offsetx) == 1) offsetx <- rep.int(offsetx, n)
    offsetx <- as.vector(offsetx)
### Todo: offset
    offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
    if(is.null(offsetz)) offsetz <- 0
    if(length(offsetz) == 1) offsetz <- rep.int(offsetz, n)
    offsetz <- as.vector(offsetz)
### standardize variables
    if(standardize){
        wt <- weights/sum(weights)
        one <- rep(1, length(Y))
        if(dim(Xnew)[2] > 0){
            xx <- Xnew
            meanx <- drop(wt %*% Xnew)
            Xnew <- scale(Xnew, meanx, FALSE) # centers x
            Xnew <- sqrt(wt) * Xnew
            normx <- sqrt(drop(one %*% (Xnew^2)))
            X[,-1] <- Xnew <- scale(xx, meanx, normx)
        }
        if(dim(Znew)[2] > 0){
            zz <- Znew
            meanz <- drop(wt %*% Znew)
            Znew <- scale(Znew, meanz, FALSE) # centers x
            Znew <- sqrt(wt) * Znew
            normz <- sqrt(drop(one %*% (Znew^2)))
            Z[,-1] <- Znew <- scale(zz, meanz, normz)
        }
    }

    ## starting values
    if(!is.null(start)) {
        valid <- TRUE
        if(!("count" %in% names(start))) {
            valid <- FALSE
            warning("invalid starting values, count model coefficients not specified")
            start$count <- rep.int(0, kx)
        }
        if(!("zero" %in% names(start))) {
            valid <- FALSE
            warning("invalid starting values, zero-inflation model coefficients not specified")
            start$zero <- rep.int(0, kz)
        }
        if(length(start$count) != kx) {
            valid <- FALSE
            warning("invalid starting values, wrong number of count model coefficients")
        }
        if(length(start$zero) != kz) {
            valid <- FALSE
            warning("invalid starting values, wrong number of zero-inflation model coefficients")
        }
        if(family == "negbin") {
            if(!("theta" %in% names(start))) start$theta <- 1
            start <- list(count = start$count, zero = start$zero, theta = as.vector(start$theta[1]))
        } else {
            start <- list(count = start$count, zero = start$zero)
        }
        if(!valid) start <- NULL
    }
    if(!is.null(lambda.count) && !is.null(lambda.zero)){
        if(length(lambda.count) != length(lambda.zero))
            stop("length of lambda.count must be the same as lambda.zero\n")
        else nlambda <- length(lambda.count)
    } 
    if(family == "negbin") theta <- rep(NA, nlambda)
    coefc <- matrix(NA, nrow=dim(X)[2], ncol=nlambda) 
    coefz <- matrix(NA, nrow=dim(Z)[2], ncol=nlambda) 
    se <- matrix(NA, nrow=dim(X)[2]+dim(Z)[2]+(family=="negbin"), ncol=nlambda) 
    ll <- rep(NA, nlambda)
    convout <- rep(NA, nlambda)
    fit0 <- zeroinfl(Y~1|1, weights=weights, dist=family)
### find solution of intercept model with family = negbin and fixed theta parameter, and replace the estimates for fit0
    ## ML estimation of theta
    if(family=="negbin" && theta.fixed) {
        if(is.null(init.theta)){
### this may not work when a model is complex, i.e., solve(fit) may not be stable in zeroinfl...
### ff is the original "formula"
            init.theta <- try(do.call("zeroinfl", list(ff, data, weights=weights, dist="negbin"))$theta)
            if(inherits(init.theta, "try-error"))
                init.theta <- do.call("glm.nb", list(formula=Y ~ X[,-1], data, weights=weights))$theta
        }
        Xtmp <- X
        Ztmp <- Z
        X <- as.matrix(X[,1])
        Z <- as.matrix(Z[,1])
        kx <- 1
        kz <- 1 
### find the solution of estimates for the intercept-only model with fixed init.theta 
        fit0.nb <- optim(fn = ziNegBin2, gr = gradfun,
                         par = c(fit0$coef$count, fit0$coef$zero), #if(dist == "negbin") log(start$theta) else NULL),
                         method = "Nelder-Mead", hessian = FALSE, control = list(fnscale=-1))
        fit0$coefficients <- list(count=fit0.nb$par[1], zero=fit0.nb$par[2])
        fit0$theta <- init.theta
        X <- Xtmp
        Z <- Ztmp
        kx <- NCOL(X)
        kz <- NCOL(Z)
    }
    if(is.null(penalty.factor.count)) penalty.factor.count <- rep(1, dim(Xnew)[2])
    if(is.null(penalty.factor.zero)) penalty.factor.zero <- rep(1, dim(Znew)[2])
### get lambda -begin
    if(is.null(lambda.count) && is.null(lambda.zero)){
        mui <- predict(fit0, type="count")
        probi <- predict(fit0, type="zero")
        probi <- switch(family,
                        "poisson" = probi/(probi + (1-probi) * dpois(0, mui)),
                        "geometric" = probi/(probi + (1-probi) * dnbinom(0, size = 1, mu = mui)),
                        "negbin"= probi/(probi + (1-probi) * dnbinom(0, size = fit0$theta, mu = mui))
                        )
        probi[Y1] <- 0
        if(dim(Xnew)[2] > 1){
            if(family=="poisson"){
                lambda.count <- glmreg_fit(Xnew, Y, weights=weights*(1-probi)/n, family="poisson", penalty.factor=penalty.factor.count, nlambda=nlambda, lambda.min.ratio = lambda.min.ratio, alpha=alpha.count, maxit=1, standardize=FALSE, ...)$lambda
            }
            else if(family=="geometric"){
                lambda.count <- glmreg_fit(Xnew, Y, weights=weights*(1-probi)/n, family="negbin", penalty.factor=penalty.factor.count, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio,alpha=alpha.count, theta=1, maxit=1, standardize=FALSE, ...)$lambda
                theta <- 1
            }   else if(family=="negbin"){
                if(theta.fixed){
                    lambda.count <- glmreg_fit(Xnew, Y, weights=weights*(1-probi)/n, family="negbin", penalty.factor=penalty.factor.count, theta=init.theta, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, alpha=alpha.count, maxit=1, standardize=FALSE, ...)$lambda
                }
                else{
                    lambda.count <- glmregNB(Y~Xnew, weights=weights*(1-probi)/n, penalty.factor=penalty.factor.count, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, alpha=alpha.count, rescale=FALSE, standardize=FALSE, maxit=1, maxit.theta=2, ...)$lambda
                }
            }
        }
        else lambda.count <- rep(0, nlambda)
        if(dim(Znew)[2] > 1){
            lambda.zero <- glmreg_fit(Znew, probi, weights=weights, family="binomial", penalty=penalty, penalty.factor=penalty.factor.zero, rescale=rescale, nlambda=nlambda, lambda.min.ratio=lambda.zero.min.ratio, alpha=alpha.zero, maxit=1, eps.bino=eps.bino, ...)$lambda
       # if(penalty %in% c("mnet", "snet") && !rescale) 
       #   {###it seems lambda.zero[1] too small, thus ...
       #    lmax <- 0.632 * sqrt(lambda.zero[1])
       #    lambda.zero <- seq(log(lmax), log(lambda.zero.min.ratio * lmax), length.out=nlambda)
       #    lambda.zero <- exp(lambda.zero)
       #   }
        }    else lambda.zero <- rep(0, nlambda)
    }
    else if(length(lambda.count) != length(lambda.zero))
        stop("length of lambda.count must be the same as lambda.zero\n")
    else nlambda <- length(lambda.count) 
### get lambda -end
    model_count <- list(coefficients = c(fit0$coefficients$count, rep(0, dim(Xnew)[2])), fitted.value=predict(fit0, type="count"))
    model_zero <- list(coefficients = c(fit0$coefficients$zero, rep(0, dim(Znew)[2])), fitted.value=predict(fit0, type="zero"))
    start <- list(count = model_count$coefficients, zero = model_zero$coefficients)
    vcov <- vector("list", length=nlambda)
    if(family == "negbin"){
        if(!theta.fixed)
            start$theta <- fit0$theta
        else start$theta <- init.theta
    }
    k <- 1
    while(k <= nlambda){
        if(trace) {cat("\nOuter loop: Decremental lambda\n")
                   cat(" lambda iteration", k, "\n")
               }
        if(EM & family == "poisson") {
            mui <- model_count$fitted.value
            probi <- model_zero$fitted.value
            probi <- probi/(probi + (1-probi) * dpois(0, mui))
            probi[Y1] <- 0
            if(trace && penalty=="enet"){
### initial penalized weighted NB log-likelihood value (in the EM complete loglik) for penalty="enet" with alpha=1
                eta <- model_zero$coef[1] + as.vector(Znew %*% coef(model_zero)[-1])
                mu = 1.0/(1+exp(-eta))
                w <- weights
                ll2init <- sum(w*(probi*eta - log(1+exp(eta)))) - length(Y)*lambda.zero[k]*sum(abs(coef(model_zero)[-1]))
            }
            ll_raw <- loglikfun(c(start$count, start$zero))
            ll_new <- ll_raw - pen_zipath1(start)
            j <- 1
            testH1 <- testH <- rep(NA, maxit.em)
            converged <- FALSE
            while(j <= maxit.em && !converged) {
                if(trace) cat("EM iteration ", j, "\n")
                ll_old <- ll_new
                start_old <- start
                                        #        start_old <- start
                if(dim(X)[2]==1){ ###intercept-only model, Todo: need to check
                    model_count <- glm.fit(X, Y, weights = weights * (1-probi), offset = offsetx,
                                           family = poisson(), start = start$count)
                    countcoef = model_count$coefficients
                }
                else{
                    capture.output(model_count <- glmreg_fit(Xnew, Y, weights = weights * (1-probi)/n, standardize = FALSE, penalty.factor=penalty.factor.count, lambda=lambda.count[k], alpha=alpha.count,gamma=gamma.count, rescale=rescale, maxit=maxit,
                                        #                                         start=start$count, mustart=mui, etastart=g(mui, family="poisson"),  
                                                             family = "poisson", penalty=penalty, trace=trace))
                    countcoef = c(model_count$b0, model_count$beta)
                    if(trace && penalty=="enet"){
                        eta <- model_count$b0 + as.vector(Xnew %*% model_count$beta)
                        mu <- exp(eta)
                        w <- weights*(1-probi)/n
                        ll1 <- sum(w*(-mu + Y*log(mu)-lgamma(Y+1))) - lambda.count[k]*sum(abs(model_count$beta))
                        cat("penalized weighted Poisson loglik by updated parameter ll1", ll1, "\n")
                        tmp <- which(abs(model_count$pll) > 0)
                        cat("penalized loglik value pll", model_count$pll[tmp], "\n")
                    }
                }
                if(dim(Z)[2]==1){
                    model_zero <- suppressWarnings(glm.fit(Z, probi, weights = weights, offset = offsetz,
                                                           family = binomial(link = linkstr), start = start$zero))
                    zerocoef = model_zero$coefficients
                }
                else{
                    capture.output(model_zero <- glmreg_fit(Znew, probi, weights = weights, standardize=standardize, penalty.factor=penalty.factor.zero, lambda=lambda.zero[k], alpha=alpha.zero, gamma=gamma.zero, rescale=rescale, maxit=maxit, 
                                        #                                        start=start$zero, mustart=model_zero$fitted.value, etastart=g(model_zero$fitted.value, family="binomial", eps.bino=eps.bino),
                                                            family = "binomial", penalty=penalty, trace=trace, eps.bino=eps.bino))
                    zerocoef = c(model_zero$b0, model_zero$beta)
                    if(trace && penalty=="enet"){
                        eta <- model_zero$b0 + as.vector(Znew %*% model_zero$beta)
                        mu = 1.0/(1+exp(-eta))
                        w <- weights/sum(weights)
                        ll2 <- sum(w*(probi*eta - log(1+exp(eta)))) - lambda.zero[k]*sum(abs(model_zero$beta))
                        cat("penalized binomial loglik by updated parameter ll2", ll2, "\n")
                        cat("model_zero$pll", model_zero$pll, "\n")
                        cat("model_zero$pllres", model_zero$pllres, "\n")
                    }
                }
                start <- list(count=countcoef, zero=zerocoef)
                ll_raw <- loglikfun(c(start$count, start$zero))
                ll_new <- ll_raw - pen_zipath1(start)
                if(rescale && penalty!="enet")
                    converged <- checkconv(start, start_old, type=convtype, thresh=reltol)
                else converged <- (ll_old - ll_new)/ll_old <= reltol
                probi_old <- probi
                mui <- model_count$fitted.value
                probi <- model_zero$fitted.value
                probi <- probi/(probi + (1-probi) * dpois(0, mui))
                probi[Y1] <- 0
                if(trace){ 
### compute H, expectation of log(f(z|gamma, beta, y)), cf Lambert(1992), page 4
### compute expected conditional loglik H, cf: Theory of Point Estimation, 2nd edition, (4.24) page 459
                    if(penalty=="enet"){
                        testH1[j] <- sum(log(probi[Y0])*probi_old[Y0])
                        testH[j] <- sum(log(probi[Y0]/probi_old[Y0])*probi_old[Y0])
                        H <- sum(log(probi[Y0])*probi_old[Y0]+log(1-probi[Y0])*(1-probi_old[Y0]))
                        llcc <- length(Y)*(ll1+ll2) - H
                        cat("H=", H, "\n")
                        cat("penalized incomplete ll1+ll2-H", llcc, "\n")
                    }
                    cat("penalized incomplete data ll_new=", ll_new, "\n")
                    cat("ll_new-ll_old", ll_new-ll_old, "\n")
                }
                j <- j + 1
            }
        }

        if(EM & family == "geometric") {
            mui <- model_count$fitted.value
            probi <- model_zero$fitted.value
            probi <- probi/(probi + (1-probi) * dnbinom(0, size = 1, mu = mui))
            probi[Y1] <- 0
            ll_raw <- loglikfun(c(start$count, start$zero))
            ll_new <- ll_raw - pen_zipath1(start)
            ll_old <- ll_new
            
            j <- 1
            converged <- FALSE
            while(j <= maxit.em && !converged) {
                ll_old <- ll_new
                start_old <- start
                                        #        start_old <- start
                if(dim(X)[2]==1){ ###intercept-only model, Todo: need to check
                    model_count <- suppressWarnings(glm.fit(X, Y, weights = weights * (1-probi),
                                                            offset = offsetx, family = negative.binomial(1), start = start$count))
                    countcoef = model_count$coefficients
                }
                else{
                    model_count <- glmreg_fit(Xnew, Y, weights = weights * (1-probi)/n, standardize=FALSE, penalty.factor=penalty.factor.count, lambda=lambda.count[k], alpha=alpha.count, gamma=gamma.count,rescale=rescale, maxit=maxit, #start=start$count, mustart=mui, etastart=g(mui,family="negbin"),
                                              family = "negbin", theta=1, penalty=penalty)
                    countcoef = c(model_count$b0, model_count$beta)
                }
                if(dim(Z)[2]==1){ ###intercept-only model, Todo: need to check
                    model_zero <- suppressWarnings(glm.fit(Z, probi, weights = weights,
                                                           offset = offsetz, family = binomial(link = linkstr), start = start$zero))
                    zerocoef = model_zero$coefficients
                }
                else{
                    model_zero <- glmreg_fit(Znew, probi, weights = weights, standardize=standardize, penalty.factor=penalty.factor.zero, lambda=lambda.zero[k], alpha=alpha.zero, gamma=gamma.zero, rescale=rescale, maxit=maxit,
                                        #start=start$zero, mustart=model_zero$fitted.value, etastart=g(model_zero$fitted.value, family="binomial", eps.bino=eps.bino),
                                             family = "binomial", penalty=penalty, eps.bino=eps.bino)
                    zerocoef = c(model_zero$b0, model_zero$beta)
                }
                start <- list(count=countcoef, zero=zerocoef)
                                        #        start <- list(count = model_count$coefficients, zero = model_zero$coefficients)
                ll_raw <- loglikfun(c(start$count, start$zero))
                ll_new <- ll_raw - pen_zipath1(start)
                if(rescale && penalty!="enet")
                    converged <- checkconv(start, start_old, type=convtype, thresh=reltol)
                else converged <- abs((ll_old - ll_new)/ll_old) <= reltol
                mui <- model_count$fitted.value
                probi <- model_zero$fitted.value
                probi <- probi/(probi + (1-probi) * dnbinom(0, size = 1, mu = mui))
                probi[Y1] <- 0
                if(trace) cat("j=", j, "ll_new=", ll_new,"ll_new-ll_old=", ll_new-ll_old, "\n")
                j <- j + 1
            }
        }

        if(EM & family == "negbin") {
            mui <- model_count$fitted.value
            probi <- model_zero$fitted.value
### initial penalized logistic log-likelihood value (in the EM complete loglik) for penalty="enet" with alpha=1
            probi <- probi/(probi + (1-probi) * dnbinom(0, size = start$theta, mu = mui))
            probi[Y1] <- 0
            if(trace && penalty=="enet"){
### initial penalized weighted NB log-likelihood value (in the EM complete loglik) for penalty="enet" with alpha=1
                loglik1 <- (1-probi)/n*suppressWarnings(dnbinom(Y, size = start$theta, mu = mui, log = TRUE))
                loglik1 <- sum(loglik1) - lambda.count[k]*sum(abs(start$count[-1]))
		eta <- coef(model_zero)[1] + as.vector(Znew %*% coef(model_zero)[-1])
                w <- weights
                ll2init <- sum(w*(probi*eta - log(1+exp(eta)))) - length(Y)*lambda.zero[k]*sum(abs(coef(model_zero)[-1]))
                cat("loglik1=", loglik1, "ll2init=", ll2init, "\n")
            }
            ll_raw <- loglikfun(c(start$count, start$zero, log(start$theta)))
            ll_new <- ll_raw - pen_zipath1(start)
### intercept is not penalized      
            ll_old <- ll_new
            
            ## offset handling in glm.nb is sub-optimal, hence...
            offset <- offsetx
            j <- 1
            converged <- FALSE
            while(j <= maxit.em && !converged) {
               ll_old <- ll_new
                if(trace) cat("\nEM iteration", j, "ll_old", ll_old, "\n")
                start_old <- start
                if(dim(X)[2]==1){
                    if(theta.fixed){
                        model_count <- suppressWarnings(glm.fit(X, Y, weights = weights * (1-probi),
                                                                offset = offsetx, family = negative.binomial(init.theta), start = start$count))
                        model_count$theta <- start$theta
                    }
                    else model_count <- suppressWarnings(glm.nb(Y ~ 0 + X + offset(offset), weights = weights * (1-probi),
                                                                start = start$count, init.theta = start$theta))
                    countcoef = model_count$coefficients
                }        else{ 
                    if(theta.fixed)
                        capture.output(model_count <- glmreg_fit(Xnew, Y, weights = weights * (1-probi)/n, standardize=FALSE, penalty.factor=penalty.factor.count, lambda=lambda.count[k], alpha=alpha.count, gamma=gamma.count,rescale=rescale, maxit=maxit,
                                        #                      start=start$count, mustart=mui, etastart=g(mui,family="negbin"),
                                                                 family = "negbin", theta=init.theta, penalty=penalty, trace=trace, ...))
                    else capture.output(model_count <- glmregNB(Y ~ X[,-1], standardize=FALSE, weights = weights * (1-probi)/n, penalty.factor=penalty.factor.count, lambda=lambda.count[k], alpha=alpha.count, gamma=gamma.count, rescale=rescale, maxit.theta = maxit.theta, maxit=maxit,
### not inclued line otherwise not accurate###
                                        #   start=start$count, mustart=mui, etastart=g(mui,family="negbin"),
                                                                init.theta=start$theta, penalty=penalty, trace=trace, ...))
                    countcoef = c(model_count$b0, model_count$beta)
                }
                if(dim(Z)[2]==1){
                    model_zero <- suppressWarnings(glm.fit(Z, probi, weights = weights, offset = offsetz,
                                                                   family = binomial(link = linkstr), start = start$zero))
                    zerocoef = model_zero$coefficients
                }
                else{
                    capture.output(model_zero <- glmreg_fit(Znew, probi, weights = weights, standardize=TRUE, penalty.factor=penalty.factor.zero, lambda=lambda.zero[k], alpha=alpha.zero, gamma=gamma.zero, rescale=rescale, maxit=maxit,
### not inclued line otherwise not accurate###
                                        #start=start$zero, mustart=model_zero$fitted.value, etastart=g(model_zero$fitted.value, family="binomial", eps.bino=eps.bino),
                                              family = "binomial", penalty=penalty, trace=trace, eps.bino=eps.bino, ...))
                                        #if(nlambda < 100){
                                        #       	cat("k=",k, "model_zero$satu=", model_zero$satu, "\n")
                                        #}
                    if(model_zero$satu==1){
                        if(k==1){
                            model_zero <- suppressWarnings(glm.fit(Z[,1], probi, weights = weights, offset = offsetz,
                                                                   family = binomial(link = linkstr)))
                            zerocoef = model_zero$coefficients
                            model_zero$satu <- 1
                        warnings("zero component model maybe too complex\n")
                        }
                        break
                    }
	            zerocoef = c(model_zero$b0, model_zero$beta)
                }
                if(trace && penalty=="enet"){
                    eta <- model_zero$b0 + as.vector(Znew %*% model_zero$beta)
                    w <- weights
                    ll2 <- sum(w*(probi*eta - log(1+exp(eta)))) - length(Y)*lambda.zero[k]*sum(abs(model_zero$beta))
                }
                start <- list(count=countcoef, zero=zerocoef, theta=model_count$theta)
                ll_raw <- loglikfun(c(start$count, start$zero, log(start$theta)))
                ll_new <- ll_raw - pen_zipath1(start)
                if(rescale && penalty!="enet")
                    converged <- checkconv(start, start_old, thresh=reltol)
                else converged <- abs((ll_old - ll_new)/ll_old) <= reltol
                probi_old <- probi
                mui <- model_count$fitted.values
                probi <- model_zero$fitted.values
                probi <- probi/(probi + (1-probi) * dnbinom(0, size = start$theta, mu = mui))
                probi[Y1] <- 0
                if(trace){
                    cat("EM iteration", j, "\n")
                    if(penalty=="enet"){
###note: model_count$pllres is for weights (1-probi)/n, but we need (1-probi) instead
                        ll1 <- logLik(model_count)*n - length(Y)*lambda.count[k]*sum(abs(model_count$beta))
                        ll2 <- model_zero$pllres[length(model_zero$pllres)]
                        cat("penalized weighted NB loglik by updated parameter ll1", ll1, "\n")
                        cat("penalized logistic regresion loglik by updated parameter ll2", ll2, "\n")
                        loglik1new <- (1-probi_old)*suppressWarnings(dnbinom(Y, size = start$theta, mu = mui, log = TRUE))
                        loglik1new <- sum(loglik1new) - length(Y)*lambda.count[k]*sum(abs(model_count$beta))
                        cat("penalized weighted NB loglik by updated parameter for weights used in the model loglik1new", loglik1new, "\n")
                        cat("penalized weighted NB loglik by updated parameter for weights used in the model", logLik(model_count) - length(Y)*lambda.count[k]*sum(abs(model_count$beta)), "\n")
### initial penalized logistic log-likelihood value (in the EM complete loglik) for penalty="enet" with alpha=1
                        H <- sum(log(probi[Y0])*probi_old[Y0]+log(1-probi[Y0])*(1-probi_old[Y0]))
                        cat("H with updated parameter", H, "\n")
                        llcc <- ll1+ll2 - H
                        cat("model_count$pll", model_count$pll, "\n")
                        cat("model_zero$pll", model_zero$pll, "\n")
                        cat("penalized incomplete ll1+ll2-H", llcc, "\n")
                        cat("ll_old=", ll_old, "\n")
                        cat("ll_new=", ll_new,"ll_new-ll_old=", ll_new-ll_old, "\n")
                    }
                }
                j <- j + 1
            }
        }
        if(is.null(model_zero$satu)) model_zero$satu <- 0
        if(model_zero$satu==1 && k>1) break
	ll[k] <- ll_raw
        convout[k] <- converged
        coefc[,k] <- start$count
        coefz[,k] <- start$zero
        if(family == "negbin") 
            theta[k] <- model_count$theta
### standard errors, experiment code, compare to sandwichreg function 
        need.se <- FALSE
        if(need.se){    
            if(family!="negbin")
                ball <- c(start$count, start$zero)
            else ball <- c(start$count, start$zero, log(start$theta))
### test gradient
            eta <- as.vector(X %*% start$count)
            mu <- exp(eta)
            ## binary mean
            etaz <- as.vector(Z %*% start$zero)
            muz <- linkinv(etaz)
            ## negbin size
            theta <- start$theta
            h0 <- (theta/(mu+theta))^theta
            h1 <- (muz+(1-muz)*(theta/(mu+theta))^theta)
            wres_count <- ifelse(Y1, Y-(Y+theta)*mu/(mu+theta),
                                 -(1-muz)*theta^theta*(mu+theta)^(-theta)*theta*mu/(mu+theta)/h1)
            wres_zero <- ifelse(Y1, -muz, muz*(1-muz)*(1-h0)/h1)
            wres_theta <- ifelse(Y1, digamma(theta+Y)-digamma(theta) - Y/(mu+theta)+log(theta)+1-log(mu+theta)-theta/(mu+theta),
                                 (1-muz)*h0*(log(theta)+1-log(mu+theta)-theta/(mu+theta))/h1)
            res1 <- colSums(cbind(wres_count*X, wres_zero*Z, wres_theta)) ### this is different from gradfun w.r.t. theta because of parameterization there is log(theta)
###res1 <- (cbind(wres_count*X, wres_zero*Z, wres_theta*theta)) ### this is different from gradfun w.r.t. theta because of parameterization there is log(theta)
            res2 <- gradfun(ball, total=TRUE)
### end of test gradient
            grad <- gradfun(ball, total=FALSE)
            tmp2 <- zeroinfl(formula, data, dist="negbin")
            tmp <- apply(grad, 2, mean)
            covLoglik <- 1/n*crossprod(grad)  - tmp %o% tmp 
### numerical calculation of hessian matrix can be slow from package "numDeriv"
            hes <- hessian(func=loglikfun, x=ball)
### be careful with theta for family="negbin"
            sigma <- rep(0, kx+kz+(family=="negbin"))
            for(ss in 2:kx)
                sigma[ss] <- pen2_eval(abs(start$count[ss]), lone=lambda.count[k]*alpha.count, ltwo=lambda.count[k]*(1-alpha.count), gamma=gamma.count, penalty=penalty)
            for(ss in (kx+2):(kx+kz))
                sigma[ss] <- pen2_eval(abs(start$zero[ss-kx]), lone=lambda.zero[k]*alpha.zero, ltwo=lambda.zero[k]*(1-alpha.zero), gamma=gamma.zero, penalty=penalty)
    	    s1 <- try(-solve(hes - n*diag(sigma))*n)
	    vcov[[k]] <- s1
	    if(!inherits(s1, "try-error"))
                tmp <- 1/n*s1%*%covLoglik%*%s1 
            else next
            tmp <- sqrt(diag(tmp))
### coef.standardized/se.standardized = coef.unstandardized/se.unstandardized except for the intercepts
            if(standardize){
                tmp[2:kx] <- tmp[2:kx]/normx
                tmp[(kx+2):(kx+kz)] <- tmp[(kx+2):(kx+kz)]/normz
            }
            se[,k] <- tmp
        }
### standard errors
        if(model_zero$satu==1 && k>1) break
        k <- k + 1
    }
    if(model_zero$satu==1){
	nlambda <- k-1
	coefc <- coefc[,1:nlambda]
	coefz <- coefz[,1:nlambda]
        if(family == "negbin")
            theta <- theta[1:nlambda]
        lambda.count <- lambda.count[1:nlambda]
        lambda.zero <- lambda.zero[1:nlambda]
        ll <- ll[1:nlambda]
    }
    if(standardize) {
        coefc <- as.matrix(coefc)
        if(dim(coefc)[1] > 1){
            coefc[-1,] <- coefc[-1,]/normx 
                                        #      if(dim(coefc)[1] > 2)
            coefc[1,] <- coefc[1,] - crossprod(meanx, coefc[-1,])
                                        #      else coefc[1,] <- coefc[1,] - meanx * coefc[-1,] ### check
        }
        coefz <- as.matrix(coefz)
        if(dim(coefz)[1] > 1){
            coefz[-1,] <- coefz[-1,]/normz
                                        #      if(dim(coefz)[1] > 2)
            coefz[1,] <- coefz[1,] - crossprod(meanz, coefz[-1,])
                                        #      else coefz[1,] <- coefz[1,] - meanz * coefz[-1,]
        }
    }
    if(family != "negbin")  theta <- NULL
    else names(theta) <- lambda.count
    rownames(coefc) <- names(start$count) <- colnames(X)
    rownames(coefz) <- names(start$zero) <- colnames(Z)
    if(nlambda > 1){
        if(lambda.count[1] > 0.0001)
            colnames(coefc) <- round(lambda.count, digits=4)
        else colnames(coefc) <- round(lambda.count, digits=8)
        if(lambda.zero[1] > 0.0001)
            colnames(coefz) <- round(lambda.zero, digits=4)
        else colnames(coefz) <- round(lambda.zero, digits=8)
    } 
    ## fitted and residuals
    mu <- exp(Xold %*% coefc + offsetx)
    phi <- linkinv(Zold %*% coefz + offsetz)
    Yhat <- (1-phi) * mu
    res <- sqrt(weights) * (Y - Yhat)
    ## effective observations
    nobs <- sum(weights > 0) ## = n - sum(weights == 0)
    dfc <- apply(abs(coefc) > 0, 2, sum) 
    if(family=="negbin") dfc <- 1+dfc
    dfz <- apply(abs(coefz) > 0, 2, sum) 
    aic <- -2*ll + 2*(dfc+dfz)
    bic <- -2*ll + log(n)*(dfc+dfz)
                                        #if(!need.se) covLoglik <- hes <- se <- vcov <- NA
    if(shortlist)
        rval <- list(coefficients = list(count = coefc, zero = coefz),
                     lambda.count = lambda.count,
                     lambda.zero = lambda.zero,
                     alpha.count = alpha.count,
                     alpha.zero = alpha.zero,
                     gamma.count = gamma.count,
                     gamma.zero = gamma.zero,
                     loglik = ll,
                     aic = aic,
                     bic = bic,                        
                     family = family,
                     nlambda = nlambda,
                     theta = theta)
    else
        rval <- list(coefficients = list(count = coefc, zero = coefz),
                     residuals = res,
                     fitted.values = Yhat,
                     satu = model_zero$satu,
                     nlambda=nlambda,
                     lambda.count = lambda.count,
                     lambda.zero = lambda.zero,
                     alpha.count = alpha.count,
                     alpha.zero = alpha.zero,
                     gamma.count = gamma.count,
                     gamma.zero = gamma.zero,
                     start = start,
                     weights = if(identical(as.vector(weights), rep.int(1L, n))) NULL else weights,
                     offset = list(count = if(identical(offsetx, rep.int(0, n))) NULL else offsetx,
                         zero = if(identical(offsetz, rep.int(0, n))) NULL else offsetz),
                     n = nobs,
                     df.null = nobs - 2, 
                     df.residual = dfc + dfz,
                     terms = list(count = mtX, zero = mtZ, full = mt),
                     theta = theta,
                     theta.fixed = theta.fixed,
                     loglik = ll,
                     aic = aic,
                     bic = bic,                         
                     family = family,
                     penalty=penalty,
                     link = linkstr,
                     linkinv = linkinv,
                     converged = convout,
                     call = cl,
                     formula = ff,
                     levels = .getXlevels(mt, mf),
                     contrasts = list(count = attr(X, "contrasts"), zero = attr(Z, "contrasts"))
                     )
    if(model) rval$model <- mf
    if(y) rval$y <- Y
    if(x) rval$x <- list(count = Xold, zero = Zold)
    class(rval) <- "zipath"
    return(rval)
}

zipath.control <- function(method = "BFGS", maxit.em=50, maxit = 10000, trace = FALSE, EM = TRUE, start = NULL, ...) {
    rval <- list(method = method, maxit.em = maxit.em, maxit = maxit, trace = trace, EM = EM, start = start)
    rval <- c(rval, list(...))
    rval
}

coef.zipath <- function(object, which=1:object$nlambda, model = c("full", "count", "zero"), ...) {
    model <- match.arg(model)
    rval <- object$coefficients
    rval <- switch(model,
                   "full" = list(count=rval$count[, which], zero=rval$zero[, which]), 
                   "count" = rval$count[,which],
                   "zero" = rval$zero[,which])
    rval
}

print.zipath <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    
    cat(paste("Count model coefficients (", x$family, " with log link):\n", sep = ""))
    print.default(format(x$coefficients$count, digits = digits), print.gap = 2, quote = FALSE)
    if(x$family == "negbin"){
        cat("\nTheta:", round(x$theta, digits),"\n")} else cat("\n")
    
    cat(paste("\nZero-inflation model coefficients (binomial with ", x$link, " link):\n", sep = ""))
    print.default(format(x$coefficients$zero, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
                                        #  }
    
    invisible(x)
}

summary.zipath <- function(object,...)
{
    ## residuals
    object$residuals <- residuals(object, type = "pearson")
    
    if(object$family == "negbin") {
        object$coefficients$count <- rbind(object$coefficients$count, "Log(theta)" = log(object$theta))
    }
    
    object$fitted.values <- object$terms <- object$model <- object$y <-
        object$x <- object$levels <- object$contrasts <- object$start <- NULL

    class(object) <- "summary.zipath"
    object
}

print.summary.zipath <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    cat(paste("\nCount model coefficients (", x$family, " with log link):\n", sep = ""))
    printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)
    
    cat(paste("\nZero-inflation model coefficients (binomial with ", x$link, " link):\n", sep = ""))
    printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
    if(x$family == "negbin"){
        cat("\nTheta:", round(x$theta, digits),"\n")} else cat("\n")
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "\n")
    cat("Converged:", x$converged, "\n")
    
    invisible(x)
}

### Based on predict.zeroinfl in R package pscl
predict.zipath <- function(object, newdata, which=1:object$nlambda, type = c("response", "prob", "count", "zero", "nonzero"),
                           na.action = na.pass, at = NULL, ...)
{
    type <- match.arg(type)
    if(type=="nonzero"){
        nbetacount=object$coefficients$count[,which]
        nbetazero=object$coefficients$zero[,which]
        if(length(which)>1)
            return(list(nbetacount=nonzeroCoef(nbetacount[-1,,drop=FALSE],bystep=TRUE),
                        nbetazero=nonzeroCoef(nbetazero[-1,,drop=FALSE],bystep=TRUE)))
        else return(list(nbetacount=which(abs(nbetacount[-1]) > 0), nbetazero=which(abs(nbetazero[-1]) > 0)))
    }
    predictzeroinfl1(object, newdata, which, type,
                      na.action = na.pass, at = at, ...)
}

### based on pscl package predict.zeroinfl function
predictzeroinfl1 <- function(object, newdata, which=1:object$nlambda, type = c("response", "prob", "count", "zero"),
                              na.action = na.pass, at = NULL, ...)
{
    type <- match.arg(type)
    object$coefficients$count <- object$coefficients$count[,which]
    object$coefficients$zero <- object$coefficients$zero[,which] 
    object$nlambda <- length(which)
    ## if no new data supplied
    if(missing(newdata)) {
        rval <- object$fitted.values[,which]
        if(type != "response") {
            if(!is.null(object$x)) {
                X <- object$x$count
                Z <- object$x$zero
            } else if(!is.null(object$model)) {
                X <- model.matrix(object$terms$count, object$model, contrasts = object$contrasts$count)
                Z <- model.matrix(object$terms$zero,  object$model, contrasts = object$contrasts$zero)
            } else {
                stop("predicted probabilities cannot be computed with missing newdata")
            }
            offsetx <- if(is.null(object$offset$count)) rep.int(0, NROW(X)) else object$offset$count
            offsetz <- if(is.null(object$offset$zero))  rep.int(0, NROW(Z)) else object$offset$zero
            mu <- exp(X %*% object$coefficients$count + offsetx)
            phi <- object$linkinv(Z %*% object$coefficients$zero + offsetz)
        }
    } else {
        mf <- model.frame(delete.response(object$terms$full), newdata, na.action = na.action, xlev = object$levels)
        X <- model.matrix(delete.response(object$terms$count), mf, contrasts = object$contrasts$count)
        Z <- model.matrix(delete.response(object$terms$zero),  mf, contrasts = object$contrasts$zero)
        offsetx <- model_offset_2(mf, terms = object$terms$count, offset = FALSE)
        offsetz <- model_offset_2(mf, terms = object$terms$zero,  offset = FALSE)
        if(is.null(offsetx)) offsetx <- rep.int(0, NROW(X))
        if(is.null(offsetz)) offsetz <- rep.int(0, NROW(Z))
        if(!is.null(object$call$offset)) offsetx <- offsetx + eval(object$call$offset, newdata)
        mu <- exp(X %*% object$coefficients$count + offsetx)
        phi <- object$linkinv(Z %*% object$coefficients$zero + offsetz)
        rval <- (1-phi) * mu
    }
    
    ## predicted means for count/zero component
    if(type == "count") rval <- mu
    if(type == "zero") rval <- phi

    ## predicted probabilities
    if(type == "prob") {
        if(!is.null(object$y)) y <- object$y
        else if(!is.null(object$model)) y <- model.response(object$model)
        else stop("predicted probabilities cannot be computed for fits with y = FALSE and model = FALSE")

        yUnique <- if(is.null(at)) 0:max(y) else at
        nUnique <- length(yUnique)
        rval <- vector("list", length=object$nlambda)
        for(j in 1:length(rval)){ 
            rval[[j]] <- matrix(NA, nrow = length(y), ncol = nUnique)
        }

        for(j in 1:object$nlambda){ 
            switch(object$family,
                   "poisson" = {
                       rval[[j]][, 1] <- phi[,j] + (1-phi[,j]) * exp(-mu[,j])
                       for(i in 2:nUnique) rval[[j]][,i] <- (1-phi[,j]) * dpois(yUnique[i], lambda = mu[,j])
                   },
                   "negbin" = {
                       theta <- object$theta[j]
                       rval[[j]][, 1] <- phi[,j] + (1-phi[,j]) * dnbinom(0, mu = mu[,j], size = theta[j])
                       for(i in 2:nUnique) rval[[j]][,i] <- (1-phi[,j]) * dnbinom(yUnique[i], mu = mu[,j], size = theta[j])
                   },
                   "geometric" = {
                       rval[[j]][, 1] <- phi[,j] + (1-phi[,j]) * dnbinom(0, mu = mu[,j], size = 1)
                       for(i in 2:nUnique) rval[[j]][,i] <- (1-phi[,j]) * dnbinom(yUnique[i], mu = mu[,j], size = 1)
                   })
        }
    }
    rval
}


fitted.zipath <- function(object, ...) {
    object$fitted.values
}

residuals.zipath <- function(object, type = c("pearson", "response"), ...) {

    type <- match.arg(type)
    res <- object$residuals

    switch(type,
           
           "response" = {
               return(res)
           },
           
           "pearson" = {
               mu <- predict(object, type = "count")
               phi <- predict(object, type = "zero")
               theta1 <- switch(object$family,
                                "poisson" = 0,
                                "geometric" = 1,
                                "negbin" = 1/object$theta)
               vv <- object$fitted.values * (1 + (phi + theta1) * mu)
               return(res/sqrt(vv))  
           })
}

terms.zipath <- function(x, model = c("count", "zero"), ...) {
    x$terms[[match.arg(model)]]
}

model.matrix.zipath <- function(object, model = c("count", "zero"), ...) {
    model <- match.arg(model)
    if(!is.null(object$x)) rval <- object$x[[model]]
    else if(!is.null(object$model)) rval <- model.matrix(object$terms[[model]], object$model, contrasts = object$contrasts[[model]])
    else stop("not enough information in fitted model to return model.matrix")
    return(rval)
}

### evaluation of penalty function
###consider enet: \lambda_1 |theta| + 1/2 * \lambda_2 theta^2, mnet and snet penalty
pen_eval <- function(theta, lone, ltwo, gamma,
                     penalty = c("enet", "mnet", "snet")){
    if(lone < 0 || ltwo < 0) stop("lone and ltwo should be non-negative values\n")
    penalty <- match.arg(penalty)
    if(!penalty %in% c("enet", "mnet", "snet"))
        stop("penalty type ", penalty, " is not implemented\n")
    if(penalty == "mnet"){
        if(gamma <= 1) stop("gamma should > 1 for mnet penenaly\n")
        if(abs(theta) <= (gamma * lone)) return(lone * abs(theta) - theta^2/(2*gamma) + 0.5*ltwo*theta^2)
        else return(0.5*gamma*lone^2 + 0.5*ltwo*theta^2)
    }
    if(penalty == "snet"){
        if(gamma <= 2) stop("gamma should > 2 for snet penalty\n")
        if(abs(theta) <= lone) return(lone*abs(theta) + 0.5*ltwo*theta^2)
        if(abs(theta) <= (gamma * lone)) return((gamma * lone * abs(theta) - 0.5*(theta^2+lone^2))/(gamma-1) + 0.5*ltwo*theta^2)
        return(lone^2*(gamma^2 - 1)/(2*(gamma-1)) + 0.5*ltwo*theta^2)
    }
    if(penalty == "enet")
        return(abs(theta)*lone + 0.5*ltwo*theta^2)
}
### evaluation of second derivative of penalty function
###consider enet: \lambda_1 |theta| + 1/2 * \lambda_2 theta^2, mnet and snet penalty
pen2_eval <- function(theta, lone, ltwo, gamma,
                      penalty = c("enet", "mnet", "snet")){
    if(lone < 0 || ltwo < 0) stop("lone and ltwo should be non-negative values\n")
    penalty <- match.arg(penalty)
    if(!penalty %in% c("enet", "mnet", "snet"))
        stop("penalty type ", penalty, " is not implemented\n")
    if(penalty == "mnet"){
        if(gamma <= 1) stop("gamma should > 1 for mnet penalty\n")
        if(abs(theta) <= (gamma * lone)) return(-1/gamma+ltwo)
        else return(ltwo)
    }
    if(penalty == "snet"){
        if(gamma <= 2) stop("gamma should > 2 for snet penalty\n")
        if(abs(theta) > lone && abs(theta) < (gamma * lone)) return(-1/(gamma-1) + ltwo)
        else return(ltwo)
    }
    if(penalty == "enet")
        return(ltwo)
}

### convert zeroinfl object to class zipath, thus can be used to predict newdata
conv2zipath <- function(object, family=c("poisson", "negbin", "geometric")){
    family <- match.arg(family)
    if(class(object)=="zeroinfl") class(object) <- "zipath"
    object$family <- family
    object$nlambda <- 1
    object$fitted.values <- matrix(object$fitted.values, ncol=1)
    namecount <- names(object$coefficients$count)
    namezero <- names(object$coefficients$zero)
    object$coefficients$count <- matrix(object$coefficients$count, ncol=1)
    object$coefficients$zero <- matrix(object$coefficients$zero, ncol=1)
    rownames(object$coefficients$count) <- namecount
    rownames(object$coefficients$zero) <- namezero
    object
}

predprob.zipath <- function(obj, which, at, ...) {
    if(length(which) > 1)
        stop("length of which should be 1\n")
    predict(obj, type = "prob", which=which, at=at, ...)[[1]]
}

hessianReg <- function(x, which, ...)
    UseMethod("hessianReg")

hessianReg.zipath <- function(x, which, log=FALSE, ...) {
    ## extract data
    Y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
    X <- model.matrix(x, model = "count")
    Z <- model.matrix(x, model = "zero")
    beta <- coef(x, model = "count", which=which)
    gamma <- coef(x, model = "zero", which=which)
    kx1 <- abs(beta) > 0
    kx2 <- which(kx1)
    kz1 <- abs(gamma) > 0
    kz2 <- which(kz1)
    X <- X[,kx2]
    Z <- Z[,kz2]
    if(length(kx2)==1) X <- matrix(X)
    if(length(kz2)==1) Z <- matrix(Z)
    beta <- beta[kx2]
    gamma <- gamma[kz2]
    if(x$family=="negbin")
        theta <- x$theta[which]

    offset <- x$offset
    if(is.list(offset)) {
        offsetx <- offset$count
        offsetz <- offset$zero
    } else {
        offsetx <- offset
        offsetz <- NULL
    }
    if(is.null(offsetx)) offsetx <- 0
    if(is.null(offsetz)) offsetz <- 0
    linkobj <- make.link(x$link)
                                        #  wts <- weights(x)
    weights <- weights(x)
    if(is.null(weights)) weights <- 1
    Y0 <- Y <= 0
    Y1 <- Y > 0
    eta <- as.vector(X %*% beta + offsetx)
    mu <- exp(eta)             ### \mu_i in paper
    etaz <- as.vector(Z %*% gamma + offsetz)
    muz <- linkobj$linkinv(etaz)
    if(x$family=="negbin"){
        g1 <- theta/(theta+mu) ### r_i in paper
        g1t <- g1^theta        ### s_i in paper
        g2 <- mu/(theta+mu)    ### t_i in paper
### warning: check if this works for other links: muz=exp(etaz)/(1+exp(etaz))?
        v <- muz               ### p_i in paper
        etaz.e <- exp(etaz)    ### \xi_i in paper
        h <- etaz.e + g1t      ### h_i in paper
        p <- dim(X)[2]
        q <- dim(Z)[2]
        rval <- matrix(NA, p+q+1, p+q+1)
### w.r.t (theta,theta)
        rv0 <- g1t*theta^3*(log(g1))^2*etaz.e+2*g1t*log(g1)*mu*etaz.e*(log(g1)+1)*theta^2 + mu^2*etaz.e*g1t*(2*log(g1)+1+(log(g1))^2)*theta+mu^2*etaz.e*g1t+mu^2*g1t^2
        rv0 <- rv0/(theta*(theta+mu)^2*h^2)
        rv1 <- trigamma(Y+theta)-trigamma(theta)+(Y+theta)/((theta+mu)^2)+1/theta*g2-1/(theta+mu)
        rval[p+q+1, p+q+1] <- sum(weights[Y0]*rv0[Y0])+sum(weights[Y1]*rv1[Y1]) 
        if(log) rval[p+q+1, p+q+1] <- theta^2*rval[p+q+1, p+q+1] # for log(theta)
### for X coef partial derivative w.r.t. (beta_j, beta_k)
        for(j in 1:p)
            for(k in j:p){
                rv <- ifelse(Y1, -theta*g2^2/mu*(theta+Y), 
                             -theta^2*mu*(-g1t*mu*etaz.e+etaz.e*g1t+g1t^2)/(((theta+mu)*h)^2)) 
                rval[j,k] <- sum(weights*X[,j]*X[,k]*rv)  
            }
### for X coef and Z coef (beta, gamma)
        for(j in 1:p)
            for(k in (p+1):(p+q)){
                rv0 <- theta*X[,j]*Z[,k-p]*g1^theta*g2*etaz.e/h^2
                rval[j,k] <- sum(weights[Y0]*rv0[Y0])     
            } 
### for (beta, theta)
        for(j in 1:p){
            rv0 <- -(theta^2*g1t*log(g1)*etaz.e+theta*g1t*log(g1)*mu*etaz.e+g1t*theta*mu*etaz.e+g1t*mu*etaz.e+g1t^2*mu)*mu
            rv0 <- rv0/(((theta+mu)*h)^2) 
            rv1 <- -mu*(-Y+mu)/((theta+mu)^2)
            rval[j, p+q+1] <- sum(weights[Y0]*rv0[Y0]*X[Y0,j]) + sum(weights[Y1]*rv1[Y1]*X[Y1,j])
            if(log) rval[j, p+q+1] <- theta*rval[j, p+q+1]  ### if theta was parameterized as log(theta), see gradfun
        }
###w.r.t. (gamma_j, gamma_k)
        for(j in (p+1):(p+q))
            for(k in j:(p+q)){
                rv <- ifelse(Y1, -v*(1-v), -v+v^2+etaz.e/h-(etaz.e/h)^2)
                rval[j,k] <- sum(weights*Z[,j-p]*Z[,k-p]*rv)  
            }
###w.r.t. {theta, gamma}
        for(k in (p+1):(p+q)){
            rv0 <- g1t*(log(g1)+g2)*etaz.e*Z[,k-p]
                                        #rv0 <- g1t*(log(g1)+(1/(theta+mu)-theta/(theta+mu)^2)*(theta+mu))*etaz.e*Z[,k-p]
            rv0 <- -rv0/h^2
            rval[k,p+q+1] <- sum(weights[Y0]*rv0[Y0])
            if(log) rval[k,p+q+1] <- theta * rval[k, p+q+1]  ### if theta was parameterized as log(theta)
        }
        rval[lower.tri(rval)] = t(rval)[lower.tri(rval)]
    }
                                        #    if(log)
                                        #        colnames(rval) <- c(names(coef(x)$count), names(coef(x)$zero), "log(theta)")
                                        #    else colnames(rval) <- c(names(coef(x)$count), names(coef(x)$zero), "theta")
                                        #else colnames(rval) <- c(names(coef(x)), "theta")
                                        #    rownames(rval) <- colnames(rval)
    return(rval)
}

### derived from sandwich function in R package sandwich
sandwichReg <- function(x, breadreg.=breadReg, meatreg.=meatReg, which, log=FALSE, ...){
    if(missing(which)) which <- 1:x$nlambda
    RET <- vector("list", length=x$nlambda)
    for(k in which){
        if(is.function(breadReg)) breadreg. <- breadreg.(x, k, log=log)
        if(is.function(meatReg)) meatreg. <- meatreg.(x, k, log=log)
        n <- x$n
        RET[[k]] <- 1/n*(breadreg. %*% meatreg. %*% breadreg.)   
    }
    RET
}

meatReg <- function(x, which, ...)
    UseMethod("meatReg")

meatReg.zipath <- function(x, which, log=FALSE, ...){
    obj1 <- x
    k <- which
    count <- as.matrix(coef(x)$count)[,k]
    zero <- as.matrix(coef(x)$zero)[,k]
    obj1$coefficients <- list(count=matrix(count, ncol=1), zero=matrix(zero, ncol=1))
    obj1$nlambda <- 1
    if(x$family=="negbin")
        obj1$theta <- x$theta[k]
    obj1$fitted.values <- x$fitted.values[,k]
    grad <- estfunReg(obj1, log=log)
    kx1 <- abs(count) > 0
    kx2 <- which(kx1)
    kz1 <- abs(zero) > 0
    kz2 <- which(kz1)
    if(x$family=="negbin") npar <- c(kx2, length(kz1)+kz2, dim(grad)[2])
    else npar <- c(kx2, length(kz1)+kz2)
    grad <- grad[,npar] 
    n <- NROW(grad)
    tmp <- apply(grad, 2, mean)
    RET <- 1/n*crossprod(grad)  - tmp %o% tmp ### note there is an additional 2nd item for penalized regression. in ordinary regression, the 2nd item=0 
                                        #        rownames(RET) <- colnames(RET)
    return(RET)
}

breadReg <- function(x, which, ...)
    UseMethod("breadReg")

breadReg.zipath <- function(x, which, log=FALSE, ...){
    family <- x$family
    n <- x$n
    kx <- dim(x$coefficients$count)[1]    
    kx1 <- abs(x$coefficients$count[,which]) > 0
    kx2 <- which(kx1)
                                        #    kz <- dim(x$coefficients$zero)[1]    
    kz1 <- abs(x$coefficients$zero[,which]) > 0
    kz2 <- which(kz1)
    sigmax <- rep(0, sum(kx1))
    sigmaz <- rep(0, sum(kz1)+(family=="negbin"))
                                        #sigma <- rep(0, sum(kx1)+sum(kz1)+(family=="negbin"))
                                        #sigma <- rep(0, kx+kz+(family=="negbin"))
    k <- which
    start <- x$coefficient
    if(length(kx2) > 1)
        for(ss in 2:length(kx2))
            sigmax[ss] <- pen2_eval(abs(start$count[kx2[ss],k]), lone=x$lambda.count[k]*x$alpha.count, ltwo=x$lambda.count[k]*(1-x$alpha.count), gamma=x$gamma.count, penalty=x$penalty)
    if(length(kz2) > 1)
        for(ss in 2:length(kz2))
            sigmaz[ss] <- pen2_eval(abs(start$zero[kz2[ss],k]), lone=x$lambda.zero[k]*x$alpha.zero, ltwo=x$lambda.zero[k]*(1-x$alpha.zero), gamma=x$gamma.zero, penalty=x$penalty)
    sigma <- c(sigmax, sigmaz)
    hes <- hessianReg(x, which, log=log)
    s1 <- try(-solve(hes - n*diag(sigma))*n)
    return(s1)
}

estfunReg <- function(x, ...)
    UseMethod("estfunReg")

estfunReg.zipath <- function(x, which=1, log=FALSE, ...) {
    ## extract data
    if(length(which)!=1) stop("argument which should be one integer only\n")
    Y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
    X <- model.matrix(x, model = "count")
    Z <- model.matrix(x, model = "zero")
    beta <- coef(x, model = "count")
    gamma <- coef(x, model = "zero")
    if(x$nlambda > 1){
     beta <- beta[,which]
     gamma <- gamma[,which]
    }
    if(x$family=="negbin")
      theta <- x$theta[which]
    offset <- x$offset
    if(is.list(offset)) {
        offsetx <- offset$count
        offsetz <- offset$zero
    } else {
        offsetx <- offset
        offsetz <- NULL
    }
    if(is.null(offsetx)) offsetx <- 0
    if(is.null(offsetz)) offsetz <- 0
    linkobj <- make.link(x$link)
    wts <- weights(x)
    if(is.null(wts)) wts <- 1
    Y1 <- Y > 0
    eta <- as.vector(X %*% beta + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% gamma + offsetz)
    muz <- linkobj$linkinv(etaz)

    ## density for y = 0
    clogdens0 <- switch(x$family,
                        "poisson" = -mu,
                        "geometric" = dnbinom(0, size = 1, mu = mu, log = TRUE),
                        "negbin" = dnbinom(0, size = theta, mu = mu, log = TRUE))
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)

    ## working residuals  
    wres_count <- switch(x$family,
                         "poisson" = ifelse(Y1, Y - mu, -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(mu))),
                         "geometric" = ifelse(Y1, Y - mu * (Y + 1)/(mu + 1), -exp(-log(dens0) +
                             log(1 - muz) + clogdens0 - log(mu + 1) + log(mu))),
                         "negbin" = ifelse(Y1, Y - mu * (Y + theta)/(mu + theta), -exp(-log(dens0) +
                             log(1 - muz) + clogdens0 + log(theta) - log(mu + theta) + log(mu))))
    wres_zero <- ifelse(Y1, -1/(1-muz) * linkobj$mu.eta(etaz),
                        (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
    wres_theta <- ifelse(Y1, digamma(Y + theta) - digamma(theta) +
                         log(theta) - log(mu + theta) + 1 - (Y + theta)/(mu + theta),
                         exp(-log(dens0) + log(1 - muz) + clogdens0) *
                         (log(theta) - log(mu + theta) + 1 - theta/(mu + theta)))

    ## compute gradient from data
    if(log)
        rval <- cbind(wres_count * wts * X, wres_zero * wts * Z, theta*wres_theta)
    else rval <- cbind(wres_count * wts * X, wres_zero * wts * Z, wres_theta)
    if(x$family!="negbin")
        colnames(rval) <- c(names(coef(x)$count), names(coef(x)$zero))
    else{
                                        #       	  if(log)
                                        #      colnames(rval) <- c(names(coef(x)$count), names(coef(x)$zero), "log(theta)")
                                        #   else colnames(rval) <- c(names(coef(x)$count), names(coef(x)$zero), "theta")
    }
                                        #  rownames(rval) <- rownames(X)
    return(rval)
}

se <- function(x, which, log=TRUE, ...)
    UseMethod("se")

se.zipath <- function(x, which, log=TRUE, ...){
    if(is.null(which)) stop("which has no numerical value\n")
    if(length(which) > 1) stop("which is an integer value\n")
    tmp <- sandwichReg(x, which=which, log=log)
    sqrt(diag(tmp[[which]]))
}


