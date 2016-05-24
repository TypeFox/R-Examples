#' Estimate average hazard ratios from k independent samples
#'
#' This function is a simple wrapper for \code{\link{ahrWKM}} and \code{\link{ahrAJ}}.
#' 
#' @title avgHR
#' @param L time-limit specifying time-interval [0,L] on which average hazard ratios will be calculated
#' @param data data frame (see \code{data} argument to \code{\link{ahrWKM}} (if \code{method == "wkm" || "km"}) or \code{\link{ahrAJ}} (if \code{method == "aj"}))
#' @param method method used for estimating survival functions (default: Kaplan-Meier estimator)
#' @param ... additional arguments passed to \code{\link{ahrWKM}} or \code{\link{ahrAJ}}
#' @return An object of class '"ahr"'
#' @seealso \code{\link{ahrWKM}}, \code{\link{ahrAJ}}
#' @references J.~D. Kalbfleisch and R.~L. Prentice. Estimation of the average hazard ratio. \emph{Biometrika}, 68(1):105--112, Apr. 1981.
#' @export
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100))
#' fit <- avgHR(2, data.frame(Y=Y, D=D, Z=Z), formula=Surv(Y, D) ~ Z)
avgHR <- function(L, data, method="km", ...) {
    if(method %in% c("km", "wkm")) ahrWKM(L=L, data=data, ...)
    else if(method == "aj") ahrAJ(L=L, data=data, ...)
    else stop(paste("Unknown method:", method))
}

#' Estimate average hazard ratios from k independent samples based on the weighted Kaplan-Meier (WKM) estimator
#'
#' @title ahrWKM
#' @param L time-limit specifying time-interval [0,L] over which average hazard ratios will be calculated
#' @param formula an object of class '"formula"' specifying the conditional survival model
#' @param data data frame containing the variables in formula
#' @param null.theta vector specifying the null hypothesis for the average hazard ratios (H_0: theta = null.theta)
#' @param contrast vector of contrasts to test H_0: contrast * (theta - null.theta) = 0
#' @param multi.test calculate multivariate test statistic if TRUE
#' @param cov if TRUE calculate covariance matrix estimator (direct)
#' @param bootstrap if > 0 then use bootstrap to estimate covariance matrix (ignore if cov is TRUE)
#' @param alpha exponent of the weight function
#' @param left.limit if TRUE use left-continuous interpolation of WKM estimates instead of right-continuous interpolation
#' @param rr.subset logical vector defining subset of observations to use for response rate estimation (default: use all observations)
#' @return An object of class '"ahr"'
#' @references J.~D. Kalbfleisch and R.~L. Prentice. Estimation of the average hazard ratio. \emph{Biometrika}, 68(1):105--112, Apr. 1981.
#' @export
#' @seealso \code{\link{wkm}}
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100)) # treatment indicator
#' fit <- ahrWKM(2, Surv(Y, D) ~ Z, data.frame(Y=Y, D=D, Z=Z))
#' fit
#'
#' ## the same as above, but estimate covariance matrix using bootstrap
#' \dontrun{fitBS <- ahrWKM(2, Surv(Y, D) ~ Z, data.frame(Y=Y, D=D, Z=Z), cov=FALSE,
#'                          bootstrap=1000)}
#' 
ahrWKM <- function(L, formula, data, null.theta=NULL, contrast=NULL, multi.test=FALSE, cov=TRUE, bootstrap=0, alpha=1, left.limit=FALSE, rr.subset=rep(TRUE, nrow(data))) {
    if(!is.null(formula)) data <- parseFormula(formula, data)
    
    wkm.param <- list(alpha=alpha, var=cov, cov=cov, left.limit=left.limit, rr.subset=rr.subset)

    fit <- ahrSurv(L, data, null.theta, contrast, multi.test, cov, bootstrap, wkm, wkm.param)
    fit <- c(fit, logHR(fit))
    class(fit) <- "ahr"
    fit
}

#' Estimate average hazard ratios from k independent samples based on the Kaplan-Meier estimator
#'
#' @title ahrKM
#' @param L time-limit specifying time-interval [0,L] over which average hazard ratios will be calculated
#' @param formula an object of class '"formula"' specifying the conditional survival model
#' @param data data frame containing the variables in formula
#' @param null.theta vector specifying the null hypothesis for the average hazard ratios (H_0: theta = null.theta)
#' @param contrast vector of contrasts to test H_0: contrast * (theta - null.theta) = 0
#' @param multi.test calculate multivariate test statistic if TRUE
#' @param cov if TRUE calculate covariance matrix estimator (direct)
#' @param bootstrap if > 0 then use bootstrap to estimate covariance matrix (ignore if cov is TRUE)
#' @param left.limit if TRUE use left-continuous interpolation of WKM estimates
#' @return An object of class '"ahr"'
#' @references  J.~D. Kalbfleisch and R.~L. Prentice. Estimation of the average hazard ratio. \emph{Biometrika}, 68(1):105--112, Apr. 1981.
#' @export
#' @seealso \code{\link{survfit}}
#' @examples
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' Y <- pmin(T, C)
#' D <- T <= C
#' Z <- rep(c(0,1), c(100, 100)) # treatment indicator
#' fit <- ahrKM(2, Surv(Y, D) ~ Z, data.frame(Y=Y, D=D, Z=Z))
#' fit
#'
#' ## the same as above, but estimate covariance matrix using bootstrap
#' \dontrun{fitBS <- ahrKM(2, Surv(Y, D) ~ Z, data.frame(Y=Y, D=D, Z=Z), cov=FALSE,
#'                          bootstrap=1000)}
ahrKM <- function(L, formula, data, null.theta=NULL, contrast=NULL, multi.test=FALSE, cov=TRUE, bootstrap=0, left.limit=FALSE) {
    if(!is.null(formula)) data <- parseFormula(formula, data)
    
    wkm.param <- list(alpha=1, var=cov, cov=FALSE, left.limit=left.limit, rr.subset=rep(TRUE, nrow(data)))
    
    fit <- ahrSurv(L, data, null.theta, contrast, multi.test, cov, bootstrap, wkm, wkm.param, TRUE)
    fit <- c(fit, logHR(fit))
    class(fit) <- "ahr"
    fit
}

#' Estimate average hazard ratios from k independent samples based on the Aalen-Johansen estimator of the empirical transition probabilities (NOTE: variance estimation not yet implemented)
#'
#' @title ahrAJ
#' @param L time-limit specifying time-interval [0,L] over which average hazard ratios will be calculated
#' @param target string specifying the target transition, for which the Aalen-Johansen estimator is to be calculated
#' @param states list of state names
#' @param transitions matrix of possible transitions
#' @param censoring name of censoring 'state'
#' @param data data frame containing variables id, time, from, to (see \code{\link{etm}}) and Trt (factor giving treatment groups)
#' @param null.theta vector specifying the null hypothesis for the average hazard ratios
#' @param contrast vector of contrasts to test H_0: contrast * (theta - null.theta) = 0
#' @param multi.test calculate multivariate test statistic if TRUE
#' @param cov if TRUE calculate covariance matrix estimator (direct)
#' @param bootstrap number of bootstrap samples to draw for variance estimation (default: 0 = no bootstrap, direct variance estimation). This parameter is ignored if cov=TRUE
#' @return An object of class '"ahr"'
#' @seealso \code{\link{aj}}
#' @references J.~D. Kalbfleisch and R.~L. Prentice. Estimation of the average hazard ratio. \emph{Biometrika}, 68(1):105--112, Apr. 1981.
#' @export
#' @examples
#' ## competing risks
#' Trt <- factor(rep(c(0,1), c(100, 100)))
#' T <- c(rexp(100, 1), rexp(100, 2))
#' C <- c(rexp(100, 1), rexp(100, 2))
#' r <- c(rbinom(100, 2, 0.5), rbinom(100, 2, 0.4))
#' r[(r == 0) | (T > C)] <- "cens"
#' data <- data.frame(id=1:200, time=pmin(T,C), from=rep(0, 200), to=r, Trt=Trt)
#' tra <- matrix(FALSE, nrow=3, ncol=3)
#' tra[1, 2:3] <- TRUE
#' # estimate average subdistribution hazard ratio up to L=2 for event type 1
#' fit <- ahrAJ(2, target="0 1", states=c("0", "1", "2"), transitions=tra, censoring="cens",
#'              data=data, cov=TRUE)
#' fit
ahrAJ <- function(L, target, states, transitions, censoring, data, null.theta=NULL, contrast=NULL, multi.test=FALSE, cov=FALSE, bootstrap=0) {
    
    msm <- list(target=target, states=states, transitions=transitions, censoring=censoring, s=0, t=L)

    if(cov) {
        if(isCmprsk(msm$transitions)) msm$covariance <- TRUE
        else stop("Variance estimation only supported for competing risks models")
        if(bootstrap > 0) warning("Bootstrap parameter ignored since 'cov' is TRUE")
    } else msm$covariance <- FALSE

    fit <- ahrSurv(L, data, null.theta, contrast, multi.test, cov, bootstrap, aj, msm, FALSE)

    fit <- c(fit, logHR(fit))
    class(fit) <- "ahr"
    fit
}

#' Print ahr object
#'
#' @title print.ahr
#' @param x an object of class '"ahr"'.
#' @param digits minimal number of significant digits.
#' @param ... further arguments passed to or from other methods.
#' @method print ahr
#' @export
print.ahr <- function(x, digits=3, ...) {
    cat("\n")
    cat("\t Average Group to Total Hazard Ratios ", paste0("(L=", x$L, ")"), "\n")
    cat("\n")

    cat(paste(x$k, "treatment groups:"), levels(x$groups), "\n")

    cat("\n")
    
    if(is.null(x$cov.theta)) tmpcov <- rep(NA, x$k) ## covariance matrix missing
    else if(is.null(dim(x$cov.theta))) tmpcov <- x$cov.theta  ## scalar value
    else tmpcov <- diag(x$cov.theta)

    if(is.null(x$Z.theta)) Ztheta <- NA
    else Ztheta <- round(x$Z.theta, digits)
    
    tmp <- cbind(round(x$theta, digits), round(sqrt(tmpcov), digits), Ztheta)
    dimnames(tmp) <- list(levels(x$groups), c("Theta", "Std.dev", "Z"))

    print(tmp)

    cat("\n\n")
    cat("\t Average Hazard Ratios ", paste0("(L=", x$L, ")"), "\n")
    cat("\n")

    if(is.null(x$cov.hr)) tmpcov <- rep(NA, x$k-1)
    else if(is.null(dim(x$cov.hr))) tmpcov <- x$cov.hr
    else tmpcov <- diag(x$cov.hr)

    cat("\n")

    if(is.null(x$Z.hr)) Zhr <- NA
    else Zhr <- round(x$Z.hr, digits)
    
    tmp <- cbind(round(x$hr, digits), round(sqrt(tmpcov), digits), Zhr)
    dimnames(tmp) <- list(levels(x$groups)[-1], c("HR", "Std.dev", "Z"))

    print(tmp)  
    
    cat("\n")

    if(!is.null(x$contrast)) cat
    
    invisible(x)
}

## INTERNAL
##
## title ahrSurv
## param data data frame containing the variables
## param L time-limit specifying time-interval [0,L] over which average hazard ratios will be calculated
## param null.theta vector specifying the null hypothesis for the average hazard ratios
## param contrast vector of contrasts to test H_0: contrast * (theta - null.theta) = 0
## param multi.test calculate multivariate test statistic if TRUE
## param cov if TRUE calculate asymptotic covariance matrix estimator (direct estimation)
## param bootstrap if > 0 then use bootstrap to estimate covariance matrix (ignored if cov is TRUE)
## param surv.fit.fun function which calculates survival function from data
## param surv.fit.param list of parameters passed to surv.fit.fun
## param log.iis asymptotic covariance function of logarithm of survival function has independent increments structure (TRUE for (weighted) Kaplan-Meier estimator, FALSE for Aalen-Johansen estimator)
## return an object of class '"ahr"'
ahrSurv <- function(L, data, null.theta=NULL, contrast=NULL, multi.test=FALSE, cov=TRUE, bootstrap=0, surv.fit.fun=wkm, surv.fit.param=list(alpha=1, cov=TRUE), log.iis=FALSE) {
    
    ahr.obj <- ahrFit(data, L, cov, surv.fit.fun, surv.fit.param, log.iis)

    if(!cov && bootstrap > 0) ahr.obj$cov.theta <- ahrBootstrap(ahr.obj, data, bootstrap, surv.fit.fun, surv.fit.param)$cov.theta
    if(cov && bootstrap > 0) warning("Bootstrap parameter ignored since 'cov' is TRUE")

    k <- ahr.obj$k
    n <- ahr.obj$n
    theta <- ahr.obj$theta
    cov.theta <- ahr.obj$cov.theta

    if(is.null(null.theta)) null.theta <- rep(1/k, k)
    ahr.obj$null.theta <- null.theta

    if(!is.null(cov.theta)) {
        ## Wald-type test
        Z.theta <- (theta - null.theta) / sqrt(diag(cov.theta))
        
        ahr.obj$cov.theta <- cov.theta
        ahr.obj$Z.theta <- Z.theta
    
        ## linear contrast test
        if(!is.null(contrast)) {
            ahr.obj$contrast <- contrast
            ahr.obj$Z.contrast <- contrast %*% (theta - null.theta) / sqrt(contrast %*% cov.theta %*% contrast)
        }        
        
        if(multi.test) {
            ahr.obj$multi.test <- multi.test
            ahr.obj$Z.multi <- (theta - null.theta) %*% MASS::ginv(cov.theta) %*% (theta - null.theta) ## has chi^2 distribution with k-1 degrees of freedom
        }
    }
    ahr.obj
}

## INTERNAL
## bootstrap estimation of covariance matrix
ahrBootstrap <- function(ahr.obj, data, B, surv.fit.fun, surv.fit.param) {

    k <- ahr.obj$k
    lev <- levels(ahr.obj$groups)
    L <- ahr.obj$L
    sel <- list()
    n <- numeric(k)

    if(is.null(data$id)) data$id <- 1:nrow(data)
    
    ## sel <- lapply(1:k, function(i) unique(data$id[data$Trt == lev[i]]))
    ## n <- do.call(c, lapply(sel, function(x) length(x)))
    for(i in 1:k) {
        sel[[i]] <- unique(data$id[data$Trt == lev[i]])
        n[i] <- length(sel[[i]])
    }

    surv.fit.param$var <- FALSE
    surv.fit.param$cov <- FALSE
    
    f <- function() {
        ## stratified resampling
        ids <- do.call(c, lapply(1:k, function(i) sample(sel[[i]], n[i], replace=TRUE)))
    
        ## draw bootstrap sample
        Bdata <- data[data$id %in% ids,]
        Bdata$id <- 1:nrow(Bdata) ## ids are not unique after sampling WITH replacement
        
        ## calculate ahr from bootstrap sample
        ahrFit(Bdata, L, FALSE, surv.fit.fun, surv.fit.param)$theta
    }

    res <- replicate(B, f())

    list(theta=rowMeans(res), cov.theta=cov(t(res)))
}

## INTERNAL
ahrFit <- function(data, L, cov, surv.fit.fun, surv.fit.param, log.iis=FALSE) {
    
    grps <- levels(data$Trt)
    k <- length(grps)
    if(k <= 1) stop("Need at least two groups!")

    times <- getTimes(L, data)
    n.times <- length(times)

    n <- NULL
    fit <- list()

    trt.sub <- data$Trt[surv.fit.param$rr.subset]
    
    ## estimate transition probabilities in each treatment group
    for(i in grps) {
        cur.data <- data[data$Trt == i,]

        par <- surv.fit.param
        par$rr.subset <- surv.fit.param$rr.subset[trt.sub == i]
        
        n <- c(n, nrow(cur.data))
        fit[[length(fit)+1]] <- surv.fit.fun(times, cur.data, par)
    }

    p <- n / sum(n)

    tmp <- lapply(1:k, function(i) fit[[i]]$S)
    
    ## dummy to ensure length(tmp[-l]) > 1
    tmp[[k+1]] <- rep.int(1, n.times)
    
    xi <- sapply(1:k, function(l) do.call(function(...) mapply(prod, ...), tmp[-l]))
    G <- xi[,k] * fit[[k]]$S                
    x <- sapply(1:(k-1), function(i) stepIntegrate(xi[,i], fit[[i]]$S))
    
    theta <- -x / (1 - G[n.times])
    theta[k] <- 1 - sum(theta)

    ahr.obj <- list(L=L, surv.fit.param=surv.fit.param, k=k, n=n, p=p,
                    surv.fit=fit, times=times, n.times=n.times, groups=data$Trt, strata=data$W, log.iis=log.iis)
    class(ahr.obj) <- "ahr"
    
    ahr.obj$theta <- theta    
    if(cov) ahr.obj$cov.theta <- ahrVar(theta, G, xi, times, p, fit, log.iis) / sum(n)

    ahr.obj
}

## INTERNAL
ahrVar <- function(theta, G, xi, times, p, fit, log.iis=FALSE) {
    k <- length(fit)
    GL <- G[length(G)]
    n.times <- length(times)
    snt <- 1:n.times

    ## logCOV = COV / S
    ## COV == 0 <=> S = 0
    get.A <- function(i, j) {
        if(log.iis) ft <- xi[,j] * fit[[i]]$logV
        else ft <- xi[,j] * fit[[i]]$logCOV[ , ncol(fit[[i]]$logCOV)]
        GL * stepIntegrate(ft, fit[[j]]$S) / p[i]
    }
    
    get.B <- function(i, j, kk) {
        Sj <- fit[[j]]$S
        if(log.iis) {
            logV <- fit[[i]]$logV
            tmp <- simplify2array(lapply(snt, function(l) stepIntegrate(xi[,j] * logV[pmin(snt, l)], Sj)))
        } else {
            logCOV <- fit[[i]]$logCOV
            tmp <- simplify2array(lapply(snt, function(l) stepIntegrate(xi[,j] * logCOV[, l], Sj)))
        }
        stepIntegrate(xi[,kk] * tmp, fit[[kk]]$S) / p[i]
    }
    
    get.C <- function(i) {
        if(log.iis) GL^2 * fit[[i]]$logV[n.times] / p[i]
        else GL^2 * fit[[i]]$logCOV[n.times, n.times] / p[i]
    }

    A <- matrix(0, nrow=k, ncol=k)
    B <- array(0, dim=rep.int(k, 3))
    C <- numeric(k)
    sel <- 1:k
    for(i in sel) {
        C[i] <- get.C(i)
        for(j in sel) {
            if(i != j) A[i,j] <- get.A(i,j)
            for(l in j:k) {
                B[i,j,l] <- get.B(i,j,l)
                B[i,l,j] <- B[i,j,l]
            }
        }
    }

    VG <- sum(C)
    VxG <- numeric(k)
    Sigma <- matrix(0, nrow=k, ncol=k)

    for(i in sel) {
        si <- sel[-i]
        Vxii <- sum(B[i, si, si]) + sum(B[si, i, i]) + C[i] - 2*sum(A[i, si])
        VxG[i] <- sum(A[i, si] - A[si, i]) - C[i]
        Sigma[i,i] <- Vxii + 2*theta[i]*VxG[i] + theta[i]^2 * VG
    }

    for(i in sel[-k]) {
        si <- sel[-i]
        for(j in sel[(i+1):k]) {
            Vxij <- -sum(B[i, si, j]) - sum(B[j, sel[-j], i]) + sum(B[sel[-c(i,j)], i, j]) + A[i,j] + A[j,i]
            Sigma[i,j] <- Vxij + theta[j]*VxG[i] + theta[i]*VxG[j] + theta[i]*theta[j]*VG
            Sigma[j,i] <- Sigma[i,j]
        }
    }

    ## must be zero:
    ##print(rowSums(Sigma))
    ##print(colSums(Sigma))
    
    Sigma / (1 - GL)^2
}

## INTERNAL
logHR <- function(ahr.obj, null.hr=1, contrast=NULL) {
    k <- ahr.obj$k
    theta <- ahr.obj$theta
    cov.theta <- ahr.obj$cov.theta
    null.loghr <- log(null.hr)
    
    hr <- theta[-1] / theta[1]
    loghr <- log(hr)

    res <- list(hr=hr, loghr=loghr, null.hr=null.hr, null.loghr=null.loghr)
    
    if(!is.null(cov.theta)) {
        ## covariance matrix cov.beta has dimension k-1
        if(k == 2) {
            cov.loghr <- cov.theta[1,1] / (theta[1] * theta[2])^2
            Z.loghr <- (loghr - null.loghr) / sqrt(cov.loghr)

            cov.hr <- cov.loghr * hr^2
            Z.hr <- (hr - null.hr) / sqrt(cov.hr)
        } else {
            J <- matrix(0, nrow=k-1, ncol=k)
            J[,1] <- -1 / theta[1]
            J[,-1] <- diag(1 / theta[-1])
            cov.loghr <- J %*% cov.theta %*% t(J)
            Z.loghr <- sum(solve(cov.loghr, loghr - null.loghr) * (loghr - null.loghr))  ## asympt. chi^2 distribution with k-1 degrees of freedom

            J2 <- diag(loghr)
            cov.hr <-  J2 %*% cov.loghr %*% J2
            Z.hr <- sum(solve(cov.hr, hr - null.hr) * (hr - null.hr))  ## asympt. chi^2 distribution with k-1 degrees of freedom
        }       

        res$Z.loghr <- Z.loghr
        res$cov.loghr <- cov.loghr

        res$Z.hr <- Z.hr
        res$cov.hr <- cov.hr
        
        if(!is.null(contrast)) { ## asympt. standard normal distribution
            res$contrast <- contrast
            res$Z.contrast.loghr <- sum(contrast * (loghr - null.loghr)) / sqrt(contrast %*% cov.loghr %*% contrast)
            res$Z.contrast.hr <- sum(contrast * (hr - null.hr)) / sqrt(contrast %*% cov.hr %*% contrast)
        }
    }
    res
}
