#######################################################################
# revised version of file MASS/polr.q:
#   [ original copyright (C) 1994-2005 W. N. Venables and B. D. Ripley ]
#   [ original 18266 bytes, Sep  5 2006 09:03 ]
# 
# Modifications by Jonathan Wand as follows
# ON 2006-11-03: JW
# - Response variable:
#   = now requires N x 2 matrix as LHS variable, e.g., cbind( from.y1, to.y2 )
#   = response variable need not be a factor, indeed I get rid of factors if they are used
#   = allows for interval values in dependent variable (from.y1 != to.y2)
#   = allows for HOLES in dependent variable (non-observed response categories)
#     with fitted values returning zero probabilities for non-observed response categories
#     (this is a bit intricate to engineer)
# - altered default method list: probit is now default, Hess=TRUE
# ON 2006-11-04: JW
# - added 'theta' to return list, and used in vcov.cpolr [ I disagree with use of zeta as-is in original ]
# ON 2008-05-01: JW
# - replace stop() with warning() and a NULL return when nlev <= 2
# - consider adding probit'y subroutine in this case in future, but skipped for now
#######################################################################
cpolr <- function(formula, data, weights, start, ..., subset,
                  na.action, contrasts = NULL, Hess = TRUE,
                  model = TRUE,
                  method = c("probit", "logistic", "cloglog", "cauchit"),
                  debug = 0)
{
  if (debug > 0) {cat("In cpolr\n") }
  
    logit <- function(p) log(p/(1 - p))

    fmin <- function(beta) {
#      cat("fmin",pc,q,"\n")
      
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
      
        eta <- offset
        if (pc > 0)
            eta <- eta + drop(x %*% beta[1:pc])

#      print(length(C[,1]))
#      print(length(eta))
#      print(length(gamm[C[,1]]))
#      print(gamm)
#      print(table(C[,1]))
      
        pr <- pfun(gamm[ C[,2] + 1] - eta) - pfun(gamm[ C[,1] ] - eta)
        if (all(pr > 0))
            -sum(wt * log(pr))
        else Inf
    }

    gmin <- function(beta)
    {
        jacobian <- function(theta) { ## dgamma by dtheta matrix
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0 , k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
        }
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
        eta <- offset
        if(pc > 0) eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[ C[,2] +1] - eta) - pfun(gamm[ C[,1] ] - eta)
        p1 <- dfun(gamm[ C[,2] +1] - eta)
        p2 <- dfun(gamm[ C[,1] ] - eta)
        g1 <- if(pc > 0) t(x) %*% (wt*(p1 - p2)/pr) else numeric(0)
        xx <- .polrY1*p1 - .polrY2*p2
        g2 <- - t(xx) %*% (wt/pr)
        g2 <- t(g2) %*% jacobian(theta)
        if(all(pr > 0)) c(g1, g2) else rep(NA, pc+q)
    }

    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    pfun <- switch(method, logistic = plogis, probit = pnorm) #  loglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm) # cloglog = dgumbel, cauchit = dcauchy)
    if(is.matrix(eval.parent(m$data)))
        m$data <- as.data.frame(data)
#    cat("D1",dim(m$data),"\n")
    m$start <- m$Hess <- m$method <- m$debug <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts") # will get dropped by subsetting
    if(xint > 0) {
        x <- x[, -xint, drop=FALSE]
        pc <- pc - 1
    } else warning("an intercept is needed and assumed")

    if (debug > 0) cat("X1",dim(x),"\n")
    if (debug > 1) print(x)

    wt <- model.weights(m)
    if(!length(wt)) wt <- rep(1, n)
    offset <- model.offset(m)
    if(length(offset) <= 1) offset <- rep(0, n)
#    cat("X2",dim(x),"\n")

    ###################################################################################
    ## For anchors:
    ## - two column response, non-factor, possible missing columns
    ## - remapping of holes done here...
    ###################################################################################
    C <- model.response(m)
    ## if not 2 columns, then make it so...
    if (NCOL(C) == 1) {
      C <- cbind(C,C)
    }
#    ## get rid of factors if they exist
#    if (is.factor( C ) ) {
#      C <- matrix(  unclass( factor( c(C[,1],C[,2]))) , ncol=2)
#    }
    idx <- C[,1] == C[,2]
#    if(!is.factor(y)) stop("response must be a factor")
#    lev <- levels(y)
    lev.orig <- sort(unique(as.numeric(C[idx,])))
    n.lev.orig <- length(lev.orig)
    max.lev.orig  <- max(lev.orig)
    lev.missing <- !(1:max.lev.orig %in% lev.orig)
    lev <- 1:sum(!lev.missing) ## collapse out any missing categories, and renumber
#    print(lev)
#    print(lev.orig)
#    print(lev.missing)
    
    if(length(lev) <= 2) {
#      cat("cpolr: response must have 3 or more levels, skipping.\n")
      warning("cpolr response must have 3 or more levels, skipping.\n")
      return(NULL)
    }
    if( any(lev.missing) ) {
      cat("\nThere is at least one missing category:",which(lev.missing),"\n\n")
      lev.map <- as.numeric(!lev.missing) ## zero: hole
      lev.map[!lev.missing] <- lev ## put in new mapping values

      lev.map.ub <- lev.map.lb <- lev.map
      ## now tricky stuff... fill in zeroes
      ## for C[,2] use, get upper bound.. which is the last smaller category
      lev.last <- 1
      for (ii in 1:length(lev.map)) {
        if (lev.map[ii] != 0) lev.last <- lev.map[ii]
        lev.map.ub[ii] <- lev.last
      }
      ## for C[,2] use, get upper bound.. which is the last smaller category
      lev.last <- max(lev.map)
      for (ii in seq(length(lev.map),1,-1)) {
        if (lev.map[ii] != 0) lev.last <- lev.map[ii]
        lev.map.lb[ii] <- lev.last
      }
      
      C.orig <- C ## keep a copy      
      C[,1]  <- lev.map.lb[ C[,1] ] ## remap
      C[,2]  <- lev.map.ub[ C[,2] ] ## remap

      if (debug > 0) { cat("lev.map\n"); print(lev.map) }
      if (debug > 1) { cat("C.orig, C\n"); print( cbind( C.orig[,1], C[,1], C.orig[,2], C[,2])) }

    } else {
      if (debug > 1) { cat("C\n"); print(C ) }
    }
###################################################################################
    
#    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) ==  C[,2] 
    .polrY2 <- col(Y) ==  C[,1]  - 1
    if(missing(start)) {
        # try logistic/probit regression on 'middle' cut
        q1 <- length(lev) %/% 2
        y1 <- (C[,1] > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <-
            switch(method,
                   "logistic"= glm.fit(X, y1, wt, family = binomial(), offset = offset),
                   "probit" = glm.fit(X, y1, wt, family = binomial("probit"), offset = offset),
                   ## this is deliberate, a better starting point
                   "cloglog" = glm.fit(X, y1, wt, family = binomial("probit"), offset = offset),
                   "cauchit" = glm.fit(X, y1, wt, family = binomial("cauchit"), offset = offset))
        if(!fit$converged)
            stop("attempt for find suitable starting values failed")
        coefs <- fit$coefficients
        if(any(is.na(coefs))) {
            warning("design appears to be rank-deficient, so dropping some coefs")
            keep <- names(coefs)[!is.na(coefs)]
            coefs <- coefs[keep]
            x <- x[, keep[-1], drop = FALSE]
            pc <- ncol(x)
        }
        spacing <- logit((1:q)/(q+1)) # just a guess
        if(method != "logit") spacing <- spacing/1.7
        gammas <- -coefs[1] + spacing - spacing[q1]
        thetas <- c(gammas[1], log(diff(gammas)))
        start <- c(coefs[-1], thetas)
    } else if(length(start) != pc + q)
	stop("'start' is not of the correct length")
    res <- optim(start, fmin, gmin, method="BFGS", hessian = Hess, ...)
    beta <- res$par[seq_len(pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1],exp(theta[-1])))
    deviance <- 2 * res$value
    niter <- c(f.evals=res$counts[1], g.evals=res$counts[2])
    names(zeta) <- paste(lev.orig[-length(lev.orig)], lev.orig[-1], sep="|")
    if(pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    } else {
        eta <- rep(0, n)
    }

    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow=TRUE) - eta), , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)

    ## rebuild fitted values if there were holes in responses
    if (any(lev.missing)) {
      fitted2 <- matrix(0,n,max.lev.orig) 
      for (i in 1:max.lev.orig) {
        if ( lev.map[i] == 0) next
        ## get the estimated/collapsed category fitted value and put it in the right place
        fitted2[,i] <- fitted[ , lev.map[i] ]
      }
      fitted <- fitted2
      dimnames(fitted) <- list(row.names(m), 1:max.lev.orig)
    }
    
    ## 
    fit <- list(coefficients = beta, zeta = zeta, theta=theta, deviance = deviance,
                fitted.values = fitted,
                lev = lev.orig,
                lev.missing = lev.missing,
                terms = Terms,
                df.residual = sum(wt) - pc - q,
                edf = pc + q,
                n = sum(wt),
                nobs = sum(wt),
                call = match.call(), method = method,
		convergence = res$convergence, niter = niter)
    if(Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if(model) fit$model <- m
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- c("cpolr","polr")
    return(fit)
}

vcov.cpolr <- function(object, ...)
{
  ## for the ginv function...
	# Now loaded by Dependencies
  	# require(MASS)
  if(is.null(object$Hessian)) {
    cat("\nRe-fitting cpolr() to get Hessian\n\n")
    flush.console
    object <- update(object, Hess=TRUE,
                     start=c(object$coef, object$theta))
  }
  structure(ginv(object$Hessian), dimnames = dimnames(object$Hessian))
}



