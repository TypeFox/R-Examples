svyolr<-function(formula, design,...) UseMethod("svyolr",design)

##
## Much of this is taken from MASS polr, thus the GPL-2 license.
##


svyolr.svyrep.design<-function(formula,design,...,return.replicates=FALSE,
                               multicore=getOption("survey.multicore")){
 	require(MASS)
 	environment(formula)<-environment()
 	df<-model.frame(design)
 	pwt<-weights(design,"sampling")
        if (multicore && !require("parallel", quietly=TRUE))
          multicore <- FALSE
        
 	rval<-suppressWarnings(polr(formula,data=df,...,Hess=TRUE,model=FALSE,
                                    weights=pwt))
 	start<-c(rval$coefficients,rval$zeta)
 	rw<-weights(design,"analysis")
        if (multicore){
          betas<-do.call(cbind,mclapply(1:ncol(rw), function(i){
            w<-rw[,i]
            environment(formula)<-environment()
            m<-polr(formula,data=df,Hess=FALSE, start=start, model=FALSE, weights=w)
            c(m$coefficients, m$zeta)
          }))
        } else {
          betas<-apply(rw,2,function(w) {
            environment(formula)<-environment()
            m<-polr(formula,data=df,Hess=FALSE, start=start, model=FALSE, weights=w)
            c(m$coefficients, m$zeta)
          })
        }
 	rval$var<-svrVar(t(betas),design$scale,design$rscales,mse=design$mse, coef=start)
        rval$df.residual<-degf(design)-length(rval$coefficients)
        rval$deviance<-rval$deviance/mean(pwt)
 	class(rval)<-"svyolr"
 	rval$call<-sys.call()
        rval$call[[1]]<-as.name(.Generic)
 	if (return.replicates) rval$replicates<-t(betas)
 	rval
 	}


pgumbel<-
function (q, loc = 0, scale = 1, lower.tail = TRUE) 
{
    q <- (q - loc)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 
        1 - p
    else p
}
dgumbel<-function (x, loc = 0, scale = 1, log = FALSE) 
{
    x <- (x - loc)/scale
    d <- log(1/scale) - x - exp(-x)
    if (!log) 
        exp(d)
    else d
}


svyolr.survey.design2<-function (formula, design,  start, ...,  na.action=na.omit, 
 method = c("logistic", "probit", "cloglog", "cauchit")) 
{
    logit <- function(p) log(p/(1 - p))
    fmin <- function(beta) {
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 
            100)
        eta <- offset
        if (pc > 0) 
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        if (all(pr > 0)) 
            -sum(wt * log(pr))
        else Inf
    }
    gmini <- function(beta) {
        jacobian <- function(theta) {
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0, k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
        }
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 
            100)
        eta <- offset
        if (pc > 0) 
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        p1 <- dfun(gamm[y + 1] - eta)
        p2 <- dfun(gamm[y] - eta)
        g1 <- if (pc > 0) 
            x * (wt * (p1 - p2)/pr)
        else numeric(0)
        xx <- .polrY1 * p1 - .polrY2 * p2
        g2 <- - xx * (wt/pr)
        g2 <- g2 %*% jacobian(theta)
        if (all(pr > 0)) 
            cbind(g1, g2)
        else NA+cbind(g1,g2)
    }
    gmin<-function(beta){
      colSums(gmini(beta))
    }
    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)

    pfun <- switch(method, logistic = plogis, probit = pnorm, 
        cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm, 
        cloglog = dgumbel, cauchit = dcauchy)


    m<-model.frame(formula,model.frame(design),na.action=na.pass)
    Terms <- attr(m, "terms")
    m<-na.action(m)
    nao<-attr(m,"na.action")
    if(length(nao)) {
      design<-design[-nao,]
    }

    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts")
    if (xint > 0) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1
    }
    else warning("an intercept is needed and assumed")

    wt <- weights(design)


    offset <- model.offset(m)
    if (length(offset) <= 1) 
        offset <- rep(0, n)
    y <- model.response(m)
    if (!is.factor(y)) 
        stop("response must be a factor")
    lev <- levels(y)
    if (length(lev) <= 2) 
        stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    if (missing(start)) {
        q1 <- length(lev)%/%2
        y1 <- (y > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <- switch(method, logistic = glm.fit(X, y1, wt/mean(wt), family = quasibinomial(), 
            offset = offset), probit = glm.fit(X, y1, wt/mean(wt), family = quasibinomial("probit"), 
            offset = offset), cloglog = glm.fit(X, y1, wt/mean(wt), family = quasibinomial("probit"), 
            offset = offset), cauchit = glm.fit(X, y1, wt/mean(wt), family = quasibinomial("cauchit"), 
            offset = offset))
        if (!fit$converged) 
            stop("attempt to find suitable starting values failed")
        coefs <- fit$coefficients
        if (any(is.na(coefs))) {
            warning("design appears to be rank-deficient, so dropping some coefs")
            keep <- names(coefs)[!is.na(coefs)]
            coefs <- coefs[keep]
            x <- x[, keep[-1], drop = FALSE]
            pc <- ncol(x)
        }
        spacing <- logit((1:q)/(q + 1))
        if (method != "logit") 
            spacing <- spacing/1.7
        gammas <- -coefs[1] + spacing - spacing[q1]
        thetas <- c(gammas[1], log(diff(gammas)))
        start <- c(coefs[-1], thetas)
    }
    else if (length(start) != pc + q) 
        stop("'start' is not of the correct length")
    res <- optim(start, fmin, gmin, method = "BFGS", hessian = TRUE, 
        ...)
    beta <- res$par[seq_len(pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1], exp(theta[-1])))
    deviance <- 2 * res$value/mean(wt)
    niter <- c(f.evals = res$counts[1], g.evals = res$counts[2])
    names(zeta) <- paste(lev[-length(lev)], lev[-1], sep = "|")
    if (pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    }
    else {
        eta <- rep(0, n)
    }
    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow = TRUE) - eta), 
        , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance, 
        fitted.values = fitted, lev = lev, terms = Terms, df.residual = sum(wt) - 
            pc - q, edf = pc + q, n = sum(wt), nobs = sum(wt), 
        method = method, convergence = res$convergence, 
        niter = niter)

    dn <- c(names(beta), names(zeta))
    H <- res$hessian
    dimnames(H) <- list(dn, dn)
    fit$Hessian <- H

    fit$call<-sys.call()
    fit$call[[1]]<-as.name(.Generic)
    
    inffun<- gmini(res$par)%*%solve(H)
    fit$var<-svyrecvar(inffun, design$cluster, 
                     design$strata, design$fpc,
                     postStrata = design$postStrata)
    fit$df.residual<-degf(design)-length(beta)


    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- "svyolr"
    fit
}


vcov.svyolr<-function(object,...) object$var

print.svyolr<-function (x, ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    if (length(coef(x))) {
        cat("\nCoefficients:\n")
        print(coef(x,intercept=FALSE), ...)
    }
    else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts:\n")
    print(x$zeta, ...)
    invisible(x)
  }


coef.svyolr<-function(object,intercepts=TRUE,...) {
	if(intercepts)
	    c(object$coefficients, object$zeta)
	else
	    object$coefficients
	}

summary.svyolr<-function (object, digits = max(3, .Options$digits - 3), correlation = FALSE, 
    ...) 
{
    cc <- coef(object)
    pc <- length(coef(object, FALSE))
    q <- length(object$zeta)
    coef <- matrix(0, pc + q, 3, dimnames = list(names(cc), c("Value", 
        "Std. Error", "t value")))
    coef[, 1] <- cc
    vc <- vcov(object)
    z.ind <- (pc + 1):(pc + q)
    gamma <- object$zeta
    theta <- c(gamma[1], log(diff(gamma)))
    jacobian <- function(theta) {
        k <- length(theta)
        etheta <- exp(theta)
        mat <- matrix(0, k, k)
        mat[, 1] <- rep(1, k)
        for (i in 2:k) mat[i:k, i] <- etheta[i]
        mat
    }
    J <- jacobian(theta)
    vc[z.ind, z.ind] <- J %*% vc[z.ind, z.ind] %*% t(J)
    coef[, 2] <- sd <- sqrt(diag(vc))
    coef[, 3] <- coef[, 1]/coef[, 2]
    object$coefficients <- coef
    object$pc <- pc
    object$digits <- digits
    if (correlation) 
        object$correlation <- (vc/sd)/rep(sd, rep(pc + q, pc + 
            q))
    class(object) <- "summary.svyolr"
    object
}

print.summary.svyolr<-function (x, digits = x$digits, ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    coef <- format(round(x$coefficients, digits = digits))
    pc <- x$pc
    if (pc > 0) {
        cat("\nCoefficients:\n")
        print(x$coefficients[seq_len(pc), , drop = FALSE], quote = FALSE, 
            ...)
    }
    else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts:\n")
    print(coef[(pc + 1):nrow(coef), , drop = FALSE], quote = FALSE, 
        ...)
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("(", mess, ")\n", sep = "")
    if (!is.null(correl <- x$correlation)) {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        print(correl[-1, -ncol(correl)], quote = FALSE, ...)
    }
    invisible(x)
}

model.frame.svyolr<-function(formula, ...){
	  mcall <- match.call(svyolr, formula$call)
	  design<- eval(mcall$design)
	  formula<-eval(mcall$formula)
	  mf<-model.frame(formula,model.frame(design))
	  w<-weights(design, type="sampling")
	  if (is.null(naa<-attr(mf,"na.action")))
	     mf[["(weights)"]]<-w
	  else
	     mf[["(weights)"]]<-w[-naa]
	  mf	
	}

