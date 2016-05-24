gqc <- function(formula,
	data,
	category,
	par = list(),
	zlimit = Inf,
	fixed = list(),
	opt = c("nlminb", "optim"),
	lower = -Inf,
	upper = Inf,
	control=list())
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    response <- model.response(mf)
    x <- model.matrix(Terms, mf)
    xint <- match("(Intercept)", colnames(x), nomatch=0L)
    if(xint > 0) x <- x[, -xint, drop=FALSE]
    dimen <- ncol(x)
    if(dimen < 2) stop("too few number of dimentions")
    if(missing(category)) category <- response
    if(nrow(data) != length(category)) 
        stop("nrow(data) and length(category) are different")
    #Set up the initial parameters, find one if missing
    if(missing(par)) par <- NULL
    if(is.null(par$pnoise)) par$pnoise <- 10
    if(is.null(par$cnoise)) par$cnoise <- 100
    if(is.null(par$coeff) | is.null(par$bias))
    {
    	mc <- mcovs.formula(formula, data=data, pooled=FALSE)
    	par <- qdb(means=mc$means, covs=mc$covs)
    }
    ncoeffs <- length(par$coeffs)
    if(missing(fixed)) fixed <- list(pnoise=FALSE, cnoise=FALSE, coeffs=FALSE, bias=FALSE)
    if(is.null(fixed$pnoise)) fixed$pnoise <- FALSE
    if(is.null(fixed$cnoise)) fixed$cnoise <- FALSE
    if(is.null(fixed$coeffs)) fixed$coeffs <- rep(FALSE, ncoeffs)
    if(length(fixed$coeffs) != ncoeffs)
        fixed$coeffs <- rep(fixed$coeffs, ncoeffs)[1:ncoeffs]
    if(is.null(fixed$bias)) fixed$bias <- FALSE
    fixidx <- unlist(fixed[c("pnoise","cnoise","coeffs","bias")])
    if(all(fixidx)) stop("no free parameters to fit.")
    initpar <- unlist(par)[!fixidx]
    
    if(missing(lower)){
        lower = c(.1, .1, rep(-Inf, length(unlist(par))-2))
        lower = lower[!fixidx]
    }
    if(missing(upper)){
	upper = c(5000, 5000, rep( Inf, length(unlist(par))-2))
	upper = upper[!fixidx]
    }
    .optfun <- function(par, response, x, zlimit, skeleton, fixed){
        params <- unlist(skeleton)
        params[!fixed] <- par
        params <- relist(params, skeleton=skeleton)
        -logLik(params, response, x, zlimit)
    }
    #Fit the model
    opt <- match.arg(opt)
    if(opt == "nlminb"){
        optRes <- nlminb(start=initpar, objective=.optfun,
            gradient = NULL, hessian = NULL,
            response = response, x = x, zlimit=zlimit, skeleton=par, fixed=fixidx,
            lower = lower, upper = upper, control = control)
        optRes$value <- optRes$objective
    } else {
        optRes <- optim(par=initpar, fn=.optfun,
            response = response, x = x, zlimit=zlimit, skeleton=par, fixed = fixidx,
            method = "L-BFGS-B", lower = lower, upper = upper, control = control)
    }
    if(optRes$convergence){
        msg <- paste(opt, " problem, convergence error code = ",
            optRes$convergence, "\n  message = ", optRes$message,
            sep='')
        warning(msg)
    }
    fit <- NULL
    fit$terms <- Terms
    cl <- match.call()
    cl[[1L]] <- as.name("gqc")
    fit$call <- cl
    fit$contrasts <- attr(x, "contrasts")
    fit$model <- mf
    fit$category <- category
    fit$initpar <- par
    #format the fitted parameters
    params <- unlist(par)
    params[!fixidx] <- optRes$par
    fit$par <- relist(params, skeleton=par)
    fit$logLik <- -optRes$value
    attr(fit$logLik, "df") <- length(optRes$par) - 1
    class(fit$logLik) <- "logLik"
    class(fit) <- "gqc"
    fit
}

gqcStruct <- function(pnoise, cnoise, coeffs, bias)
{
    par <- NULL
    par$pnoise <- pnoise
    par$cnoise <- cnoise
    par$coeffs <- coeffs
    par$bias <- bias
    class(par) <- c("gqcStruct", "list")
    par
}

print.gqc <- function(x, digits = max(5, getOption("digits") - 3), ...)
{
    cat('\nPerceptual Noise:\n');
    print(x$par$pnoise, digits=digits);

    cat('\nCriterial Noise:\n');
    print(x$par$cnoise, digits=digits);

    cat('\nBoundary Parameters:\n');
    print(c(x$par$coeff, x$par$bias), digits=digits);

    cat('\nNegative Loglikelihood:\n');
    cat(format(-x$logLik, digits=digits),
    " (df=",format(attr(x$logLik,"df")),")\n",sep="")

    cat('\nAIC score:\n');
    print(AIC(x), digits=digits)
}

logLik.gqcStruct <- function(object, response, x, zlimit = Inf, ...)
{
    x <- as.matrix(x)
    dimen <- ncol(x)
    nr <- nrow(x)
    if(nrow(x) != length(response)) stop("nrow(x) and length(response) are different")
    lev <- levels(g <- as.factor(response))
    ncoeff <- length(coeffs <- object$coeffs)
    idx <- 1:dimen
    if(dimen < 2 | dimen > 3) stop("Unsupported number of dimenensions.")
    if((sum(idx) + dimen) != ncoeff) stop("Inappropriate number of coeffs")

    A <- diag(coeffs[idx], nrow=dimen, ncol=dimen)
    A[lower.tri(A)] <- A[upper.tri(A)] <- (.5 * coeffs[(dimen+1):(ncoeff-dimen)])
    b <- coeffs[(ncoeff-dimen+1):(ncoeff)]
    c0 <- as.vector(object$bias)
    pnoise <- object$pnoise
    cnoise <- object$cnoise
    # Compute mean of h(x) for each x (Eqn.11, Ashby p.462).
    ## obtain index for product terms
    comb <- cbind(rbind(idx,idx),combn(idx,2))
    tmp1 <- pnoise*sum(diag(A))
    meanhxs <- tmp1 + cbind(x[,comb[1,]] * x[,comb[2,]], x, 1) %*% c(coeffs, c0)
    # Compute variance of h(x) for each x (Eqn.12, Ashby p.462)
    # Take advantage of fact that perceptual noise matrix is diagonal with
    # pnoise the same in all three dimenensions
    bvec <- matrix(rep(b,nr), ncol=nr, nrow=dimen, byrow=FALSE)
    varhxs <- 2*tmp1*tmp1 + pnoise*colSums((bvec + 2*A %*% t(x))^2)
    #Compute z-scores for each data point
    z <- drop( meanhxs / sqrt(varhxs+cnoise) )
    #Truncate the large z-scores
    z[z < -zlimit] <- -zlimit
    z[z >  zlimit] <-  zlimit
    log_A_probs <- pnorm(z[g == lev[1]], lower.tail=FALSE, log.p=TRUE)
    log_B_probs <- pnorm(z[g != lev[1]], lower.tail=TRUE, log.p=TRUE)
    res <- sum(log_A_probs, log_B_probs)
    attr(res, "df") <- length(unlist(object)) - 1
    class(res) <- "logLik"
    res
}

logLik.gqc <- function(object, ...)
{
    val <- object$logLik
    class(val) <- "logLik"
    val
}

extractAIC.gqc <- function(fit, scale, k = 2, ...)
{
    loglik <- fit$logLik
    edf <- attr(loglik, "df")
    c(edf, -2 * loglik + k * edf)
}