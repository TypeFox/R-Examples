glc <- function(formula,
	data,
	category,
	par = list(),
	zlimit = Inf,
	covstruct=c("unstructured", "scaledIdentity", "diagonal", "identity"),
	fixed = list(),
	opt = c("nlminb", "optim"),
	lower=-Inf,
	upper=Inf,
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
    if(missing(category)) category <- response
    if(nrow(data) != length(category)) 
        stop("nrow(data) and length(category) are different")
    if(missing(par)) par <- NULL
    if(is.null(par$noise)) par$noise <- 10
    if(is.null(par$coeffs) | is.null(par$bias))
    {
        mcovs <- mcovs.default(x = x, grouping = category, pooled = TRUE)
        covstruct <- match.arg(covstruct)
    	par <- ldb(means=mcovs$means, covs=mcovs$covs, covstruct=covstruct, noise=par$noise)
    }
    if(missing(fixed)) fixed <- list(noise=FALSE, coeffs=FALSE, bias=FALSE)
    if(is.null(fixed$noise)) fixed$noise <- FALSE
    if(is.null(fixed$coeffs)) fixed$coeffs <- FALSE
    if(is.null(fixed$bias)) fixed$bias <- FALSE
    if(dimen == 1) {
        fixed$coeffs <- TRUE
        fixidx <- unlist(fixed[c("noise","coeffs","bias")])
        new_par <- glcStruct(noise=par$noise, coeffs=par$coeffs, bias=par$bias)
        lim <- extendrange(range(-par$coeffs*x))
        if(missing(lower)) lower <- c(.001, 0, lim[1])[!fixidx]
        if(missing(upper)) upper <- c(500, 0, lim[2])[!fixidx]
    } else {
        fixed$angle <- fixed$coeffs
        fixidx <- unlist(fixed[c("noise","angle","bias")])
        # Convert to search format: list(noise, angle, bias)
        new_par <- old2new_par(par)
        if(missing(lower)) lower <- c(.001, -Inf, -Inf)[!fixidx]
        if(missing(upper)) upper <- c(500, Inf, Inf)[!fixidx]
    }
    if(all(fixidx)) stop("no free parameters to fit.")
    initpar <- unlist(new_par[!fixidx])
    # Set up the optimization function
    .optfun <- function(par, response, x, zlimit, skeleton, fixed){
        freepar <- relist(par, skeleton=skeleton[!fixed])
        fixedpar <- skeleton[fixed]
        params <- c(freepar, fixedpar)[names(skeleton)]
        class(params) <- c("glcStruct","list")
        -logLik(params, response, x, zlimit)
    }
    # Search for the parameters
    opt_method <- match.arg(opt)
    if(opt_method == "nlminb"){
        optRes <- nlminb(start=initpar, objective = .optfun,
            gradient = NULL, hessian = NULL,
            response=response, x=x, zlimit=zlimit, skeleton=new_par, fixed=fixidx,
            lower=lower, upper=upper, control=control)
        optRes$value <- optRes$objective
    } else {
        optRes <- optim(par=initpar, fn = .optfun,
            response=response, x=x, zlimit=zlimit, skeleton=new_par, fixed=fixidx,
            method="L-BFGS-B", lower=lower,upper=upper, control=control)
    }
    if(optRes$convergence){
        msg <- paste(opt_method, " problem, convergence error code = ",
            optRes$convergence, "\n  message = ", optRes$message,
            sep='')
        warning(msg)
     }
    # Format the output
    fit <- NULL
    fit$terms <- Terms
    cl <- match.call()
    cl[[1L]] <- as.name("glc")
    fit$call <- cl
    fit$contrasts <- attr(x, "contrasts")
    fit$model <- mf
    fit$category <- category
    fit$convergence <- optRes$convergence
    fit$initpar <- par
    #format the fitted parameters
    freepar <- relist(optRes$par, skeleton=new_par[!fixidx])
    fixedpar <- new_par[fixidx]
    params <- c(freepar, fixedpar)
    params <- unlist(params[names(new_par)])
    if(dimen == 1) {
        fit$par <- relist(params, skeleton=new_par)
    } else {
        fit$par <- new2old_par(relist(params, skeleton=new_par))
    }
    fit$logLik <- -optRes$value
    #fit$RMSE <- RMSE.glcStruct(object = fit$par, response=response, x=x, zlimit=zlimit)
    attr(fit$logLik, "df") <- length(optRes$par)
    class(fit$logLik) <- "logLik"
    class(fit) <- "glc"
    fit
}

glcStruct <- function(noise, coeffs, bias)
{
    par <- NULL
    #denom <- ifelse(length(coeffs) > 1, max(svd(coeffs)$d), 1)
    denom <- ifelse(length(coeffs) > 1, sqrt(sum(coeffs^2)), 1)
    par$noise <- noise / denom
    par$coeffs <- coeffs / denom
    par$bias <- bias / denom
    class(par) <- c("glcStruct", "list")
    par
}

print.glc <- function(x, digits = getOption("digits"), ...)
{
    cat('\nNoise:\n');
    print(x$par$noise, digits=digits);

    cat('\nBoundary Parameters:\n');
    print(c(x$par$coeff, x$par$bias), digits=digits);

    cat('\nNegative Loglikelihood:\n');
    cat(format(-x$logLik, digits=digits),
    " (df=",format(attr(x$logLik,"df")),")\n",sep="")
    
    cat('\nAIC score:\n');
    print(AIC(x), digits=digits)
}

predict.glc <- function(object, newdata, seed = NULL, ...)
{
    if(!inherits(object, "glc")) stop("x not of class \"glc\"")
    if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)                     # initialize the RNG if necessary
    if(is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
	set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    noise <- object$par$noise
    coeffs <- object$par$coeffs
    bias <- object$par$bias
    label <- levels(factor(object$category))
    if(missing(newdata)){
        newdata <- model.matrix(delete.response(object$terms), object$model)
        xint <- match("(Intercept)", colnames(newdata), nomatch=0L)
        if(xint > 0) newdata <- newdata[, -xint, drop=FALSE]
    }
    newdata <- as.matrix(newdata)
    if(any(!is.finite(newdata))) 
        stop("infinite, NA or NaN values in 'newdata'")
    n <- nrow(newdata)
    if(ncol(newdata) != length(attr(object$terms,"term.labels")))
        stop("mismatching ncol(newdata)")
    snoise <- rnorm(n=n, mean=0, sd=noise)
    
    discrimvals <- newdata %*% coeffs + bias + snoise
    response <- vector(mode="numeric",length=n)
    response[discrimvals <= 0] <- label[1]
    response[discrimvals >  0] <- label[2]
    response
}

coef.glcStruct <- function(object, ...)
{
    if (!inherits(object, "glcStruct")) 
        stop("object not of class \"glcStruct\"")
    b <- object$coeffs
    c0 <- object$bias
    if(length(b) == 1){
        coeffs <- -c0/b
        names(coeffs) <- "(Intercept)"
    } else {
        if(is.null(b)) b <- angle2cart(object$angle)
        dimen <- length(b)
        if(!is.null(names(b))) xlabels <- names(b)
        else if(dimen < 4) xlabels <- c("x","y","z")[1:(dimen)]
        else xlabels <- paste("x", 1:(dimen), sep="")
        coeffs <- -c(c0, b[-dimen]) / b[dimen]
        names(coeffs) <- c("(Intercept)", xlabels[-dimen])
    }
    coeffs
}

coef.glc <- function(object, ...)
{
    if (!inherits(object, "glc")) stop("x not of class \"glc\"")
    if(is.null(names(object$par$coeffs))){
        Terms <- object$terms
        xlabels <- attr(Terms,"term.labels")
        names(object$par$coeffs) <- xlabels
    }
    coef(object$par)
}

logLik.glcStruct <- 
    function(object,
        response, 
        x, 
        zlimit = Inf, ...)
{
    dimen <- ncol(as.matrix(x))
    x <- as.matrix(cbind(x,1))
    if(nrow(x) != length(response)) 
        stop("nrow(x) and length(response) are different")
    g <- as.factor(response)
    lev <- levels(g)
    
    if(is.null(object$angle)) a <- object$coeffs
    else a <- angle2cart(object$angle)
    z_coefs <- c(a, object$bias) / object$noise
    z <- x %*% as.matrix(z_coefs)
    # Truncate the large z-scores
    z[z < -zlimit] <- -zlimit
    z[z > zlimit] <- zlimit
    log_A_probs <- pnorm(z[g==lev[1]], lower.tail=FALSE, log.p=TRUE)
    log_B_probs <- pnorm(z[g!=lev[1]], lower.tail=TRUE, log.p=TRUE)
    res <- sum(log_A_probs, log_B_probs)
    attr(res, "df") <- dimen + 1
    class(res) <- "logLik"
    res
}

logLik.glc <- function(object, ...)
{
    val <- object$logLik
    class(val) <- "logLik"
    val
}

extractAIC.glc <- function(fit, scale, k = 2, ...)
{
    loglik <- fit$logLik
    edf <- attr(loglik, "df")
    c(edf, -2 * loglik + k * edf)
}

# RMSE.glcStruct <- 
#     function(object,
#         response, 
#         x, 
#         zlimit = Inf, ...)
# {
#     dimen <- ncol(as.matrix(x))
#     x <- as.matrix(cbind(x,1))
#     if(nrow(x) != length(response)) 
#         stop("nrow(x) and length(response) are different")
#     g <- as.factor(response)
#     lev <- levels(g)
#     
#     if(is.null(object$angle)) a <- object$coeffs
#     else a <- angle2cart(object$angle)
#     z_coefs <- c(a, object$bias) / object$noise
#     z <- x %*% as.matrix(z_coefs)
#     # Truncate the large z-scores
#     z[z < -zlimit] <- -zlimit
#     z[z > zlimit] <- zlimit
#     resp <- rep(0, length.out=nrow(x))
#     resp[g==lev[1]] <- 1
#     #calculate p(A|item)
#     prob <- pnorm(z, lower.tail=FALSE)
#     sqrt(mean((resp - prob)^2))
# }