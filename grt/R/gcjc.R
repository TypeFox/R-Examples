gcjc <- function(formula,
    data,
    category,
    par,
    config = 1,
	zlimit = Inf,
	fixed = list(),
    equal.noise = TRUE,
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
    #biases <- as.vector(rowMeans(x))
    if(missing(par)) par <- list(noise=NULL,bias=NULL)
    if(!("gcjcStruct" %in% class(par)))
    {
        if(is.null(par$noise)) par$noise <- 10
        if(is.null(par$bias)){
            means <- sapply(split(as.data.frame(x), category), colMeans)
            par$bias <- as.vector(rowMeans(means))
        }
        par <- gcjcStruct(noise=par$noise,bias=par$bias, config)
        names(par) <- attr(Terms,"term.labels")
    }
    if(missing(fixed)) fixed <- list(noise=NULL, bias=NULL)
    if(is.null(fixed$noise)) fixed$noise <- c(FALSE, TRUE)
    if(is.null(fixed$bias)) fixed$bias <- c(FALSE, FALSE)
    if(length(fixed$noise) == 1) fixed$noise <- rep(fixed$noise, 2)
    if(equal.noise) fixed$noise[2] <- TRUE
    if(length(fixed$bias) == 1) fixed$bias <- rep(fixed$bias, 2)
    
    fixed <- list(
        list(noise=fixed$noise[1], coeffs=c(TRUE, TRUE), bias=fixed$bias[1]),
        list(noise=fixed$noise[2], coeffs=c(TRUE, TRUE), bias=fixed$bias[2])
        )
    names(fixed) <- names(par)
    fixidx <- unlist(fixed)
    new_par <- par
    initpar <- unlist(new_par)[!fixidx]

    lims <- apply(x,2,function(x){extendrange(range(x))})
    if(missing(lower)) lower <- unlist(gcjcStruct(noise=.001,bias=lims[1,],config))[!fixidx]
    if(missing(upper)) upper <- unlist(gcjcStruct(noise=500,bias=lims[2,],config))[!fixidx]
    #make sure that the lower and upper bouds are indeed lower and upper
    lims <- apply(rbind(lower,upper),2,range)
    lower <- lims[1,]
    upper <- lims[2,]
    
    .optfun <- function(par, response, x, zlimit, skeleton, fixed, equal.noise){
        ulpar <- unlist(skeleton)
        ulpar[!fixed] <- par
        params <- relist(ulpar, skeleton=skeleton)
        if(equal.noise)
            for(i in 1:length(params)){ params[[i]]$noise <- params[[1L]]$noise }
        class(params) <- c("gcjcStruct","list")
        -logLik.gcjcStruct(params,response,x,zlimit)
    }
    # Search for the parameters
    opt_method <- match.arg(opt)
    if(opt_method == "nlminb"){
        optRes <- nlminb(start=initpar, objective = .optfun,
            gradient = NULL, hessian = NULL,
            response=response, x=x, zlimit=zlimit, skeleton=new_par, fixed=fixidx, equal.noise=equal.noise,
            lower=lower, upper=upper, control=control)
        optRes$value <- optRes$objective
    } else {
        optRes <- optim(par=initpar, fn = .optfun,
            response=response, x=x, zlimit=zlimit, skeleton=new_par, fixed=fixidx, equal.noise=equal.noise,
            lower=lower, upper=upper, method="L-BFGS-B", control=control)
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
    cl[[1L]] <- as.name("goo")
    fit$call <- cl
    fit$contrasts <- attr(x, "contrasts")
    fit$model <- mf
    fit$category <- category
    fit$convergence <- optRes$convergence
    fit$initpar <- par
    #format the fitted parameters
    params <- unlist(par)
    params[!fixidx] <- optRes$par
    params <- relist(params, skeleton=par)
    if(equal.noise)
        for(i in 1:length(params)){ params[[i]]$noise <- params[[1L]]$noise }
    fit$par <- params
    fit$logLik <- -optRes$value
    #fit$RMSE <- RMSE.gcjcStruct(object=params, response=response, x=x, zlimit=zlimit)
    fit$convergence <- optRes$convergence
    attr(fit$logLik, "df") <- length(optRes$par)
    class(fit$logLik) <- "logLik"
    class(fit) <- "gcjc"
    fit
}

gcjcStruct <- function(noise, bias, config=c(1,2,3,4))
{
    if(length(noise) == 1)
        noise <- rep(noise, length.out=2)
    if(length(bias) != 2)
        stop("length(bias) must be 2")
        
    signs <- list(c(-1,1,-1,1), c(1,-1,-1,1), c(1,-1,1,-1), c(-1,1,1,-1))[[config]]
    par <- list(
        list(noise=noise[1], coeffs= signs[1]*c(1,0), bias=signs[2]*bias[1]),
        list(noise=noise[2], coeffs= signs[3]*c(0,1), bias=signs[4]*bias[2])
    )
    class(par) <- c("gcjcStruct", "list")
    par
}

print.gcjc <- function(x, digits = getOption("digits"), ...)
{
    cat('\nNoise:\n');
    print(sapply(x$par, "[[", "noise"))
    
    cat('\nBoundary Parameters :\n');
    print(unlist(coef(x)))
    
    cat('\nNegative Loglikelihood:\n');
    cat(format(-x$logLik, digits=digits),
    " (df=",format(attr(x$logLik,"df")),")\n",sep="")
    
    cat('\nAIC score:\n');
    print(AIC(x), digits=digits)
}

extractAIC.gcjc <- function(fit, scale, k = 2, ...)
{
    loglik <- fit$logLik
    edf <- attr(loglik, "df")
    c(edf, -2 * loglik + k * edf)
}

coef.gcjc <- function(object, ...)
{
  .coef <- function(obj)
  {
    val <- sum(obj$bias * obj$coeffs * -1)
    names(val) <- "(Intercept)"
    val
  }
  lapply(object$par, .coef)
}

logLik.gcjcStruct <- 
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
    #define a fuction to calculate probability
    .prob <- function(obj, x, zlimit)
    {
      a <- obj$coeffs
      if(is.null(obj$angle)) a <- obj$coeffs
      else a <- angle2cart(obj$angle)
      z_coefs <- c(a, obj$bias) / obj$noise
      z <- as.matrix(x) %*% as.matrix(z_coefs)
      # Truncate the large z-scores
      z[z < -zlimit] <- -zlimit
      z[z > zlimit] <- zlimit
      pnorm(z, lower.tail=FALSE, log.p=FALSE)
    }
    #do the rest
    pr <- sapply(object, .prob, x, zlimit)
    prA <- apply(pr, 1, prod)
    pr <- c( prA[g==lev[1]], 1 - prA[g!=lev[1]] )
    res <- sum(log(pr))
    attr(res, "df") <- dimen + 1
    class(res) <- "logLik"
    res
  }

logLik.gcjc <- function(object, ...)
{
  val <- object$logLik
  class(val) <- "logLik"
  val
}
