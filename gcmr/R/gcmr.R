# Return a function which computes an approximation of log(c(p(y[1]),p(y[2]|y[1]),...))
llik <- function(m, irep, magic=NA) {
    n <- m$n
    not.na <- m$not.na
    y <- m$y[not.na,]
    x <- m$x[not.na,,drop=FALSE]
    z <- m$z[not.na,,drop=FALSE]
    offset <- list(mean = m$offset$mean[not.na,,drop=FALSE],
                   precision = m$offset$precision[not.na,,drop=FALSE])
    is.int <- m$marginal$type == "integer"
    ind.lik <- !is.null(m$cormat$independent) 
    dp <- m$marginal$dp
    ibeta <- m$ibeta
    igamma <- m$igamma
    vchol <- m$cormat$chol
    if (missing(irep)) irep <- max(m$options$nrep)
    seed <- m$options$seed
    id <- if(is.null(m$cormat$id)) rep(1,n) else m$cormat$id
    lstrata <- rle(id)$length
    ipar <- c(if(is.int) 1 else 0,1,n,irep,length(lstrata),lstrata)
    magic <- rep(magic,m$n)
    theta <- m$fixed
    ifree <- is.na(theta)
    cache <- new.env()
    function( theta.free ) {
        theta[ifree] <- theta.free
        beta <- theta[ibeta]
        if (!identical(cache$beta,beta)) {
            dp <- dp(y,x,z,offset,beta)
            assign("beta",beta,envir=cache)
            assign("dp",dp,envir=cache)
        } else {
            dp <- get("dp",envir=cache)
        }
        if ( is.null(dp) ) return( magic )
        if ( ind.lik ) return( log(dp[,1]) )
        gamma <- theta[igamma]
        if (!identical(cache$gamma,gamma)) {
            q <- vchol( gamma , not.na )
            assign("gamma",gamma,envir=cache)
            assign("q",q,envir=cache)
        } else {
            q <- get("q",envir=cache)
        }
        if ( is.null(q) ) return ( magic )
        if (is.int) set.seed(seed)
        lk <- .Call(gcmrcomp, ipar, unlist(q), dp)
        if ( all(is.finite(lk)) ) lk else magic
    }    
} 

# Workhorse which is called from gcmr.fit and profile.gcmr
truefit <- function(x) {
    # saving/restoring the random seed
    if ( x$marginal$type == "integer" ) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
    }
    ifree <- is.na( x$fixed )
    theta <- x$fixed
    theta[ifree] <- x$estimate[ifree]
    low <- x$lower[ifree]
    up <- x$upper[ifree]
    big <- -sqrt(.Machine$double.xmax)
    nrep <- if(x$marginal$type=="numeric") 1 else x$options$nrep
    for ( i in 1:length(nrep) ) {
      log.lik <- llik(x,nrep[i],big)
      if( length( theta[ifree] ) > 0 ){
          ans <- suppressWarnings( x$options$opt(theta[ifree] , log.lik , low, up) )
        theta[ifree] <- ans$estimate
      }
    }
    names(theta) <- names(x$estimate)
    x$estimate <- theta
    x$maximum <- if( length( theta[ifree] ) > 0 ) ans$maximum else -sum( log.lik( theta ) )
    x$convergence <- ans$convergence
    x
}

# Add/replace an estimate of jacobian and hessian to an gcmr object
# Hessian is approximated by finite differences
# in a rotated space in which the true hessian is near
# the identity matrix. The rotation is based on
# jacobian crossprod
jhess <- function(x, options=x$options, only.jac = FALSE) {
    if ( !inherits( x , "gcmr" ) ) stop("First argument must be a gcmr object")
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    theta <- x$estimate
    ifree <- is.na( x$fixed )
    theta.free <- theta[ifree]
    low <- x$lower[ifree]
    up <- x$upper[ifree]
    xlik <- llik(x)
    log.lik <- function(th) xlik(pmax(low,pmin(up,th)))
    eps <- .Machine$double.eps^(1/4)
    relStep <- 0.1
    maxtry <- 10
    delta <- ifelse(abs(theta.free)<1, eps, eps*theta.free)
    di <- function(i,delta) {
        x1 <- x2 <- theta.free
        x1[i] <- x1[i] - delta[i]
        x2[i] <- x2[i] + delta[i]
        (log.lik(x2)-log.lik(x1))/(2*delta[i])
    }
    while (1) {
        x$jac <- sapply(seq_along(theta.free),di,delta)
        if( all(is.finite(x$jac)) ) break
        delta <- delta/2
        maxtry <- maxtry - 1
        if (maxtry<0) stop("impossible to compute a finite jacobian")
    }
    if (!only.jac) {
        a <- svd(x$jac)
        a$d <- pmax(a$d,sqrt(.Machine$double.eps)*a$d[1])
        x$hessian <- nlme::fdHess(rep(0,length(theta.free)),
                                  function(tx) sum(log.lik(theta.free+a$v%*%(tx/a$d))),
                                  minAbsPar=1,.relStep=relStep)$Hessian
        x$hessian = (x$hessian+t(x$hessian))/2
        x$hessian <- a$v%*%(outer(a$d,a$d)*(x$hessian))%*%t(a$v)
    }
    x$options$no.se <- FALSE
    x
}

gradcheck <- function(x, options=x$options) {
    j <- jhess(x,options,only.jac=TRUE)$jac
    g <- colSums(j)
    sum(g*solve(crossprod(j),g))
}

gcmr.opt <- function(start, loglik, lower, upper) {
    fn.opt <- function(x){
        if( any(x <= lower || x >= upper) ) NA
        else  -sum(loglik(x))
    }
    ans <- optim(start, fn.opt, method="BFGS")
    if(ans$convergence) warning(paste("optim exits with code",ans$convergence))
    list(estimate=ans$par,maximum=ans$value,convergence=ans$convergence)
}

gcmr.options <- function(seed=round(runif(1,1,100000)), nrep=c(100,1000),
                         no.se=FALSE, method=c("BFGS", "Nelder-Mead", "CG"), ...) {

    method <- match.arg(method)
    control <- list(...)
    opt <- function(start, loglik, lower, upper) {
        fn.opt <- function(x){
            if( any(x <= lower || x >= upper) ) NA
            else  -sum(loglik(x))
        }
        ans <- optim(start, fn.opt, method=method, control=control)
        if(ans$convergence) warning(paste("optim exits with code",ans$convergence))
    list(estimate=ans$par,maximum=ans$value,convergence=ans$convergence)
    }
    list(seed=seed,nrep=nrep,no.se=no.se,opt=opt,method=method)
}

gcmr <- function(formula, data, subset, offset, marginal, cormat,
                 start, fixed, options=gcmr.options(...), model=TRUE, ...){

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if (missing(data)) data <- environment(formula)
    marginal <- eval(cl$marginal,data,parent.frame()) 
    cormat <- eval(cl$cormat,if (missing(subset)) data else data[eval(mf$subset,data),], parent.frame()) 
    if (is.function(marginal))
        marginal <- marginal()
    if (!inherits(marginal, "marginal.gcmr")) stop("Unknown marginal")
    if (is.function(cormat))
        cormat <- cormat()
    if (!inherits(cormat, "cormat.gcmr")) stop("Unknown cormat")
    m <- match(c("formula", "data", "subset", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass
    ############################################################
    ## Next lines "inherited" from betareg (Cribari-Neto and Zeleis, 2010, JSS)
    ## with minor modifications needed for gcmr
    ## For complete reference to betareg see citation("betareg")
    ############################################################
    oformula <- as.formula(formula)
    formula <- as.Formula(formula)
    if(length(formula)[2L] < 2L) {
        formula <- as.Formula(formula(formula), ~ 1)
        simple_formula <- TRUE
    } else {
        if(length(formula)[2L] > 2L) {
            formula <- Formula(formula(formula, rhs = 1:2))
            warning("formula must not have more than two RHS parts")
        }
        simple_formula <- FALSE
    }
    mf$formula <- formula
    
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    mt <- terms(formula, data = data)
    mtX <- terms(formula, data = data, rhs = 1L)
    mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
    Y <- model.response(mf, "any")
    ## as in glm
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- model.matrix(mtX, mf)
    Z <- model.matrix(mtZ, mf)

    expand_offset <- function(offset) {
        if (is.null(offset)) offset <- 0
        if (length(offset) == 1) offset <- rep.int(offset, NROW(Y))
        as.vector(offset)
    }
    offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
    offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
    if (!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
    offset <- list(mean = offsetX, precision = offsetZ)
    
    fit <- gcmr.fit(X, Y, if(!simple_formula) Z, offset, marginal, cormat, start, fixed, options)
    
    fit$call <- cl
    fit$formula <- oformula
    fit$terms <- list(mean = mtX, precision = mtZ, full = mt)
    fit$levels <- list(mean = .getXlevels(mtX, mf), precision = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
    fit$contrasts <- list(mean = attr(X, "contrasts"), precision = attr(Z, "contrasts"))
    if(model) fit$model <- mf

    class(fit) <- "gcmr"
    return(fit)
}

gcmr.fit <- function(x=rep(1,NROW(y)), y, z=NULL, offset=NULL, marginal, cormat, start, fixed, options=gcmr.options()) {
  if( is.function( marginal ) )
    marginal <- marginal()
  if( is.function( cormat ) )
    cormat <- cormat()
  ## arguments check
  if ( !inherits( marginal , "marginal.gcmr" ) ) stop("Unknown marginal")
  if ( !inherits( cormat , "cormat.gcmr" ) ) stop("Unknown cormat")
  ## gcmr object
  nb <- marginal$npar(x, z)
  ng <- cormat$npar
  if (is.null(offset))
    offset <- rep.int(0, NROW(y))
  if (!is.list(offset))
    offset <- list(mean = offset, precision = rep.int(0, NROW(y)))
  else{
    if (is.null(offset$mean))
      offset$mean <- rep.int(0, NROW(y))
    if (is.null(offset$precision))
      offset$precision <- rep.int(0, NROW(y))
  }
  not.na <- apply(cbind(y,x,z,offset$mean,offset$precision),1,function(x) !any(is.na(x)))
  m <- structure(list(y=as.matrix(y), x=as.matrix(x), z=if(!is.null(z)) as.matrix(z) else NULL,
                      offset=list(mean=as.matrix(offset$mean), precision=as.matrix(offset$precision)),
                      not.na=not.na, n=sum(not.na), marginal=marginal, cormat=cormat, ibeta=1:nb,
                      igamma=if (ng) (nb+1):(nb+ng) else NULL, nbeta=nb, ngamma=ng, call=match.call()), class="gcmr")
  if ( missing(start) ) {
    lambda <- marginal$start(m$y[not.na,],m$x[not.na,,drop=FALSE],m$z[not.na,,drop=FALSE],
                             list(mean = m$offset$mean[not.na,], precision = m$offset$precision[not.na,]))
    tau <- cormat$start()
    m$estimate <- c(lambda,tau)
    ll <- attr(lambda,"lower")
    lt <- attr(tau,"lower")
    m$lower <- c(if(is.null(ll)) rep(-Inf,length(lambda)) else ll,
                 if(is.null(lt)) rep(-Inf,length(tau)) else lt)
    ll <- attr(lambda,"upper")
    lt <- attr(tau,"upper")
    m$upper <- c(if(is.null(ll)) rep(Inf,length(lambda)) else ll,
                 if(is.null(lt)) rep(Inf,length(tau)) else lt)
  } else {
    m$estimate <- start
    ll <- attr(start,"lower")
    m$lower <- if(is.null(ll)) rep(-Inf,length(start)) else ll
    ll <- attr(start,"upper")
    m$upper <- if(is.null(ll)) rep(Inf,length(start)) else ll
  }
  if (length(m$estimate) != nb+ng) stop("mismatch in the number of initial parameters")
  if (length(m$lower) != nb+ng) stop("length of lower different from the number of the parameters")
  if (length(m$upper) != nb+ng) stop("length of upper different from the number of the parameters")
  if ( missing(fixed) ) {
    m$fixed <- rep( NA , length(m$estimate) )
  } else {
    if (length(fixed) != length(m$estimate) ) stop("fixed has a wrong length")
    m$fixed <- fixed
  }
  ifree <- is.na( m$fixed )
  if( length( m$estimate[ifree] ) == 0 ) ## if all parameters are fixed, skip se
    options$no.se <- TRUE
  m$options <- options
  ## compute estimate
  m <- truefit(m)
  fv <- marginal$fitted.val(x=m$x[not.na,,drop=FALSE],
                            z=m$z[not.na,,drop=FALSE],
                            offset=list(mean = m$offset$mean[not.na,],
                                precision = m$offset$precision[not.na,]),
                            lambda=m$estimate[1:NCOL(x)])
  m$fitted.values <- rep(NA, NROW(y))
  m$fitted.values[not.na] <- fv
  ## and s.e. and return
  if (m$options$no.se) m else jhess(m)
}

coef.gcmr <- function(object,...) object$estimate

logLik.gcmr <- function(object,...) {
    ans <- -object$maximum
    attr(ans,"df") <- sum( is.na(object$fixed) )
    class(ans) <- "logLik"
    ans
}

terms.gcmr <- function(x, model = c("mean", "precision"), ...) {
  x$terms[[match.arg(model)]]
}

model.frame.gcmr <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  if(is.Formula(formula$formula)) formula$call$formula <- formula$formula <-
    formula(formula$formula, collapse = TRUE)
  formula$terms <- formula$terms$full
  NextMethod()
}

## almost as in betareg
model.matrix.gcmr <- function(object, model = c("mean", "precision"), ...) {
    model <- match.arg(model)
    model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
}

## exactly as in betareg
update.gcmr <- function (object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if(is.null(call)) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) call$formula <- formula(update(Formula(formula(object)), formula.))
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}

pinv <- function(h,tol) {
    h <- svd(h)
    idx <- h$d > sqrt(.Machine$double.eps)*h$d[1]
    ans <- h$u[,idx,drop=FALSE]%*%( (1/h$d[idx])*t(h$u[,idx,drop=FALSE]))
    attr(ans,"df") <- sum(idx)
    ans
}

estfun.gcmr <- function(x, ...) {
    if (is.null(x$cormat$id)) x$jac
    else t(sapply(split(1:x$n, x$cormat$id[x$not.na]), function(i) colMeans(x$jac[i,])))
}

bread.gcmr <- function(x,...) pinv(x$hessian)*x$n

coeftest.gcmr <- function(x, vcov. = NULL, df = Inf, ...)
    coeftest.default(x, vcov. = vcov., df = df, ...)

vcov.gcmr <- function(object, ...) {
    if ( !inherits( object , "gcmr" ) ) stop("First argument must be a gcmr object")
    if (is.null(object$hessian) || is.null(object$jac)) object <- jhess(object)
    v <- pinv(object$hessian)
    ans <- matrix(0,length(object$estimate),length(object$estimate))
    ifree <- is.na(object$fixed)
    ans[ifree,ifree] <- v
    colnames(ans) <- rownames(ans) <- names(object$estimate)
    ans
}

print.gcmr <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(x$convergence) {
      cat("model did not converge\n")
  } else {
      if( x$nbeta ) {
          cat("Marginal model parameters:\n")
          print.default(format(x$estimate[x$ibeta], digits = digits), print.gap = 2, quote = FALSE)
          cat("\n")
      } else cat("No parameters in the marginal model\n\n")
    if( x$ngamma ){
        cat("Gaussian copula parameters:\n")
        print.default(format(x$estimate[x$igamma] , digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
    } else cat("No parameters in the Gaussian copula\n\n")
  }
  invisible(x)
}

## summary
summary.gcmr <- function(object, ...)
{

    cf <- object$estimate
    se <- sqrt(diag(vcov(object)))
    cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
    ## set to NA se, z and p-value for parameters with fixed values
    cf[which(!is.na(object$fixed)), -1] <- NA
    colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    cf <- list(marginal = cf[seq.int(length.out = object$nbeta), , drop = FALSE], 
               copula = cf[seq.int(length.out = object$ngamma) +  object$nbeta, , drop = FALSE])
    rownames(cf$marginal) <- names(object$estimate)[seq.int(length.out = object$nbeta)]
    rownames(cf$copula) <- names(object$estimate)[seq.int(length.out = object$ngamma) + object$nbeta]
    object$coefficients <- cf
    object$aic <- AIC(object)

    ## delete some slots
     object$fitted.values <- object$y <- object$x <- object$z <- object$offset <- object$levels <- object$ibeta <- object$igamma <- object$fixed <- object$hessian <- object$jac <- object$marginal <- object$cormat <- object$not.na <- object$options <- NULL

    class(object) <- "summary.gcmr"
    object
}

print.summary.gcmr <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    if(x$convergence) {
        cat("model did not converge\n")
    } else {
        if( x$nbeta ) {
            cat("\nCoefficients marginal model:\n")
            printCoefmat(x$coefficients$marginal, digits = digits, signif.legend = FALSE)
        } else cat("\nNo coefficients in the marginal model\n")

        if( x$ngamma ) {
            cat("\nCoefficients Gaussian copula:\n")
            printCoefmat(x$coefficients$copula, digits = digits, signif.legend = FALSE)
        } else cat("\nNo coefficients in the Gaussian copula\n")
        
        if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
            cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

        cat("\nlog likelihood = ", formatC(x$maximum, digits = max(5L, digits + 1L)),
        ",  AIC = ", format(x$aic, digits = max(4L, digits + 1L)), "\n", sep = "")
    }
    invisible(x)
}


## note: method=interval is hidden
residuals.gcmr <- function (object, type=c("conditional","marginal"),
                            method=c("random","mid"), ...) {
    type <- match.arg(type)
    method <- match.arg(method)
    is.int <- object$marginal$type == "integer"
    cond <- is.null(object$cormat$independent) && (type=="conditional") 
    if (is.int) {
        ##before saving/setting seed to obtain different randomization 
        if (method=="random") u <- runif(object$n) 
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
        set.seed(object$options$seed)
    }
    a <- object$marginal$dp(object$y[object$not.na,],
                            object$x[object$not.na,,drop=FALSE],
                            object$z[object$not.na,,drop=FALSE],
                            list(mean = object$offset$mean[object$not.na,],
                                 precision = object$offset$precision[object$not.na,]),
                            object$estimate[object$ibeta])
    if (cond) {
        q <- object$cormat$chol(object$estimate[object$igamma],object$not.na)
        id <- if(is.null(object$cormat$id)) rep(1,object$n) else object$cormat$id
        lstrata <- rle(id)$length
        ipar <- c(if(is.int) 1 else 0, 0,
                  object$n,max(object$options$nrep),length(lstrata),lstrata)
        a <- .Call(gcmrcomp,ipar,unlist(q),a)
    }
    if (is.int) {
        if (method=="random") {
            ans <- rep(NA,length(object$not.na))
            ans[object$not.na] <- qnorm(a[,2]-u*a[,1])            
        } else if (method=="interval") {
            ans <- matrix(NA,length(object$not.na),2)
            ans[object$not.na,] <- qnorm(cbind(a[,2],a[,2]-a[,1]))
        } else {
            ans <- rep(NA,length(object$not.na))
            ans[object$not.na] <- qnorm(a[,2]-a[,1]/2)            
        }
    } else {
            ans <- rep(NA,length(object$not.na))
            ans[object$not.na] <- if (cond) a[,2] else qnorm(a[,2])            
    }
    ans
}


profile.gcmr <- function(fitted , which , low , up, npoints = 10 , display = TRUE , alpha = 0.05, progress.bar = TRUE, ... ) {
    if ( !inherits( fitted , "gcmr" ) ) stop("first argument must be a gcmr object")
    if ( !( which %in% fitted$ibeta ) )
        stop("profile likelihood is computed only for mean or dispersion parameters")
    if(missing(low) || missing(up)) se <- sqrt(diag(vcov(fitted)))
    if(missing(low)) low <- coef(fitted)[which]-3*se[which]
    if(missing(up)) up <- coef(fitted)[which]+3*se[which]
    low <- max(low, fitted$lower[which])
    up <- min(up, fitted$upper[which])
    points <- seq( low , up , length = npoints )
    ## progress bar... 
    pb <- txtProgressBar(min=0, max=npoints)  
    prof <- function(x,which,fitted) {
        fitted$fixed[which] <- x
        truefit(fitted)$maximum
    }
    loglik <- rep(NA, npoints)
    for(i in 1:npoints){
        loglik[i] <- prof(points[i], which, fitted)
        if(progress.bar) setTxtProgressBar(pb, i) ## progress bar updated
    }
    close(pb)  ## progress bar closed! 
    points <- c( points , fitted$par[which] )
    ord <- order( points)
    points <- points[ord]
    loglik <- -c(loglik,fitted$maximum)[ord]
    if ( display ) {
        npoints <- seq( low , up , length = 200 )
        plot( npoints , splinefun( points, loglik )(npoints) , type = "l" ,
             xlab = names(fitted$estimate)[which], ylab = "log-likelihood profile", ...)
        grid()
        abline( h = max(loglik) - qchisq( 1 - alpha , 1 )/2 , lty = "dashed" )
    }
    invisible(list(points=points,profile=loglik))
}

plot.gcmr <- function(x, which = if(!time.series) 1:4 else c(1, 3, 5, 6), caption = c("Residuals vs indices of obs.", "Residuals vs linear predictor", "Normal plot of residuals", "Predicted vs observed values", "Autocorrelation plot of residuals", "Partial ACF plot of residuals"), main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), level = 0.95, col.lines = "gray", time.series = inherits(x$cormat, "arma.gcmr"), ...){

    if (!inherits(x, "gcmr")) 
        stop("use only with \"gcmr\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 7)) 
        stop("'which' must be in 1:6")
    res <- residuals(x)
    fitted.val <- x$fitted.values
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    Main <- rep("", 6)
    Main[which] <- rep(main, length.out = sum(show))
    one.fig <- prod(par("mfcol")) == 1
    op <- par(ask = ask, ...)
    on.exit(par(op))
    if (show[1]) {
        plot(1:length(res), res, xlab = "Obs. number", ylab = "Quantile residual", main = Main[1], type = if(time.series) "l" else "p", ...)
        abline(h = 0, col = col.lines, ...)
        mtext(caption[1], 3, 0.25)
    }
    if (show[2]) {
        plot(fitted(x), res, xlab = "Linear predictor", ylab = "Quantile residual", main = Main[2], ...)
        abline(h = 0, col = col.lines,...)
        mtext(caption[2], 3, 0.25)
    }
    if (show[3]) {
        qqPlot(res, envelope = level, grid = FALSE, xlab = "Normal quantiles", ylab = "Sorted quantile residuals", main = Main[3], col.lines = col.lines)
        mtext(caption[3], 3, 0.25)
    }
    if (show[4]) {
        Y <- if(NCOL(x$y)==1) x$y else x$y[,1]/(x$y[,1]+x$y[,2])
        plot(Y, fitted.val, xlab = "Observed values", ylab = "Predicted values", main = Main[4], ...)
        abline(0, 1, col = col.lines, ...)
        mtext(caption[4], 3, 0.25)
    }
    if (show[5]) {
        plot( acf(res, na.action = na.pass, plot = FALSE), ci.col = col.lines, main = Main[5], ...)
        mtext(caption[5], 3, 0.25)
    }
    if (show[6]) {
        plot( pacf(res, na.action = na.pass, plot = FALSE), ci.col = col.lines, main = Main[6], ...)
        mtext(caption[6], 3, 0.25)
    }
    invisible()
}
 

