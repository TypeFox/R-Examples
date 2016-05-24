"step.length" <- function(corrector, lambda2, min.lambda, max.arclength, frac.arclength, add.newvars, backshoot, h0=NULL, eps=.Machine$double.eps)
  {
    active <- corrector$active
    force.active <- corrector$force.active
    lambda <- corrector$lambda - min.lambda
    Lambda2 <- c(0, rep(lambda2, length(active)-1))
    C <- corrector$corr
    b <- corrector$b[active]
    xw <- corrector$xw
    xwx <- t(xw) %*% xw[ ,active]
    xwx.in <- solve(xwx[active, ] + diag(Lambda2))
    db <- drop(xwx.in %*% c(rep(0, length(force.active)), sign(C[active[-force.active]])))
    newa <- NULL
    if (!backshoot) {
      Cmax <- max(abs(C))
      inactive <- seq(C)[-active]
      ninact <- length(inactive)
      a <- drop(xwx[inactive, ] %*% db)
      gam <- c((Cmax - C[inactive])/(1 - a), (Cmax + C[inactive])/(1 + a))
      ha <- min(gam[gam > eps], lambda)
      hd <- -b[-force.active]/db[-force.active]
      h <- min(hd[hd > eps], ha)
      if (h==ha & h < lambda) {
        ii <- which(gam > eps)
        nii <- length(ii)
        oo <- order(gam[ii])
        oo <- oo[1:min(nii, add.newvars)]
        oo <- ii[oo]
        oo <- unique(ifelse(oo <= ninact, oo, oo-ninact))
        newa <- inactive[oo]
      }
      h <- min(h*frac.arclength, max.arclength/sum(abs(db)))
    }
    else {
      hd <- b[-force.active]/db[-force.active]
      ii <- hd > eps & hd < -h0
      if (any(ii)) h <- -max(hd[ii])
      else h = 0
    }
    list(h = -h, db = db, newa = newa, arclength = h*sum(abs(db)))
  }

"predictor1" <- function(b, step)
  {
    b - step$db * step$h
  }

"corrector1" <- function(x, y, family, weight, offset, active, tmpa, force.active, lambda, lambda2, b0, a0, bshoot.threshold, relax.lambda, trace, function.precision, no.iter = FALSE, eps = .Machine$double.eps)
  {
    variance <- family$variance
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    dev.resids <- family$dev.resids
    p <- length(tmpa)
    fk <- length(force.active)
    k <- p - fk
    if (!no.iter) {
      param <- c(b0[tmpa], a0, 0, -b0[tmpa[-force.active]]-a0, b0[tmpa[-force.active]]-a0)
      xa <- x[ ,tmpa,drop=FALSE]
      nobs <- nrow(xa)
      xstate <- rep(2, p+3*k+1)
      xstate[param==0] <- 0
      dstr <- switch(family$family, gaussian=0, binomial=1, poisson=2)
      lenz <- 10+(p+3)*nobs
      zsmall <- rep(0, lenz)
      zsmall[1:6] <- c(nobs, lambda, lambda2, dstr, k, function.precision)
      zsmall[10 + seq((p+3)*nobs)] <- c(as.vector(xa), y, weight, offset)
      sol <- .Fortran("solution",
                      k = as.integer(k),
                      n = as.integer(p+k),
                      nb = as.integer(p+3*k+1),
                      ne = as.integer(p+3*k),
                      hs = as.integer(xstate),
                      xn = as.double(param),
                      zsmall = as.double(zsmall),
                      lenz = as.integer(lenz),
                      inform = as.integer(0))
      b0[tmpa] <- sol$xn[1:p]
      if (k > 0) a0 <- sol$xn[(p+1):(p+k)]
      if (sol$inform != 0) cat("\nconvergence warning in corrector step\n")
    }
    Lambda2 <- c(0, rep(lambda2, length(b0)-1))
    eta <- drop(x%*%b0) + offset
    mu <- linkinv(eta)
    mu.eta.eta <- mu.eta(eta)
    w <- (weight * mu.eta.eta^2/variance(mu))^0.5
    z <- (y-mu)/mu.eta.eta
    xw <- x * w
    wz <- w * z
    corr <- drop(t(wz) %*% xw) - Lambda2*b0
    if (k > 0) {
      i <- which(abs(corr) >= lambda*(1-relax.lambda))
      newa <- i[!i%in%tmpa]
      newactive <- i[!i%in%active]
      i <- which(abs(b0[active[-force.active]]) < eps)
      inactive <- active[-force.active][i]
      active <- active[!active%in%inactive]
      active <- c(active, newactive)
      b0[-active] <- 0
    }
    else {
      lambda <- max(abs(corr[-active]))
      c1 <- which(abs(corr) == lambda)
      newa <- newactive <- c1
      inactive <- NULL
      active <- c(active, c1)
    }
    df <- length(active) - length(newactive)
    dev <- sum(dev.resids(y, mu, weight))
    backshoot <- ifelse(any(abs(b0[newactive]) > bshoot.threshold), TRUE, FALSE)
    list(b=b0, a=a0, active=active, force.active=force.active, newactive=newactive, newa=newa, inactive=inactive, corr=corr, lambda=lambda, xw=xw, df=df, dev=dev, backshoot=backshoot)
  }

"glmpath" <- function(x, y, data, nopenalty.subset = NULL, family = binomial, weight = rep(1, n), offset = rep(0, n), lambda2 = 1e-5, max.steps = 10*min(n, m), max.norm = 100*m, min.lambda = (if (m >= n) 1e-6 else 0), max.vars = Inf, max.arclength = Inf, frac.arclength = 1, add.newvars = 1, bshoot.threshold = 0.1, relax.lambda = 1e-8, standardize = TRUE, function.precision = 3e-13, eps = .Machine$double.eps, trace = FALSE)
  {
    call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    else if (family$family=="gaussian") family <- gaussian()
    else if (family$family=="binomial") family <- binomial()
    else if (family$family=="poisson") family <- poisson()    
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (!missing(data)) {
      x <- data$x
      y <- data$y
    }
    n <- nrow(x)
    m <- ncol(x)    
    no.iter <- FALSE
    if (family$family=="gaussian" & family$link=="identity") no.iter <- TRUE
    if (family$family=="binomial") {
      if (any(y==-1)) {
        y[y==-1] <- 0
        cat("y=-1 converted to y=0.\n")
      }
    }
    if (length(weight) != n) {
      stop("Length of the weight vector != sample size")
    }
    if (any(weight < 0)) {
      stop("Negative weights are not allowed.")
    }
    if (length(offset) != n) {
      stop("Length of the offset vector != sample size")
    }
    if (frac.arclength>1 | frac.arclength<=0) {
      frac.arclength <- 1
      cat("frac.arclength should be in (0,1]. frac.arclength is now set to 1.\n")
    }
    if (max.arclength<Inf & frac.arclength<1) {
      frac.arclength <- 1
      cat("frac.arclength<1 can be used only if max.arclength=Inf. frac.arclength is now set to 1.\n")
    }
    n.repeat <- n.repeat1 <- ceiling(1/frac.arclength)
    one <- rep(1, n)
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    if (standardize) {
      sdx <- sqrt(drop(one %*% (x^2))/(n-1))
      ignore <- sdx < eps
      if (any(ignore)) sdx[ignore] <- 1
    }
    else sdx <- rep(1, m)
    x <- scale(x, FALSE, sdx)
    x <- cbind(1, x)
    if (is.null(dimnames(x)[[2]])) {
      xnames <- c("Intercept", paste("x",seq(m),sep=""))
    }
    else xnames <- c("Intercept", dimnames(x)[[2]][-1])
    if (length(nopenalty.subset) >= m) {
      stop("cannot force in all the variables")
    }
    else nopenalty.subset <- c(1, nopenalty.subset+1)
    force.active <- c(1:length(nopenalty.subset))
    lam.vec <- step.len <- rep(0, max.steps)
    bmat.pred <- bmat.corr <- cmat <- matrix(0, max.steps, m+1)
    new.df <- df <- dev <- rep(0, max.steps)
    new.A <- rep(FALSE, max.steps)
    actions <- vector("list", length=max.steps)
    backshoot <- FALSE
    b <- rep(0, (m+1))
    corrector <- corrector1(x, y, family, weight, offset, nopenalty.subset, nopenalty.subset, force.active, 0, lambda2, b, NULL, bshoot.threshold, relax.lambda, function.precision, trace)
    if (trace) cat("The model begins with", xnames[nopenalty.subset])
    k <- 1
    b <- bmat.pred[k, ] <- bmat.corr[k, ] <- corrector$b
    cmat[k, ] <- corrector$corr
    lam.vec[k] <- lambda <- corrector$lambda
    new.df[k] <- df[k] <- corrector$df
    dev[k] <- corrector$dev
    new.A[k] <- TRUE
    active <- corrector$active
    actions[[k]] <- active[-1] - 1
    names(actions[[k]]) <- xnames[active[-1]]
    if (trace) {
      cat(paste("\nLambda=",lambda,"lets the first factor in.\nStep",k,":"))
      cat(paste("\t",xnames[active[-force.active]],"added"))
    }
    if (max.steps <= 1) {
      stop("Increase max.steps.")
    }
    if (max.norm <= sum(abs(b))) {
      stop("Increase max.norm.")
    }
    if (lambda <= min.lambda) {
      stop("Decrease min.lambda")
    }
    if (max.vars <= 1) {
      stop("Increase max.vars.")
    }
    while(TRUE) {
      if (!backshoot) {
        arclength <- max.arclength
        if (new.A[k]) {
          frac.arclength1 <- frac.arclength
        }
        else {
          frac.arclength1 <- 1
          if (n.repeat1 > 1 & max.arclength == Inf) {
            arclength <- step$arclength
          }
        }
        k <- k + 1
        if (trace) cat(paste("\nStep",k,":"))
        step <- step.length(corrector, lambda2, min.lambda, arclength, frac.arclength1, add.newvars, backshoot)
        b[active] <- predictor1(b[active], step)
        bmat.pred[k, ] <- b
        step.len[k-1] <- h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        tmpa <- c(active, step$newa)
        a <- abs(b[tmpa[-force.active]])
      }
      else {
        if (trace) cat(paste("\nStep",k,":"))
        step <- step.length(corrector, lambda2, min.lambda, Inf, 1, add.newvars, backshoot, h)
        step.len[k-1] <- h + step$h
        h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        a <- abs(b[tmpa[-force.active]])
      }
      corrector <- corrector1(x, y, family, weight, offset, active, tmpa, force.active, lambda, lambda2, b, a, bshoot.threshold, relax.lambda, trace, function.precision, no.iter)
      newa <- corrector$newa
      while(length(newa) > 0) {
        if (trace) cat(paste("\nRepeating step",k,":"))
        tmpa <- c(tmpa, newa)
        a <- abs(b[tmpa[-force.active]])
        corrector <- corrector1(x, y, family, weight, offset, active, tmpa, force.active, lambda, lambda2, b, a, bshoot.threshold, relax.lambda, trace, function.precision, no.iter)
        newa <- corrector$newa 
      }
      newaction <- c(corrector$newactive, -corrector$inactive)
      if (length(newaction) > 0 & length(corrector$active) <= n) {
        if (corrector$backshoot & !backshoot) {
          if (trace) cat("\nOvershooting occurred: increasing lambda again")
          backshoot <- TRUE
          n.repeat1 <- 1
        }
        else {
          active <- corrector$active
          b <- corrector$b
          actions[[k]] <- sign(newaction)*(abs(newaction) - 1)
          names(actions[[k]]) <- xnames[abs(newaction)]
          new.df[k] <- corrector$df
          new.A[k] <- TRUE
          if (trace) {
            na <- newaction[newaction > 0]
            ina <- -newaction[newaction < 0]
            if (length(na) > 0) cat(paste("\t",xnames[na],"added"))
            if (length(ina) > 0) cat(paste("\t",xnames[ina],"dropped"))
          }
          backshoot <- FALSE
          n.repeat1 <- n.repeat
        }
      }
      else if (length(corrector$active) <= n) {
        active <- corrector$active
        b <- corrector$b
        backshoot <- FALSE
        n.repeat1 <- max(n.repeat1-1, 1)
      }
      if (!backshoot) {
        bmat.corr[k, ] <- b 
        cmat[k, ] <- corrector$corr
        df[k] <- corrector$df
        dev[k] <- corrector$dev
        if (lambda <= min.lambda | k == max.steps | length(corrector$active) > min(n, max.vars) | sum(corrector$a) >= max.norm) {
          if (trace & lambda <= min.lambda) cat("\nLambda=",min.lambda,"\n")
          else if (trace & k == max.steps) cat(paste("\nMaximum steps (",max.steps,") taken.\n"))                     
          else if (length(corrector$active) > min(n, max.vars)) {
            k <- k - 1
            if (trace) cat("\nNumber of active variables has reached its maximum.\n")
          }
          else if (trace) cat("\n|beta| >=",max.norm,"\n")
          break
        }
      }
    }
    bmat.pred <- bmat.pred[1:k, ]
    bmat.corr <- bmat.corr[1:k, ]
    cmat <- cmat[1:k, ]
    bmat.pred[ ,-1] <- scale(bmat.pred[ ,-1], FALSE, sdx)
    bmat.corr[ ,-1] <- scale(bmat.corr[ ,-1], FALSE, sdx)
    bmat.pred[ ,1] <- bmat.pred[ ,1] - bmat.pred[ ,-1,drop=FALSE] %*% meanx
    bmat.corr[ ,1] <- bmat.corr[ ,1] - bmat.corr[ ,-1,drop=FALSE] %*% meanx
    dimnames(bmat.pred) <- dimnames(bmat.corr) <- dimnames(cmat) <- list(seq(k), xnames)
    df <- df[1:k]
    dev <- dev[1:k]
    if (family$family=="gaussian") Dev <- sum(weight)*(log(2*pi*dev/sum(weight)) + 1) + 2
    else Dev <- dev
    aic <- Dev + 2*df
    bic <- Dev + log(n)*df
    if (length(nopenalty.subset) > 1) {
      names(nopenalty.subset) <- xnames[nopenalty.subset]
      nopenalty.subset <- nopenalty.subset[-1]-1
    }
    else nopenalty.subset <- NULL
    object <- list(call = call, lambda=lam.vec[1:k], lambda2=lambda2, step.length=abs(step.len[1:(k-1)]), corr = cmat, new.df = new.df[1:k], df = df, deviance = dev, aic = aic, bic = bic, b.predictor = bmat.pred, b.corrector = bmat.corr, new.A = new.A[1:k], actions = actions[1:k], meanx = meanx, sdx = sdx, xnames = xnames, family = family, weight = weight, offset = offset, nopenalty.subset = nopenalty.subset, standardize = standardize)
    class(object) <- "glmpath"
    object
  }

"predict.glmpath" <- function(object, newx, newy, s, type = c("link", "response", "loglik", "coefficients"), mode = c("step", "norm.fraction", "norm", "lambda.fraction", "lambda"), weight = NULL, offset = NULL, eps = .Machine$double.eps, ...)
  {
    mode <- match.arg(mode)
    type <- match.arg(type)
    if (missing(newx) & type != "coefficients") {
        warning("no newx argument; type switched to coefficients")
        type <- "coefficients"
    }
    if (missing(newy) & type == "loglik") {
        warning("no newy argument; type switched to coefficients")
        type <- "coefficients"
    }
    if (type!="coefficients") {
        if (is.vector(newx)) newx <- matrix(newx, nrow=1)
    }
    if (any(object$offset != 0) & type != "coefficients") {
        if (is.null(offset)) {
          stop("offset missing")
        }
        if (length(offset) != nrow(newx)) {
          stop("Length of the offset vector != sample size in newx")
        }
    }
    b <- object$b.corrector
    std.b <- scale(b[ ,-1], FALSE, 1/object$sdx)
    if (!is.null(object$nopenalty.subset))
      std.b <- std.b[ ,-object$nopenalty.subset, drop=FALSE]
    k <- nrow(b)
    steps <- seq(k)
    if (missing(s)) {
        s <- steps[object$new.A]
        if (mode != "step") {
          warning("no s argument; mode switched to step")
          mode <- "step"
        }
    }
    sb <- switch(mode, step = {
        if (any(s < 1) | any(s > k)) 
            stop("Argument s out of range")
        steps
    }, norm.fraction = {
        if (any(s > 1) | any(s < 0)) 
            stop("Argument s out of range")
        bnorm <- apply(abs(std.b), 1, sum)
        bnorm / bnorm[k]
    }, norm = {
        bnorm <- apply(abs(std.b), 1, sum)
        if (any(s > bnorm[k]) | any(s < bnorm[1])) 
            stop("Argument s out of range")
        bnorm
    }, lambda.fraction = {
        if (any(s > 1) | any(s < 0))
            step("Argument s out of range")
        lam <- object$lambda
        lam[lam < eps] <- eps
        lam <- log(lam)
        (lam - min(lam)) / (max(lam) - min(lam))
    }, lambda = {
        lam <- object$lambda
        if (any(s > lam[1]) | any(s < lam[k]))
          stop("Argument s out of range")
        lam
    })
    sfrac <- (s - sb[1])/(sb[k] - sb[1])
    sb <- (sb - sb[1])/(sb[k] - sb[1])
    usb <- unique(sb)
    useq <- match(usb, sb)
    sb <- sb[useq]
    b <- b[useq, ]
    coord <- approx(sb, seq(sb), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newb <- ((sb[right] - sfrac) * b[left, , drop = FALSE] + 
        (sfrac - sb[left]) * b[right, , drop = FALSE])/(sb[right] - sb[left])
    newb[left == right, ] <- b[left[left == right], ]    
    if (type!="coefficients") {
      if (missing(offset)) offset <- rep(0, nrow(newx))
      fit <- cbind(1, newx) %*% t(newb) + offset
      if (type!="link") {
        fit <- object$family$linkinv(fit)
        if (type=="loglik") {
          if (is.null(weight)) weight <- rep(1, length(newy))
          loglik <- matrix(0, length(newy), length(s))
          for (k in 1:length(s)) loglik[ ,k] <- -object$family$dev.resids(newy, fit[ ,k], weight)/2
          fit <- loglik
        }
      }
      dimnames(fit) <- list(seq(nrow(newx)), s)
    }
    else {
      fit <- newb
      dimnames(fit) <- list(s, object$xnames)
    }
    attr(fit, "s") <- s
    attr(fit, "fraction") <- sfrac
    attr(fit, "mode") <- mode
    fit
}

