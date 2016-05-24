quickInteraction <- function(by){
  if(is.list(by)){
    n.arg <- length(by)
    f <- 0L
    uf <- 0L
    for(i in rev(1:n.arg)){
      y <- by[[i]]
      y <- as.numeric(y)
      uy <- unique(na.omit(y))
      y <- match(y,uy,NA)
      l <- length(uy)
      f <- f*l + y - 1
      uf <- unique(na.omit(f))
      f <- match(f,uf,NA)
      uf <- seq(length(uf))
    }
  }
  else {
    by <- as.numeric(by)
    uf <- unique(na.omit(by))
    f <- match(by,uf,NA)
    uf <- seq(length(uf))
  }
  return(structure(f,unique=uf))
}

constInSets <- function(X,sets){
    ans <- integer(0)
    for(i in 1:ncol(X)){
        v <- tapply(X[,i],sets,var)
        if(all(v[is.finite(v)]==0)) ans <- c(ans,i)
    }
    names(ans) <- colnames(X)[ans]
    ans
}

mclogit <- function(
                formula,
                data=parent.frame(),
                random=NULL,
                subset,
                weights=NULL,
                offset=NULL,
                na.action = getOption("na.action"),
                model = TRUE, x = FALSE, y = TRUE,
                contrasts=NULL,
                start.theta=NULL,
                control=mclogit.control(...),
                ...
            ){
# Assumptions:
#   left hand side of formula: cbind(counts,  choice set index)
#   right hand side of the formula: attributes
#   intercepts are removed!

    call <- match.call(expand.dots = TRUE)

    if(length(random))
        random <- setupRandomFormula(random)

    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "offset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    na.action <- attr(mf,"na.action")
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    Y <- as.matrix(model.response(mf, "any"))
    if(ncol(Y)<2) stop("need response counts and choice set indicators")
    sets <- Y[,2]
    sets <- match(sets,unique(sets))
    Y <- Y[,1]
    if (is.null(weights)){
        prior.weights <- rep(1,length(Y))
        N <- rowsum(Y,sets,na.rm=TRUE)
        weights <- N[sets]
        }
    else{
        prior.weights <- weights
        N <- rowsum(weights*Y,sets,na.rm=TRUE)
        weights <- N[sets]
        }
    N <- sum(N)
    Y <- Y/weights
    Y[weights==0] <- 0
    
    X <- model.matrix(mt,mf,contrasts)
    contrasts <- attr(X, "contrasts")
    xlevels <- .getXlevels(mt,mf)
    icpt <- match("(Intercept)",colnames(X),nomatch=0)
    if(icpt) X <- X[,-icpt,drop=FALSE]
    const <- constInSets(X,sets)
    if(length(const)){
        warning("removing ",paste(names(const),collapse=",")," from model")
        X <- X[,-const,drop=FALSE]
    }
    fit <- mclogit.fit(Y,sets,weights,X,
                        control=control,
                        offset = offset)
    null.dev <- fit$null.deviance
    if(length(random)){ ## random effects
        if(length(all.vars(random$covariates))){
          warning("random slopes not yet implemented, ignoring covariates")
          random$covarates <- ~1
        }
        mfr <- match.call(expand.dots = FALSE)
        mr <- match(c("formula", "data", "subset", "weights", "na.action"), names(mfr), 0)
        mfr <- mfr[c(1, mr)]
        environment(random$all.vars) <-environment(formula)
        mfr$formula <- random$all.vars
        mfr$drop.unused.levels <- TRUE
        mfr[[1]] <- as.name("model.frame")
        mfr <- eval(mfr, parent.frame())
        Z <- reDesignMatrix(random,mfr,use=if(length(na.action))-na.action else NULL)
        fit <- mclogit.fit.rePQL(Y,sets,weights,X,Z,
                              start=fit$coef,
                              start.theta=start.theta,
                              control=control,
                              offset = offset)
     }
    if(x) fit$x <- X
    if(x && length(random)) fit$z <- Z
    if(!y) {
        fit$y <- NULL
        ftt$s <- NULL
    }
    fit$null.deviance <- null.dev
    fit <- c(fit,list(call = call, formula = formula,
                      terms = mt,
                      random = random,
                      data = data,
                      contrasts = contrasts,
                      xlevels = xlevels,
                      na.action = na.action,
                      prior.weights=prior.weights,
                      weights=weights,
                      model=mf,
                      N=N))
    if(length(random))
        class(fit) <- c("mclogitRandeff","mclogit","lm")
    else
        class(fit) <- c("mclogit","lm")
    fit
}

mclogit.fit <- function(
      y,
      s,
      w,
      X,
      start=NULL,
      offset=NULL,
      control=mclogit.control()
      ){
    nvar <- ncol(X)
    deviance <- Inf
    eta <- mclogitLinkInv(y,s,w)
    nobs <- length(eta)
    if (is.null(offset))
        offset <- rep.int(0, nobs)

    pi <- mclogitP(eta,s)
    if(length(start))
      last.coef <- start
    else last.coef <- NULL
    converged <- FALSE
    for(iter in 1:control$maxit){
        y.star <- eta - offset + (y-pi)/pi
        yP.star <- y.star - rowsum(pi*y.star,s)[s]
        XP <- X - rowsum(pi*X,s)[s,,drop=FALSE]
        ww <- w*pi
        good <- ww > 0
        wlsFit <- lm.wfit(x=XP[good,,drop=FALSE],y=yP.star[good],w=ww[good])
        coef <- wlsFit$coefficients
        
        eta <- c(X%*%coef) + offset
        pi <- mclogitP(eta,s)
        last.deviance <- deviance
        dev.resids <- ifelse(y>0,
                2*w*y*(log(y)-log(pi)),
                0)
        deviance <- sum(dev.resids)
          ## check for divergence
          boundary <- FALSE
          if(!is.finite(deviance)){
            if(is.null(last.coef))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
             warning("step size truncated due to divergence", call. = FALSE)
             ii <- 1
             while (!is.finite(deviance)){
                if(ii > control$maxit)
                  stop("inner loop; cannot correct step size")
                ii <- ii + 1
                coef <- (coef + last.coef)/2
                eta <- c(X %*% coef) + offset
                pi <- mclogitP(eta,s)
                dev.resids <- ifelse(y>0,2*w*y*(log(y)-log(pi)),0)
                deviance <- sum(dev.resids)
             }
             boundary <- TRUE
             if (control$trace)
                  cat("Step halved: new deviance =", deviance, "\n")
          } ## inner loop
        if(control$trace)
            cat("\nIteration",iter,"- Deviance =",deviance)
        crit <- abs(deviance-last.deviance)/abs(0.1+deviance)
        if(crit < control$eps){
          converged <- TRUE
          if(control$trace)
            cat("\nconverged\n")
          break
        }
    }
    if (!converged) warning("algorithm did not converge")
    if (boundary) warning("algorithm stopped at boundary value")
    eps <- 10*.Machine$double.eps
    if (any(pi < eps) || any(1-pi < eps))
        warning("fitted probabilities numerically 0 occurred")

    XP <- X - rowsum(pi*X,s)[s,,drop=FALSE]
    ww <- w*pi
    Information <- crossprod(XP,ww*XP)
        
    ntot <- length(y)
    pi0 <- mclogitP(offset,s)
    null.deviance <- sum(ifelse(y>0,
                    2*w*y*(log(y)-log(pi0)),
                    0))
    resid.df <- length(y)#-length(unique(s))
    model.df <- ncol(X)
    resid.df <- resid.df - model.df
    ll <- mclogit.logLik(y,pi,w)
    return(list(
        coefficients = drop(coef),
        linear.predictors = eta,
        working.residuals = (y-pi)/pi,
        working.weights = w,
        response.residuals = y-pi,
        residual.df = resid.df,
        model.df = model.df,
        fitted.values = pi,
        deviance=deviance,
        ll=ll,
        deviance.residuals=dev.resids,
        null.deviance=null.deviance,
        iter = iter,
        y = y,
        s = s,
        offset = offset,
        converged = converged,
        control=control,
        covmat=solve(Information)
        ))
}

tr <- function(x) sum(diag(x))




mclogit.fit.rePQL <- function(
      y,
      s,
      w,
      X,
      Z,
      start=NULL,
      start.theta=NULL,
      offset=NULL,
      control=mclogit.control()
      ){
    #crossprod <- Matrix:::crossprod
    nvar <- ncol(X)
    deviance <- Inf
    eta <- mclogitLinkInv(y,s,w)
    nobs <- length(eta)
    if (is.null(offset))
      offset <- rep.int(0, nobs)
    lev.ics <- attr(Z,"col.indices")
    nlev <- length(lev.ics)
    sw <- c(tapply(w,s,"[",1))
    sqrt.sw <- sqrt(sw)
    nvalid <- length(s)
    coef <- start
    eta <- c(X%*%coef) + offset
    pi <- mclogitP(eta,s)
    
    dev.resids <- ifelse(y>0,2*w*y*(log(y)-log(pi)),0)
    deviance <- sum(dev.resids)
    converged <- FALSE
    lev.seq <- seq(length(lev.ics))

    ## Score Test for Variance Parameters
        
      ww <- w*pi
      P.sw <- sqrt(w)*pi*t(fac2sparse(s))
      Zw <- crossprod(P.sw,Z)
      ZWZ <- crossprod(Z,ww*Z) - crossprod(Zw)
      u <- crossprod(Z,w*(y-pi)) 
      l <- numeric(nlev)
      K <- matrix(nrow=nlev,ncol=nlev)
      for(i in lev.seq){
        ii <- lev.ics[[i]]
        ZWZ.i <- ZWZ[ii,ii]
        l[i] <- crossprod(u[ii]) - tr(ZWZ.i)
        K[i,i] <- sum(ZWZ.i*ZWZ.i)
        for(j in lev.seq[lev.seq > i]){
          jj <- lev.ics[[j]]
          ZWZ.ij <- ZWZ[ii,jj]
          K[i,j] <- K[j,i] <- sum(ZWZ.ij*ZWZ.ij)
        }
      }
    chisq.theta <- l^2/diag(K)
   ## Starting values for variance parameters
    if(length(start.theta)){
      if(length(start.theta)<length(lev.ics))
        stop("insufficient starting values for variance parameters")
      if(length(start.theta)>length(lev.ics))
        stop("to many starting values for variance parameters")
      last.theta <- theta <- start.theta
    }
    else{
      if(all(l < 0)) stop("insufficient residual variance; set some variance parameters to zero")
      theta <- solve(K,l)
      if(any(theta<0)){
        warning("Negative initial estimates found; correcting",call.=FALSE)
        cat("\ntheta=",theta)
        l[l < 0] <- 0
        cat("\nl=",l)
        ii <- 0
        bb <- 1
        I <- diag(x=diag(K))
        while(any(theta<0)){
          ii <- ii + 1
          cat("\nTrial ",ii)
          bb <- bb/2
          aa <- 1 - bb
#           if(ii > control$maxit)
#               stop("insufficient residual variance; set some variance parameters to zero")
          theta <- solve(aa*I+bb*K,l)
          cat("\ntheta=",theta)
        }
      }
      last.theta <- theta
    }
    if(control$trace)
      cat("\nInitial estimate of theta: ",theta)
    b <- rep(0,ncol(Z))
    rept <- sapply(lev.ics,length)
    iSigma <- as(Diagonal(x=rep(1/theta,rept)),"sparseMatrix")
    
    ## Extended IWLS and Fisher-scoring for variance parameters
    converged <- FALSE
    for(iter in 1:control$maxit){

      y.star <- eta - offset + (y-pi)/pi
      P.sw <- sqrt(w)*pi*t(fac2sparse(s))
      Py.star <- crossprod(P.sw,y.star)
      ww <- w*pi
      wX <- ww*X
      wZ <- ww*Z
      PX <- crossprod(P.sw,X)
      PZ <- crossprod(P.sw,Z)
      XWX <- crossprod(X,wX) - crossprod(PX)
      wZ <- crossprod(P.sw,Z)
      ZWZ <- crossprod(Z,ww*Z) - crossprod(PZ)    
      iZWZiSigma <- solve(ZWZ + iSigma)
      ZWX <- crossprod(Z,wX) - crossprod(PZ,PX)
      
      w.y.star <- ww*y.star
      XWy<- crossprod(X,w.y.star) - crossprod(PX,Py.star)
      ZWy<- crossprod(Z,w.y.star) - crossprod(PZ,Py.star)
      
      
      
      XiVX <- XWX - crossprod(ZWX,iZWZiSigma%*%ZWX)
      XiVy <- XWy - crossprod(ZWX,iZWZiSigma%*%ZWy)
      
      last.coef <- coef
      coef <- as.matrix(solve(XiVX,XiVy))
      
      last.b <- b
      b <- iZWZiSigma%*%(ZWy - ZWX%*%coef)
      #print(range(b))
      
      eta <- as.vector(X%*%coef) + as.vector(Z%*%b)
      pi <- mclogitP(eta,s)
      
      dev.resids <- ifelse(y>0,2*w*y*(log(y)-log(pi)),0)
      last.deviance <- deviance
      deviance <- sum(dev.resids)

      ## Fisher-scoring step for variance parameters
      grad.theta <- numeric(nlev)
      Info1.theta <- H2.theta <- matrix(0,ncol=nlev,nrow=nlev)
      for(i in lev.seq){
        ii <- lev.ics[[i]]
        m.i <- length(ii)
        A.i <- as.matrix(iZWZiSigma[ii,ii])
        grad.theta[i] <- (crossprod(b[ii]) -theta[i]*m.i + tr(A.i))/theta[i]^2
        Info1.theta[i,i] <- m.i/theta[i]^2 
        H2.theta[i,i] <- sum((A.i)^2)/theta[i]^4
        for(j in lev.seq[lev.seq > i]){
          jj <- lev.ics[[j]]
          A.ij <- as.matrix(iZWZiSigma[ii,jj])
          H2.theta[i,j] <- H2.theta[j,i] <- sum((A.ij)^2)/theta[i]^2/theta[j]^2
        }
      }
      Info.theta <- Info1.theta - H2.theta
      #print(grad.theta)
      #print(Info.theta)
      if(all(eigen(Info.theta)$values>0))
        diff.theta <- solve(Info.theta,grad.theta)
      else
        diff.theta <- solve(Info1.theta,grad.theta)
      theta <- theta + diff.theta
      #message("theta=",theta)
      if(any(theta<0)){
        ## Handle negative variances
        warning("negative values of variance parameters occured",call.=FALSE)
      }
 
      crit <- abs(deviance-last.deviance)/abs(0.1+deviance)
      ## Checking convergence
      #crit <- (control$eps)^(-2) *max(abs(delta.theta/(theta+.1)),abs(delta.coef1/(coef1+.1)))
      #crit.theta <- max(abs(delta.theta/(theta+.1)))
      if(control$trace)
          cat("\nIteration",iter,"- Deviance =",deviance,#"theta =",theta,
          "criterion = ",abs(deviance-last.deviance)/abs(0.1+deviance)#,
          #"criterion[2] = ",max(abs(theta - last.theta))
          )

      if(crit < control$eps){
        converged <- TRUE
        if(control$trace) cat("\nconverged\n")
        break
      }
   }
   if (!converged) warning("algorithm did not converge")
    
    P.sw <- sqrt(w)*pi*t(fac2sparse(s))
    ww <- w*pi
    wX <- ww*X
    wZ <- ww*Z
    PX <- crossprod(P.sw,X)
    PZ <- crossprod(P.sw,Z)
    XWX <- crossprod(X,wX) - crossprod(PX)
    wZ <- crossprod(P.sw,Z)
    ZWZ <- crossprod(Z,ww*Z) - crossprod(PZ)    
    iZWZiSigma <- solve(ZWZ + iSigma)
    ZWX <- crossprod(Z,wX) - crossprod(PZ,PX)

    XiVX <- XWX - crossprod(ZWX,iZWZiSigma%*%ZWX)
    covmat <- solve(XiVX)
    
    Info1.theta <- H2.theta <- matrix(0,ncol=nlev,nrow=nlev)
    for(i in lev.seq){
      ii <- lev.ics[[i]]
      m.i <- length(ii)
      A.i <- as.matrix(iZWZiSigma[ii,ii])
      Info1.theta[i,i] <- m.i/theta[i]^2 
      H2.theta[i,i] <- sum((A.i)^2)/theta[i]^4
      for(j in lev.seq[lev.seq > i]){
        jj <- lev.ics[[j]]
        A.ij <- as.matrix(iZWZiSigma[ii,jj])
        H2.theta[i,j] <- H2.theta[j,i] <- sum((A.ij)^2)/theta[i]^2/theta[j]^2
      }
    }
  
    if(all(eigen(Info.theta)$values>0))
      covmat.theta <- solve(Info.theta)
    else{
      warning("Fisher Information not positive definite, using simple approximation")
      covmat.theta <- solve(Info1.theta)
    }

   

   names(theta) <- names(lev.ics)
   colnames(covmat.theta) <- rownames(covmat.theta) <- names(theta)
   reff <- b
   reff <- lapply(lev.ics,function(ii)reff[ii])
   resid.df <- length(y)#-length(unique(s))
   model.df <- ncol(X) + length(theta)
   resid.df <- resid.df-model.df
   return(list(
      coefficients = drop(coef),
      varPar = theta,
      chisq.theta = chisq.theta,
      random.effects = reff,
      linear.predictors = eta,
      fitted.values = pi,
      ll=NA,
      deviance=deviance,
      deviance.residuals=dev.resids,
      working.residuals=(y-pi)/pi,
      response.residuals=y-pi,
      residual.df = resid.df,
      model.df = model.df,
      iter = iter,
      y = y,
      s = s,
      converged = converged,
      control=control,
      covmat=covmat,
      covmat.varPar = covmat.theta,
      rank=rank
      ))
}


setupRandomFormula <- function(formula){
  fo <- delete.response(terms(formula))
  attributes(fo) <- NULL
  if(length(fo[[2]]) < 2 || as.character(fo[[2]][1])!="|")
    stop("missing '|' operator")
  covariates <- fo
  groups <- fo
  covariates[2] <- fo[[2]][2]
  groups[2] <- fo[[2]][3]
  list(
    covariates=covariates,
    groups=groups,
    all.vars=reformulate(all.vars(fo))
    )
}

reDesignMatrix <- function(random,data,use=NULL){
  covariates <- all.vars(random$covariates)
  if(length(covariates)) warning("covariates not yet implemented")
  if(!length(use)) use <- TRUE
  groups <- data[use,all.vars(random$groups),drop=FALSE]
  gnames <- names(groups)
  n <- length(groups[[1]])
  nlev <- length(groups)
  groups[[1]] <- quickInteraction(groups[[1]])
  if(nlev>1)
    for(i in 2:nlev)
      groups[[i]] <- quickInteraction(groups[c(i-1,i)])
  un <- length(attr(groups[[1]],"unique"))
  Z <- Matrix(0,nrow=n,ncol=un,dimnames=list(NULL,paste(gnames[1],seq(un),sep="")))
  ij <- cbind(1:n,groups[[1]])
  Z[ij] <- 1
  lev.ics <- list()
  lev.ics[[1]] <- seq.int(un)
  if(nlev>1)
    for(i in 2:nlev){
      un.i <- length(attr(groups[[i]],"unique"))
      Z.i <- Matrix(0,nrow=n,ncol=un.i,dimnames=list(NULL,paste(gnames[i],seq(un.i),sep="")))
      ij <- cbind(1:n,groups[[i]])
      Z.i[ij] <- 1
      Z <- cbind2(Z,Z.i)
      lev.ics.i <- seq.int(un.i) + max(lev.ics[[i-1]])
      lev.ics <- c(lev.ics,list(lev.ics.i))
    }
  rand.names <- all.vars(random$groups)
  if(length(rand.names) > 1)
    for(i in 2:length(rand.names))
      rand.names[i] <- paste(rand.names[i-1],rand.names[i],sep=":")
  names(lev.ics) <- rand.names
  structure(Z,col.indices=lev.ics)
}

mclogit.control <- function(
                            epsilon = 1e-08,
                            maxit = 25,
                            trace=TRUE
                            ) {
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace)
}

log.Det <- function(x) determinant(x,logarithm=TRUE)$modulus

mclogitP <- function(eta,s){
  expeta <- exp(eta)
  sum.expeta <- rowsum(expeta,s)
  expeta/sum.expeta[s]
}

# mclogit.dev.resids <- function(y,p,w)
#       ifelse(y>0,
#                 2*w*y*(log(y)-log(p)),
#                 0)

mclogit.logLik <- function(y,p,w) sum(w*y*log(p))
                
                
mclogitLinkInv <- function(y,s,w){
  #n.alt <- tapply(y,s,length)
  #c(log(sqrt(w)*y+1/n.alt[s])-log(w)/2)
  n <- w*y+0.5
  f <- n/(rowsum(n,s)[s])
  log(f) - ave(log(f),s)
}

print.mclogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
    cat("\nCall: ", deparse(x$call), "\n\n")
    if(length(coef(x))) {
        cat("Coefficients")
        if(is.character(co <- x$contrasts))
            cat("  [contrasts: ",
                apply(cbind(names(co),co), 1, paste, collapse="="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n\n")
    if(length(x$varPar)) {
        cat("Variance paremeters")
        cat(":\n")
        print.default(format(x$varPar, digits=digits),
                      print.gap = 2, quote = FALSE)
    }
    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)))
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
    invisible(x)
}

vcov.mclogit <- function(object,...){
    return(object$covmat)
}

weights.mclogit <- function(object,...){
  return(object$weights)
}

deviance.mclogit <- function(object,...){
  return(object$deviance)
}

summary.mclogit <- function(object,dispersion=NULL,correlation = FALSE, symbolic.cor = FALSE,...){

    ## calculate coef table

    coef <- object$coefficients
    covmat.scaled <- object$covmat 
    var.cf <- diag(covmat.scaled)
    s.err <- sqrt(var.cf)
    zvalue <- coef/s.err
    pvalue <- 2*pnorm(-abs(zvalue))

    coef.table <- array(NA,dim=c(length(coef),4))
    dimnames(coef.table) <- list(names(coef),
            c("Estimate", "Std. Error","z value","Pr(>|z|)"))
    coef.table[,1] <- coef
    coef.table[,2] <- s.err
    coef.table[,3] <- zvalue
    coef.table[,4] <- pvalue

    if(length(object$varPar)){
      covmat.vp <- object$covmat.varPar 
      varPar <- object$varPar
      var.vp <- diag(covmat.vp)
      s.err.vp <- sqrt(var.vp)
      #zvalue.vp <- varPar/s.err.vp
      #pvalue.vp <- 2*pnorm(-abs(zvalue.vp))
      zvalue.vp <- sqrt(object$chisq.theta)
      pvalue.vp <- pchisq(object$chisq.theta,1,lower.tail=FALSE)


      varPar.table <- array(NA,dim=c(length(varPar),4))
      dimnames(varPar.table) <- list(names(varPar),
            c("Estimate", "Std. Error","z value","Pr(>|z|)"))
      varPar.table[,1] <- varPar
      varPar.table[,2] <- s.err.vp
      varPar.table[,3] <- zvalue.vp
      varPar.table[,4] <- pvalue.vp
    } else varPar.table <- NULL

    ans <- c(object[c("call","terms","baseline","deviance","contrasts",
                       "null.deviance","iter","na.action","model.df","residual.df")],
              list(coefficients = coef.table,
                    varPar = varPar.table,
                    cov.coef=object$covmat,
                    cov.varPar=object$covmat.varPar))
    p <- length(coef)
    if(correlation && p > 0) {
        dd <- sqrt(diag(ans$cov.coef))
        ans$correlation <-
            ans$cov.coef/outer(dd,dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.mclogit"
    return(ans)
}

print.summary.mclogit <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
              signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

    coefs <- x$coefficients
    if(length(x$varPar)){
      varPar <- x$varPar
      rownames(varPar) <- paste("Var(",rownames(varPar),")",sep="")
      coefs <- rbind(coefs,varPar)
    }
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)
#     cat("\n")
#     cat("AIC: ", format(x$aic, digits= max(4, digits+1)),"\n\n",
#         "Number of Fisher Scoring iterations: ", x$iter,
#         "\n", sep="")
    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\nNumber of Fisher Scoring iterations: ", x$iter,
        "\nNumber of observations: ",x$N,
        "\n")
    correl <- x$correlation
    if(!is.null(correl)) {
        p <- NCOL(correl)
        if(p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if(is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
        }
    }

    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
    invisible(x)
}




fitted.mclogit <- function(object,type=c("probabilities","counts"),...){
  weights <- object$weights
  nobs <- length(weights)
  res <- object$fitted.values
  type <- match.arg(type)
  
  na.act <- object$na.action
  
  res <- switch(type,
          probabilities=res,
          counts=weights*res)
          
  if(is.null(na.act))
    res
  else 
    napredict(na.act,res)
}

predict.mclogit <- function(object, newdata=NULL,type=c("link","response"),se=FALSE,...){

  type <- match.arg(type)
  rhs <- object$formula[-2]
  if(missing(newdata))
    m <- model.frame(rhs,data=object$model)
  else
    m <- model.frame(rhs,data=newdata)
  X <- model.matrix(rhs,m,
          contasts.arg=object$contrasts,
          xlev=object$xlevels
          )
  drop <- match("(Intercept)",colnames(X))
  X <- X[,-drop,drop=FALSE]
  eta <- c(X %*% object$coef)
  if(se){
    stopifnot(dim(X) == dim(X %*% vcov(object)))
    se.eta <- sqrt(rowSums(X * (X %*% vcov(object))))
  }

  na.act <- object$na.action
  
  if(type=="response") {
    lhs <- object$formula[[2]]
    set <- lhs[[3]]
    set <- eval(set,newdata,parent.frame())
    set <- match(set,unique(set))
    exp.eta <- exp(eta)
    sum.exp.eta <- rowsum(exp.eta,set)
    p <- exp.eta/sum.exp.eta[set]
    if(se){
      if(is.null(na.act))
        list(pred=p,se.pred=p*(1-p)*se.eta)
      else
        list(pred=napredict(na.act,p),
             se.pred=napredict(na.act,p*(1-p)*se.eta))
    }
    else {
      if(is.null(na.act)) p
      else napredict(na.act,p)
    }
  }
  else if(se) {
    if(is.null(na.act))
      list(pred=eta,se.pred=se.eta) 
    else
      list(pred=napredict(na.act,eta),
           se.pred=napredict(na.act,se.eta))
    }
  else {
    if(is.null(na.act)) eta
    else napredict(na.act,eta)
  }
}

logLik.mclogit <- function(object,...){
    if (length(list(...)))
        warning("extra arguments discarded")
    val <- if(length(object$ll))
            object$ll
           else NA
    attr(val, "nobs") <- object$N
    attr(val, "df") <- object$model.df
    class(val) <- "logLik"
    return(val)
}

residuals.mclogit <- 
 function(object, type = c("deviance", "pearson", "working",
                           "response", "partial"), ...){
   type <- match.arg(type)
   
   resid <- switch(type,
                 deviance=mclogit.dev.resids(object),
                 pearson=stop("not yet implemented"),
                 working=object$working.residuals,
                 response=object$response.residuals,
                 partial=stop("not yet implemented")
                 )
   naresid(object$na.action,resid)
 }
 
mclogit.dev.resids <- function(obj){
 y <- obj$y
 s <- obj$s
 w <- obj$weights
 pi <- obj$fitted.values
 
 n <- w*y+0.5
 f <- n/(rowsum(n,s)[s])
 #sign(y-p)*sqrt(2*abs(log(f)-log(y)))
 r <- 2*(f*log(f/pi))
 r - ave(r,s)
}