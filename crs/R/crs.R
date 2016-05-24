## This function provides for a complete S3 implementation of
## regression splines with categorical factors using two approaches,
## (i) kernel smoothing, and (ii) indicator function bases. Both
## additive and tensor product bases are supported (default is
## additive, but also see the option basis="auto" that computes both
## and uses that with the smallest cross-validation
## score). Cross-validation (leave-one-out, generalized, and the AIC_c
## method of Hurvich, Simonoff, and Tsai (1998, JRSS B)) can be used
## to select (i) the degree and number of knots (`segments'+1) of the
## basis spline for each continuous predictor, (ii) bandwidth for each
## ordinal/nominal predictor, or (iii) whether or not to include each
## ordinal/nominal predictor's indicator basis. S3 methods include
## fitted, predict, residuals, plot and so forth.

## 2010 (C) Jeffrey S. Racine (racinej@mcmaster.ca).

crs <- function(...) UseMethod("crs")

## This function computes the fit and returns the fit, degree
## (vector), segments (vector), and include (vector) for categorical
## predictors. Note that degree of zero and include of zero drop the
## variable from the resulting fit.

crsEst <- function(xz,
                   y,
                   degree=NULL,
                   segments=NULL,
                   include=NULL,
                   kernel=TRUE,
                   lambda=NULL,
                   complexity=c("degree-knots","degree","knots"),
                   knots=c("quantiles","uniform","auto"),
                   basis=c("additive","tensor","glp","auto"),
                   deriv=0,
                   data.return=FALSE,
                   prune=FALSE,
                   prune.index=NULL,
                   model.return=FALSE,
                   tau=NULL,
                   weights=NULL) {

  ## Take data frame xz and parse into factors (z) and numeric (x).

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)

  if(!kernel) {
    xztmp <- splitFrame(xz)
  } else {
    xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
  }
  x <- xztmp$x
  xnames <- xztmp$xnames
  num.x <- xztmp$num.x
  z <- xztmp$z
  znames <- xztmp$znames
  num.z <- xztmp$num.z
  is.ordered.z <- xztmp$is.ordered.z
  ## The default is kernel==TRUE - this will throw an error with no
  ## categorical predictors so first check
  if(is.null(num.z) && kernel==TRUE) kernel <- FALSE
  rm(xztmp)
  if(is.null(z)) {
    include <- NULL
  }

  y <- as.numeric(y)

  ## If weights are provided and there are NA values in the data, we
  ## need only those weights corresponding to the complete cases

  if(!is.null(weights))
    weights <- na.omit(data.frame(xz,y,weights))$weights

  if(!kernel) {

    model <- predict.factor.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=cbind(degree,segments),
                                   I=include,
                                   knots=knots,
                                   basis=basis,
                                   prune=prune,
                                   tau=tau,
                                   weights=weights)

    prune.index <- model$prune.index

    if(deriv > 0) {

      deriv.mat <- xz ## copy for dimension only
      deriv.mat.lwr <- deriv.mat
      deriv.mat.upr <- deriv.mat
      l <- 1 ## num.z
      m <- 1 ## num.x
      for(i in 1:ncol(xz)) {
        if(!is.factor(xz[,i])) {
          if(deriv <= degree[m]) {
            tmp <- deriv.factor.spline(x=x,
                                       y=y,
                                       z=z,
                                       K=cbind(degree,segments),
                                       I=include,
                                       knots=knots,
                                       basis=basis,
                                       deriv.index=m,
                                       deriv=deriv,
                                       prune.index=prune.index,
                                       tau=tau,
                                       weights=weights)
          } else {
            tmp <- matrix(0,length(y),3)
          }
          deriv.mat[,i] <- tmp[,1]
          deriv.mat.lwr[,i] <- tmp[,2]
          deriv.mat.upr[,i] <- tmp[,3]
          rm(tmp)
          m <- m + 1
        } else {
          ztmp <- z
          ztmp[,l] <- factor(rep(levels(xz[,i])[1],NROW(xz)),levels=levels(xz[,i]),ordered=is.ordered(xz[,i]))

          zpred <- predict.factor.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=cbind(degree,segments),
                                         I=include,
                                         knots=knots,
                                         basis=basis,
                                         prune=prune,
                                         prune.index=prune.index,
                                         tau=tau,
                                         weights=weights)$fitted.values

          zpred.base <- predict.factor.spline(x=x,
                                              y=y,
                                              z=z,
                                              K=cbind(degree,segments),
                                              I=include,
                                              xeval=x,
                                              zeval=ztmp,
                                              knots=knots,
                                              basis=basis,
                                              prune=prune,
                                              prune.index=prune.index,
                                              tau=tau,
                                              weights=weights)$fitted.values

          deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
          deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
          deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

          l <- l + 1
        }
      }

    } else {
      deriv.mat <- NULL
      deriv.mat.lwr <- NULL
      deriv.mat.upr <- NULL
    }

  } else {

    model <- predict.kernel.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=cbind(degree,segments),
                                   lambda=lambda,
                                   is.ordered.z=is.ordered.z,
                                   knots=knots,
                                   basis=basis,
                                   model.return=model.return,
                                   tau=tau,
                                   weights=weights)

    prune.index <- NULL

    if(deriv > 0) {

      deriv.mat <- xz ## copy for dimension only
      deriv.mat.lwr <- deriv.mat
      deriv.mat.upr <- deriv.mat
      l <- 1 ## num.z
      m <- 1 ## num.x
      for(i in 1:ncol(xz)) {
        if(!is.factor(xz[,i])) {
          if(deriv <= degree[m]) {
            tmp <- deriv.kernel.spline(x=x,
                                       y=y,
                                       z=z,
                                       K=cbind(degree,segments),
                                       lambda=lambda,
                                       is.ordered.z=is.ordered.z,
                                       knots=knots,
                                       basis=basis,
                                       deriv.index=m,
                                       deriv=deriv,
                                       tau=tau,
                                       weights=weights)
          } else {
            tmp <- matrix(0,length(y),3)
          }
          deriv.mat[,i] <- tmp[,1]
          deriv.mat.lwr[,i] <- tmp[,2]
          deriv.mat.upr[,i] <- tmp[,3]
          rm(tmp)
          m <- m + 1
        } else {

          ztmp <- z
          ztmp[,l] <- rep(sort(unique(z[,l]))[1],NROW(z))
          
          zpred <- predict.kernel.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=cbind(degree,segments),
                                         lambda=lambda,
                                         is.ordered.z=is.ordered.z,
                                         knots=knots,
                                         basis=basis,
                                         model.return=model.return,
                                         tau=tau,
                                         weights=weights)$fitted.values

          zpred.base <- predict.kernel.spline(x=x,
                                              y=y,
                                              z=z,
                                              K=cbind(degree,segments),
                                              lambda=lambda,
                                              is.ordered.z=is.ordered.z,
                                              xeval=x,
                                              zeval=ztmp,
                                              knots=knots,
                                              basis=basis,
                                              model.return=model.return,
                                              tau=tau,
                                              weights=weights)$fitted.values

          deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
          deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
          deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

          l <- l + 1
        }
      }

    } else {
      deriv.mat <- NULL
      deriv.mat.lwr <- NULL
      deriv.mat.upr <- NULL
    }

  }

  if(!data.return) {
    x <- NULL
    z <- NULL
  }

  if(is.null(tau)) {
    fitted.values <- model$fitted.values[,1]
    lwr <- model$fitted.values[,2]
    upr <- model$fitted.values[,3]
  } else {
    fitted.values <- model$fitted.values[,1]
    lwr <- NULL
    upr <- NULL
  }

  return(list(fitted.values=fitted.values,
              lwr=lwr,
              upr=upr,
              df.residual=model$df.residual,
              K=cbind(degree,segments),
              degree=degree,
              segments=segments,
              complexity=complexity,
              knots=knots,
              include=include,
              lambda=lambda,
              kernel=kernel,
              basis=basis,
              num.x=num.x,
              num.z=num.z,
              is.ordered.z=is.ordered.z,
              xnames=xnames,
              znames=znames,
              deriv=deriv,
              deriv.mat=deriv.mat,
              deriv.mat.lwr=deriv.mat.lwr,
              deriv.mat.upr=deriv.mat.upr,
              model.lm=model$model,
              hatvalues=model$hatvalues,
              nobs=length(y),
              k=model$rank,
              x=x,
              z=z,
              prune=prune,
              prune.index=prune.index,
              P.hat=model$P.hat,
              tau=tau,
              weights=weights))

}

## Default method - this function takes the minimum arguments (data,
## degree of spline with one element for each column of xz having
## continuous data (presumed default is all xz continuous)).

crs.default <- function(xz,
                        y,
                        degree=NULL,
                        segments=NULL,
                        include=NULL,
                        kernel=TRUE,
                        lambda=NULL,
                        complexity=c("degree-knots","degree","knots"),
                        knots=c("quantiles","uniform","auto"),
                        basis=c("additive","tensor","glp","auto"),
                        deriv=0,
                        data.return=FALSE,
                        prune=FALSE,
                        model.return=FALSE,
                        tau=NULL,
                        weights=NULL,
                        ...) {

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)

  ## Does the following properly belong here or crsEst?

  est <- crsEst(xz=xz,
                y=y,
                degree=degree,
                segments=segments,
                include=include,
                kernel=kernel,
                lambda=lambda,
                complexity=complexity,
                knots=knots,
                basis=basis,
                deriv=deriv,
                data.return=data.return,
                prune=prune,
                model.return=model.return,
                tau=tau,
                weights=weights)

  ## Add results to estimated object.

  est$residuals <- y - est$fitted.values
  est$r.squared <- RSQfunc(y,est$fitted.values,weights)
  est$call <- match.call()
  class(est) <- "crs"

  ## Return object of type crs

  return(est)

}

## Here we define the formula and split y (always first column of the
## model frame) from xz (the remaining continuous and
## ordinal/nominal).  nomad::FALSE exhaustive search nmulti is the
## number for multiple initial points.  if it is bigger than 1, when
## nomad is true, it will call snomadRSolve, otherwise, it will call
## smultinomadRSolve See ?snomadr
## Jun 4,  2011
##1) degree.max (we have removed  basis.maxdim)
##2) segments.max (we have removed  basis.maxdim)
##3) degree.min (currently 0)
##4) segments.min (currently 1)

crs.formula <- function(formula,
                        data=list(),
                        degree=NULL,
                        segments=NULL,
                        include=NULL,
                        degree.max=10,
                        segments.max=10,
                        degree.min=0,
                        segments.min=1,
                        cv.df.min=1,
                        cv=c("nomad","exhaustive","none"),
                        cv.threshold=1000,
                        cv.func=c("cv.ls","cv.gcv","cv.aic"),
                        kernel=TRUE,
                        lambda=NULL,
                        lambda.discrete=FALSE,
                        lambda.discrete.num=100,
                        complexity=c("degree-knots","degree","knots"),
                        knots=c("quantiles","uniform","auto"),
                        basis=c("auto","additive","tensor","glp"),
                        deriv=0,
                        data.return=FALSE,
                        prune=FALSE,
                        model.return=FALSE,
                        restarts=0,
                        random.seed=42,
                        max.bb.eval=10000,
                        initial.mesh.size.real="r1.0e-01", 
                        initial.mesh.size.integer="1",
                        min.mesh.size.real=paste(sqrt(.Machine$double.eps)),
                        min.mesh.size.integer=paste(sqrt(.Machine$double.eps)),
                        min.poll.size.real=paste(sqrt(.Machine$double.eps)),
                        min.poll.size.integer=paste(sqrt(.Machine$double.eps)),
                        opts=list(),
                        nmulti=5,
                        tau=NULL,
                        weights=NULL,
                        singular.ok=FALSE,
                        ...) {

  cv <- match.arg(cv)
  cv.func <- match.arg(cv.func)
  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)

#  if(!is.null(tau)) {
#    if(!require(quantreg)) stop(" Error: you must first install the quantreg package")
#  }

  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  xz <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

  ## Set DISPLAY_DEGREE to 0 if crs.messages=FALSE and DISPLAY_DEGREE
  ## is not provided

  if(!options('crs.messages')$crs.messages && is.null(opts[["DISPLAY_DEGREE"]])) opts$"DISPLAY_DEGREE"=0

  ## If a weights vector is provided and there exists missing data
  ## then the weight vector must be parsed to contain weights
  ## corresponding to the non-missing observations only.

  rows.omit <- as.vector(attr(mf,"na.action"))
  if(!is.null(weights) && !is.null(rows.omit))
    weights <- weights[-rows.omit]

  if(!kernel) {
    xztmp <- splitFrame(xz)
  } else {
    xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
  }
  x <- xztmp$x
  xnames <- xztmp$xnames
  num.x <- xztmp$num.x
  z <- xztmp$z
  znames <- xztmp$znames
  num.z <- xztmp$num.z
  is.ordered.z <- xztmp$is.ordered.z
  ## The default is kernel==TRUE - this will throw an error with no
  ## categorical predictors so first check
  if(is.null(num.z) && kernel==TRUE) kernel <- FALSE
  rm(xztmp)
  if(is.null(z)) {
    include <- NULL
  }

  ## Check for dynamic cv and if number of combinations is not overly
  ## large use exhaustive search

  if(cv=="nomad" && is.null(num.z) && ((degree.max-degree.min)*(segments.max-segments.min))**num.x <= cv.threshold) {
    warning(" Dynamically changing search from nomad to exhaustive (if unwanted set cv.threshold to 0)")
    cv <- "exhaustive"
  }

  ## If no degree nor include nor lambda, return cubic spline
  ## (identity bases) or non-smooth model (kernel).

  if(cv.df.min < 1 || cv.df.min > length(y)-1) stop(" cv.df.min must be a positive integer less than n")
  if(!is.null(degree)&&length(degree)!=num.x) stop(" degree vector must be the same length as x")
  if(!is.null(segments)&&length(segments)!=num.x) stop(" segments vector must be the same length as x")
  if(degree.max > 100) stop(paste(" degree.max (",degree.max,") exceeds reasonable value (",100,")",sep=""))
  if(lambda.discrete && (lambda.discrete.num < 1)) stop(" lambda.discrete.num must be a positive integer")

  if(cv=="none"){

    ## When no cross-validation is selected and no defaults are set
    ## for various parameters, we set them to ad hoc defaults and warn
    ## the user to this effect.

    if(is.null(degree)&!is.null(x)) {
        warning(paste(" cv=\"none\" selected but no degree provided, using degree=rep(3,num.x): you might consider other degree settings",sep=""),immediate.=TRUE)
        degree <- rep(3,num.x)
      }
      if(is.null(segments)&!is.null(x)) {
        warning(paste(" cv=\"none\" selected but no segments provided, using segments=rep(1,num.x): you might consider other segment settings",sep=""),immediate.=TRUE)
        segments <- rep(1,num.x)
      }
      if(is.null(include)&!is.null(z)&!kernel) {
        warning(paste(" cv=\"none\" selected but no inclusion for factors indicated, using include=rep(1,num.z): you might consider other include settings",sep=""),immediate.=TRUE)
        include <- rep(1,num.z)
      }
      if(is.null(lambda)&!is.null(z)&&kernel) {
        warning(paste(" cv=\"none\" selected but no bandwidths for factors indicated, using lambda=rep(0,num.z): you might consider other lambda settings",sep=""),immediate.=TRUE)
        lambda <- rep(0,num.z)
      }

      ## With one continuous predictor all bases are identical, so
      ## simply set the basis to additive and be done (no warning
      ## necessary)

      if(basis=="auto"&&num.x==1) basis <- "additive"

      if(basis=="auto"&&num.x>1) {
        warning(paste(" cv=\"none\" selected, basis=\"auto\" changed to basis=\"additive\": you might consider basis=\"tensor\" etc.",sep=""),immediate.=TRUE)
        basis <- "additive"
      }

      if(knots=="auto"&&num.x>1) {
        warning(paste(" cv=\"none\" selected, knots=\"auto\" changed to knots=\"quantiles\": you might consider knots=\"uniform\" etc.",sep=""),immediate.=TRUE)
        knots <- "quantiles"
      }

  }

  if(kernel==TRUE&&prune==TRUE) warning(" pruning cannot coexist with categorical kernel smoothing (pruning ignored)")
  if(!is.null(tau)&&prune==TRUE) stop(" pruning is not supported for quantile regression splines")

  ## Check for cv="nomad" and complexity="degree-knots" but
  ## degree.min==degree.max or segments==segments.max

  if((cv=="nomad" && complexity=="degree-knots") && (segments.min==segments.max)) stop("NOMAD search selected with complexity degree-knots but segments.min and segments.max are equal")
  if((cv=="nomad" && complexity=="degree-knots") && (degree.min==degree.max)) stop("NOMAD search selected with complexity degree-knots but degree.min and degree.max are equal")

  ## Check for proper derivative

  if(deriv < 0) stop("derivative order must be a non-negative integer")

  ## Check for logical singular.ok

  if(!is.logical(singular.ok)) stop("singular.ok must be logical (TRUE/FALSE)")

  cv.min <- NULL
  ptm <- system.time("")

  if(!kernel) {

    ## indicator bases and B-spline bases cross-validation

    if(cv=="nomad") {

      ptm <- ptm + system.time(cv.return <- frscvNOMAD(xz=xz,
                                                       y=y,
                                                       degree.max=degree.max,
                                                       segments.max=segments.max,
                                                       degree.min=degree.min,
                                                       segments.min=segments.min,
                                                       cv.df.min=cv.df.min,
                                                       complexity=complexity,
                                                       knots=knots,
                                                       basis=basis,
                                                       cv.func=cv.func,
                                                       degree=degree,
                                                       segments=segments,
                                                       include=include,
                                                       random.seed=random.seed,
                                                       max.bb.eval=max.bb.eval,
                                                       initial.mesh.size.integer=initial.mesh.size.integer,
                                                       min.mesh.size.integer=min.mesh.size.integer,
                                                       min.poll.size.integer=min.poll.size.integer,
                                                       opts=opts,
                                                       nmulti=nmulti,
                                                       tau=tau,
                                                       weights=weights,
                                                       singular.ok=singular.ok))

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min,sqrt(.Machine$double.xmax)))) stop(" Search failed: restart with larger nmulti or smaller degree.max  (or degree if provided)")
    }  else if(cv=="exhaustive") {

      ptm <- ptm + system.time(cv.return <- frscv(xz=xz,
                                                  y=y,
                                                  degree.max=degree.max,
                                                  segments.max=segments.max,
                                                  degree.min=degree.min,
                                                  segments.min=segments.min,
                                                  complexity=complexity,
                                                  knots=knots,
                                                  basis=basis,
                                                  cv.func=cv.func,
                                                  degree=degree,
                                                  segments=segments,
                                                  tau=tau,
                                                  weights=weights,
                                                  singular.ok=singular.ok))

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min,sqrt(.Machine$double.xmax)))) stop(" Search failed: restart with smaller degree.max")
    }

  } else {

    ## kernel smooth and B-spline bases cross-validation

    if(cv=="nomad") {

      ptm <- ptm + system.time(cv.return <- krscvNOMAD(xz=xz,
                                                       y=y,
                                                       degree.max=degree.max,
                                                       segments.max=segments.max,
                                                       degree.min=degree.min,
                                                       segments.min=segments.min,
                                                       cv.df.min=cv.df.min,
                                                       complexity=complexity,
                                                       knots=knots,
                                                       basis=basis,
                                                       cv.func=cv.func,
                                                       degree=degree,
                                                       segments=segments,
                                                       lambda=lambda,
                                                       lambda.discrete=lambda.discrete,
                                                       lambda.discrete.num=lambda.discrete.num,
                                                       random.seed=random.seed,
                                                       max.bb.eval=max.bb.eval,
                                                       initial.mesh.size.real=initial.mesh.size.real,
                                                       initial.mesh.size.integer=initial.mesh.size.integer,
                                                       min.mesh.size.real=min.mesh.size.real,
                                                       min.mesh.size.integer=min.mesh.size.integer,
                                                       min.poll.size.real=min.poll.size.real,
                                                       min.poll.size.integer=min.poll.size.integer,
                                                       opts=opts,
                                                       nmulti=nmulti,
                                                       tau=tau,
                                                       weights=weights,
                                                       singular.ok=singular.ok))

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      lambda <- cv.return$lambda
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min,sqrt(.Machine$double.xmax)))) stop(" Search failed: restart with larger nmulti or smaller degree.max")

    } else if(cv=="exhaustive") {

      ptm <- ptm + system.time(cv.return <- krscv(xz=xz,
                                                  y=y,
                                                  degree.max=degree.max,
                                                  segments.max=segments.max,
                                                  degree.min=degree.min,
                                                  segments.min=segments.min,
                                                  complexity=complexity,
                                                  knots=knots,
                                                  basis=basis,
                                                  cv.func=cv.func,
                                                  degree=degree,
                                                  segments=segments,
                                                  restarts=restarts,
                                                  tau=tau,
                                                  weights=weights,
                                                  singular.ok=singular.ok))

      cv.min <- cv.return$cv.objc
      degree <- cv.return$degree
      segments <- cv.return$segments
      include <- cv.return$I
      lambda <- cv.return$lambda
      basis <- cv.return$basis
      knots <- cv.return$knots
      if(isTRUE(all.equal(cv.min,sqrt(.Machine$double.xmax)))) stop(" Search failed: restart with smaller degree.max")

    }

  }

  ptm <- ptm + system.time(est <- crs.default(xz=xz,
                                              y=y,
                                              degree=degree,
                                              segments=segments,
                                              include=include,
                                              kernel=kernel,
                                              lambda=lambda,
                                              complexity=complexity,
                                              knots=knots,
                                              basis=basis,
                                              deriv=deriv,
                                              data.return=data.return,
                                              prune=prune,
                                              model.return=model.return,
                                              tau=tau,
                                              weights=weights,
                                              ...))


  est$cv.score <- cv.min
  est$call <- match.call()
  est$formula <- formula
  est$terms <- mt
  est$xlevels <- .getXlevels(mt, mf)
  est$xz <- xz
  est$y <- y
  est$prune <- prune
  est$cv.min <- cv.min
  est$cv <- cv
  est$restarts <- restarts
  est$ptm <- ptm
  est$nmulti <- nmulti

  return(est)

}

## Method for predicting given a new data frame.

predict.crs <- function(object,
                        newdata=NULL,
                        deriv=0,
                        ...) {

  if(is.null(newdata)) {

    ## If no new data provided, return sample fit.
    fitted.values <- fitted(object)
    deriv.mat <- object$deriv.mat

    lwr <- NULL
    upr <- NULL
    deriv.mat.lwr <- NULL
    deriv.mat.upr <- NULL

  } else{

    ## Get training data from object (xz and y) and parse into factors
    ## and numeric.

    basis <- object$basis
    deriv <- object$deriv
    prune <- object$prune
    prune.index <- object$prune.index
    tau <- object$tau
    weights <- object$weights

    xz <- object$xz
    y <- object$y

    ## Divide into factors and numeric

    if(!object$kernel) {
      xztmp <- splitFrame(xz)
    } else {
      xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
    }
    x <- xztmp$x
    z <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z
    rm(xztmp)

    ## Get evaluation data (newdata) and divide into factors and
    ## numeric.

    Terms <- delete.response(terms(object))
    newdata <- model.frame(Terms,newdata,xlev=object$xlevels)

    if(!object$kernel) {
      xztmp <- splitFrame(data.frame(newdata))
    } else {
      xztmp <- splitFrame(data.frame(newdata),factor.to.numeric=TRUE)
    }
    xeval <- xztmp$x
    zeval <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z
    rm(xztmp)

    ## Compute the predicted values.

    if(!object$kernel) {

      ## Get degree vector and include vector.

      complexity <- object$complexity
      knots <- object$knots
      K <- object$K
      degree <- object$degree
      segments <- object$segments
      include <- object$include

      tmp <- predict.factor.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=K,
                                   I=include,
                                   xeval=xeval,
                                   zeval=zeval,
                                   basis=basis,
                                   knots=knots,
                                   prune=prune,
                                   prune.index=prune.index,
                                   tau=tau,
                                   weights=weights)$fitted.values

      fitted.values <- tmp[,1]
      lwr <- tmp[,2]
      upr <- tmp[,3]
      rm(tmp)

      if(deriv > 0) {

        deriv.mat <- matrix(NA,nrow=NROW(newdata),ncol=NCOL(newdata))
        deriv.mat.lwr <- deriv.mat
        deriv.mat.upr <- deriv.mat
        l <- 1 ## num.z
        m <- 1 ## num.x
        for(i in 1:ncol(newdata)) {
          if(!is.factor(newdata[,i])) {
            if(deriv <= degree[m]) {
              tmp <- deriv.factor.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=K,
                                         I=include,
                                         xeval=xeval,
                                         zeval=zeval,
                                         knots=knots,
                                         basis=basis,
                                         deriv.index=m,
                                         deriv=deriv,
                                         prune.index=prune.index,
                                         tau=tau,
                                         weights=weights)
            } else {
              tmp <- matrix(0,nrow(xeval),3)
            }
            deriv.mat[,i] <- tmp[,1]
            deriv.mat.lwr[,i] <- tmp[,2]
            deriv.mat.upr[,i] <- tmp[,3]
            rm(tmp)
            m <- m + 1
          } else {
            zevaltmp <- zeval
            zevaltmp[,l] <- factor(rep(levels(newdata[,i])[1],NROW(newdata)),levels=levels(newdata[,i]),ordered=is.ordered(newdata[,i]))
            zpred <- predict.factor.spline(x=x,
                                           y=y,
                                           z=z,
                                           K=K,
                                           I=include,
                                           xeval=xeval,
                                           zeval=zeval,
                                           knots=knots,
                                           basis=basis,
                                           prune=prune,
                                           prune.index=prune.index,
                                           tau=tau,
                                           weights=weights)$fitted.values

            zpred.base <- predict.factor.spline(x=x,
                                                y=y,
                                                z=z,
                                                K=K,
                                                I=include,
                                                xeval=xeval,
                                                zeval=zevaltmp,
                                                knots=knots,
                                                basis=basis,
                                                prune=prune,
                                                prune.index=prune.index,
                                                tau=tau,
                                                weights=weights)$fitted.values

            deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
            deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
            deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

            l <- l + 1
          }
        }

      } else {
        deriv.mat <- NULL
        deriv.mat.lwr <- NULL
        deriv.mat.upr <- NULL
      }

    } else {

      ## Get degree vector and lambda vector

      complexity <- object$complexity
      knots <- object$knots
      K <- object$K
      segments <- object$segments
      degree <- object$degree
      lambda <- object$lambda

      is.ordered.z <- object$is.ordered.z

      z <- as.matrix(z)
      zeval <- as.matrix(zeval)

      tmp <- predict.kernel.spline(x=x,
                                   y=y,
                                   z=z,
                                   K=K,
                                   lambda=lambda,
                                   is.ordered.z=is.ordered.z,
                                   xeval=xeval,
                                   zeval=zeval,
                                   knots=knots,
                                   basis=basis,
                                   tau=tau,
                                   weights=weights)$fitted.values

      fitted.values <- tmp[,1]
      lwr <- tmp[,2]
      upr <- tmp[,3]
      rm(tmp)

      if(deriv > 0) {

        deriv.mat <- matrix(NA,nrow=NROW(newdata),ncol=NCOL(newdata))
        deriv.mat.lwr <- deriv.mat
        deriv.mat.upr <- deriv.mat
        l <- 1 ## num.z
        m <- 1 ## num.x
        for(i in 1:ncol(newdata)) {
          if(!is.factor(newdata[,i])) {
            if(deriv <= degree[m]) {
              tmp <- deriv.kernel.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=K,
                                         lambda=lambda,
                                         is.ordered.z=is.ordered.z,
                                         xeval=xeval,
                                         zeval=zeval,
                                         knots=knots,
                                         basis=basis,
                                         deriv.index=m,
                                         deriv=deriv,
                                         tau=tau,
                                         weights=weights)
            } else {
              tmp <- matrix(0,nrow(xeval),3)
            }
            deriv.mat[,i] <- tmp[,1]
            deriv.mat.lwr[,i] <- tmp[,2]
            deriv.mat.upr[,i] <- tmp[,3]
            rm(tmp)
            m <- m + 1
          } else {
            zevaltmp <- zeval
            zevaltmp[,l] <- rep(sort(unique(zeval[,l]))[1],NROW(zeval))

            zpred <- predict.kernel.spline(x=x,
                                           y=y,
                                           z=z,
                                           K=K,
                                           lambda=lambda,
                                           is.ordered.z=is.ordered.z,
                                           xeval=xeval,
                                           zeval=zeval,
                                           knots=knots,
                                           basis=basis,
                                           tau=tau,
                                           weights=weights)$fitted.values

            zpred.base <- predict.kernel.spline(x=x,
                                                y=y,
                                                z=z,
                                                K=K,
                                                lambda=lambda,
                                                is.ordered.z=is.ordered.z,
                                                xeval=xeval,
                                                zeval=zevaltmp,
                                                knots=knots,
                                                basis=basis,
                                                tau=tau,
                                                weights=weights)$fitted.values

            deriv.mat[,i] <- zpred[,1]-zpred.base[,1]
            deriv.mat.lwr[,i] <- deriv.mat[,i] - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
            deriv.mat.upr[,i] <- deriv.mat[,i] + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

            l <- l + 1
          }
        }

      } else {
        deriv.mat <- NULL
        deriv.mat.lwr <- NULL
        deriv.mat.upr <- NULL
      }

    }

  }

  ## Return the predicted values.

  attr(fitted.values, "lwr") <- lwr
  attr(fitted.values, "upr") <- upr
  attr(fitted.values, "deriv.mat") <- deriv.mat
  attr(fitted.values, "deriv.mat.lwr") <- deriv.mat.lwr
  attr(fitted.values, "deriv.mat.upr") <- deriv.mat.upr

  return(fitted.values)

}

## Basic print method.

print.crs <- function(x,
                      ...) {

  cat("Call:\n")
  print(x$call)

}

## print.summary is different from print.

summary.crs <- function(object,
                        sigtest=FALSE,
                        ...) {

  cat("Call:\n")
  print(object$call)
  if(!object$kernel) {
    if(is.null(object$tau))
      cat("\nIndicator Bases/B-spline Bases Regression Spline\n",sep="")
    else
      cat("\nIndicator Bases/B-spline Bases Quantile Regression Spline\n",sep="")
  } else {
    if(is.null(object$tau))
      cat("\nKernel Weighting/B-spline Bases Regression Spline\n",sep="")
    else
      cat("\nKernel Weighting/B-spline Bases Quantile Regression Spline\n",sep="")
  }
  if(!is.null(object$tau)) cat(paste("\nQuantile estimated: tau = ",format(object$tau),sep=""),sep="")
  if(object$num.x==1){
    cat(paste("\nThere is ",format(object$num.x), " continuous predictor",sep=""),sep="")
  } else {
    cat(paste("\nThere are ",format(object$num.x), " continuous predictors",sep=""),sep="")
  }
  if(!is.null(object$num.z)) if(object$num.z==1) {
    cat(paste("\nThere is ",format(object$num.z), " categorical predictor",sep=""),sep="")
  }  else {
    cat(paste("\nThere are ",format(object$num.z), " categorical predictors",sep=""),sep="")
  }
  for(j in 1:object$num.x)
    cat(paste("\nSpline degree/number of segments for ",format(object$xnames[j]),": ",format(object$degree[j]),"/",format(object$segments[j]),sep=""),sep="")
  if(!is.null(object$include)) for(j in 1:length(object$include))
    cat(paste("\nInclusion indicator for ",format(object$znames[j]),": ",format(object$include[j]),sep=""),sep="")
  if(!is.null(object$lambda)) for(j in 1:length(object$lambda))
    cat(paste("\nBandwidth for ",format(object$znames[j]),": ",format(object$lambda[j]),sep=""),sep="")
  cat(paste("\nModel complexity proxy: ", format(object$complexity), sep=""))
  cat(paste("\nKnot type: ", format(object$knots), sep=""))
  if(object$num.x > 1) cat(paste("\nBasis type: ",format(object$basis),sep=""))
  if(!object$kernel) cat(paste("\nPruning of final model: ",format(ifelse(object$prune,"TRUE","FALSE")),sep=""))
  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  if(is.null(object$tau)) cat(paste("\nRank of model frame: ", format(object$k), sep=""))
  cat(paste("\nTrace of smoother matrix: ", format(round(sum(object$hatvalues))), sep=""))

  if(is.null(object$weights)) {
    if(is.null(object$tau)) cat(paste("\n\nResidual standard error: ", format(sqrt(sum(object$residuals^2)/object$df.residual),digits=4)," on ", format(object$df.residual)," degrees of freedom",sep=""))
    adjusted.r.squared <- 1-(1-object$r.squared)*(length(object$fitted.values)-1)/object$df.residual
    if(is.null(object$tau)) cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4),",   Adjusted R-squared: ",format(adjusted.r.squared,digits=4), sep=""))
    df1 <- round(sum(object$hatvalues))-1
    df2 <- (object$nobs-round(sum(object$hatvalues)))
    F <- (df2/df1)*(sum((object$y-mean(object$y))^2)-sum(residuals(object)^2))/sum(residuals(object)^2)
    if(is.null(object$tau)) cat(paste("\nF-statistic: ", format(F,digits=4), " on ", df1, " and ", df2, " DF, p-value: ", format(pf(F,df1=df1,df2=df2,lower.tail=FALSE),digits=4), sep=""))
    if(!is.null(object$cv.score)) cat(paste("\n\nCross-validation score: ", format(object$cv.score,digits=8), sep=""))
  } else {
    if(is.null(object$tau)) cat(paste("\n\nResidual standard error (weighted): ", format(sqrt(sum((object$residuals^2)*object$weights)/object$df.residual),digits=4)," on ", format(object$df.residual)," degrees of freedom",sep=""))
    adjusted.r.squared <- 1-(1-object$r.squared)*(length(object$fitted.values)-1)/object$df.residual
    if(is.null(object$tau)) cat(paste("\nMultiple R-squared (weighted): ", format(object$r.squared,digits=4),",   Adjusted R-squared (weighted): ",format(adjusted.r.squared,digits=4), sep=""))
    df1 <- round(sum(object$hatvalues))-1
    df2 <- (object$nobs-round(sum(object$hatvalues)))
    F <- (df2/df1)*(sum((object$y-mean(object$y))^2*object$weights)-sum(residuals(object)^2*object$weights))/sum(residuals(object)^2*object$weights)
    if(is.null(object$tau)) cat(paste("\nF-statistic (weighted): ", format(F,digits=4), " on ", df1, " and ", df2, " DF, p-value: ", format(pf(F,df1=df1,df2=df2,lower.tail=FALSE),digits=4), sep=""))
    if(!is.null(object$cv.score)) cat(paste("\n\nCross-validation score (weighted): ", format(object$cv.score,digits=8), sep=""))
  }

  if(object$cv != "none") cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))

  if(sigtest&!object$kernel) {
    cat("\n\nPredictor significance test:\n")
    crs.sigtest(object)
  }

  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))

  cat("\n\n")

}

plot.crs <- function(x,
                     mean=FALSE,
                     deriv=0,
                     ci=FALSE,
                     num.eval=100,
                     caption=list("Residuals vs Fitted",
                       "Normal Q-Q Plot",
                       "Scale-Location",
                       "Cook's Distance"),
                     xtrim = 0.0,
                     xq = 0.5,
                     plot.behavior = c("plot","plot-data","data"),
                     common.scale=TRUE,
                     persp.rgl=FALSE,
                     ...) {

  plot.behavior <- match.arg(plot.behavior)

  object <- x

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  xq <- double(ncol(object$xz)) + xq

  ## Check for proper derivative

  if(deriv < 0) stop("derivative order must be a non-negative integer")

  ## Default - basic residual plots

  if(!mean&!deriv) {

    par(mfrow=c(2,2))

    ## Residuals versus fitted

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("\rWorking...",console = console)

    plot(fitted(object),
         residuals(object),
         xlab="Fitted Values",
         ylab="Residuals",
         main=caption[[1]],
         ...)

    ## QQ plot

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("\rWorking...",console = console)

    std.res <- residuals(object)/sqrt(mean(residuals(object)^2))

    qqnorm(std.res,
           main=caption[[2]])

    qqline(std.res)

    ## Standardized versus fitted

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("\rWorking...",console = console)

    plot(fitted(object),
         sqrt(abs(residuals(object,"pearson"))),
         xlab="Fitted Values",
         ylab=as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Standardized Residuals")))),
         main=caption[[3]],
         ...)

    ## Cook's distance \frac{\hat\epsilon_t^2
    ## h_{tt}}{\hat\sigma^2(1-h_{tt})^2k}. Note that this is not
    ## computed for kernel method (what is df.residual when there are
    ## multiple models, etc.?)


    if(!object$kernel) {

      console <- printClear(console)
      console <- printPop(console)
      console <- printPush("\rWorking...",console = console)

      sigmasq <- sum(residuals(object)^2)/object$df.residual
      cook <- (residuals(object)^2*object$hatvalues)/(sigmasq*(1-object$hatvalues)^2*object$k)
      plot(cook,
           type = "h",
           main=caption[[4]],
           xlab = "Observation Number",
           ylab = "Cook's Distance",
           ...)

    }

    console <- printClear(console)
    console <- printPop(console)

  }

  ## Mean

  if(mean) {

    ## Information required to compute predictions

    basis <- object$basis
    prune <- object$prune
    prune.index <- object$prune.index

    xz <- object$xz
    y <- object$y

    ## Divide into factors and numeric

    if(!object$kernel) {
      xztmp <- splitFrame(xz)
    } else {
      xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
    }
    x <- xztmp$x
    z <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z
    rm(xztmp)

    ## Get degree vector and lambda vector

    complexity <- object$complexity
    knots <- object$knots
    K <- object$K
    degree <- object$degree
    segments <- object$segments
    include <- object$include
    lambda <- object$lambda
    is.ordered.z <- object$is.ordered.z
    tau <- object$tau
    weights <- object$weights

    ## End information required to compute predictions

    if(!persp.rgl) {

      mg <- list()

      for(i in 1:NCOL(object$xz)) {

        if(!is.factor(object$xz[,i])) {
          newdata <- matrix(NA,nrow=num.eval,ncol=NCOL(object$xz))
          neval <- num.eval
        } else {
          newdata <- matrix(NA,nrow=length(levels(object$xz[,i])),ncol=NCOL(object$xz))
          neval <- length(levels(object$xz[,i]))
        }

        newdata <- data.frame(newdata)

        if(!is.factor(object$xz[,i])) {
          xlim <- trim.quantiles(object$xz[,i],xtrim)
          newdata[,i] <- seq(xlim[1],xlim[2],length=neval)
        } else {
          newdata[,i] <- factor(levels(object$xz[,i]),levels=levels(object$xz[,i]),ordered=is.ordered(object$xz[,i]))
        }

        for(j in (1:NCOL(object$xz))[-i]) {
          if(!is.factor(object$xz[,j])) {
            newdata[,j] <- rep(uocquantile(object$xz[,j],prob=xq[j]),neval)
          } else {
            newdata[,j] <- factor(rep(uocquantile(object$xz[,j],prob=xq[j]),neval),levels=levels(object$xz[,j]),ordered=is.ordered(object$xz[,j]))
          }
        }

        newdata <- data.frame(newdata)
        names(newdata) <- names(object$xz)

        if(!object$kernel) {
          xztmp <- splitFrame(data.frame(newdata))
        } else {
          xztmp <- splitFrame(data.frame(newdata),factor.to.numeric=TRUE)
        }
        xeval <- xztmp$x
        zeval <- xztmp$z
        is.ordered.z <- xztmp$is.ordered.z
        rm(xztmp)

        ## Compute the predicted values.

        if(!object$kernel) {

          tmp <- predict.factor.spline(x=x,
                                       y=y,
                                       z=z,
                                       K=K,
                                       I=include,
                                       xeval=xeval,
                                       zeval=zeval,
                                       basis=basis,
                                       knots=knots,
                                       prune=prune,
                                       prune.index=prune.index,
                                       tau=tau,
                                       weights=weights)$fitted.values

          fitted.values <- tmp[,1]
          lwr <- tmp[,2]
          upr <- tmp[,3]
          rm(tmp)

        } else {

          z <- as.matrix(z)
          zeval <- as.matrix(zeval)

          tmp <- predict.kernel.spline(x=x,
                                       y=y,
                                       z=z,
                                       K=K,
                                       lambda=lambda,
                                       is.ordered.z=is.ordered.z,
                                       xeval=xeval,
                                       zeval=zeval,
                                       knots=knots,
                                       basis=basis,
                                       tau=tau,
                                       weights=weights)$fitted.values

          fitted.values <- tmp[,1]
          lwr <- tmp[,2]
          upr <- tmp[,3]
          rm(tmp)

        }

        if(!ci) {

          mg[[i]] <- data.frame(newdata[,i],fitted.values)
          names(mg[[i]]) <- c(names(newdata)[i],"mean")

        } else {
          mg[[i]] <- data.frame(newdata[,i],fitted.values,lwr,upr)
          names(mg[[i]]) <- c(names(newdata)[i],"mean","lwr","upr")
        }

        console <- printClear(console)
        console <- printPop(console)

      }

      if(common.scale) {
        min.mg <- Inf
        max.mg <- -Inf
        for(i in 1:length(mg)) {
          min.mg <- min(min.mg,min(mg[[i]][,-1]))
          max.mg <- max(max.mg,max(mg[[i]][,-1]))
        }
        ylim <- c(min.mg,max.mg)
      } else {
        ylim <- NULL
      }

      if(plot.behavior!="data") {

        if(!is.null(object$num.z)||(object$num.x>1)) par(mfrow=dim.plot(NCOL(object$xz)))

        for(i in 1:NCOL(object$xz)) {

          if(!ci) {
            plot(mg[[i]][,1],mg[[i]][,2],
                 xlab=names(newdata)[i],
                 ylab=ifelse(is.null(tau),"Conditional Mean",paste("Conditional Quantile (tau = ",format(tau),")",sep="")),
                 ylim=ylim,
                 type="l",
                 ...)

          } else {
            if(!common.scale) ylim <- c(min(mg[[i]][,-1]),max(mg[[i]][,-1]))
            plot(mg[[i]][,1],mg[[i]][,2],
                 xlab=names(newdata)[i],
                 ylab=ifelse(is.null(tau),"Conditional Mean",paste("Conditional Quantile (tau = ",format(tau),")",sep="")),
                 ylim=ylim,
                 type="l",
                 ...)
            ## Need to overlay for proper plotting of factor errorbars
            par(new=TRUE)
            plot(mg[[i]][,1],mg[[i]][,3],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2,
                 ...)
            par(new=TRUE)
            plot(mg[[i]][,1],mg[[i]][,4],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2,
                 ...)
          }

        }

      }

    } else {


#      if(!require(rgl)) stop(" Error: you must first install the rgl package")

      if(!is.null(object$num.z)) stop(" Error: persp3d is for continuous predictors only")
      if(object$num.x != 2) stop(" Error: persp3d is for cases involving two continuous predictors only")

      newdata <- matrix(NA,nrow=num.eval,ncol=2)
      newdata <- data.frame(newdata)

      xlim <- trim.quantiles(object$xz[,1],xtrim)
      ylim <- trim.quantiles(object$xz[,2],xtrim)

      x1.seq <- seq(xlim[1],xlim[2],length=num.eval)
      x2.seq <- seq(ylim[1],ylim[2],length=num.eval)

      x.grid <- expand.grid(x1.seq,x2.seq)
      newdata <- data.frame(x.grid[,1],x.grid[,2])
      names(newdata) <- names(object$xz)

      z <- matrix(predict(object,newdata=newdata),num.eval,num.eval)

      mg <- list()

      mg[[1]] <- data.frame(newdata,z)

      if(plot.behavior!="data") {

        num.colors <- 1000
        colorlut <- topo.colors(num.colors)
        col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]

        open3d()

        par3d(windowRect=c(900,100,900+640,100+640))
        rgl.viewpoint(theta = 0, phi = -70, fov = 80)

        persp3d(x=x1.seq,y=x2.seq,z=z,
                xlab=names(object$xz)[1],ylab=names(object$xz)[2],zlab="Y",
                ticktype="detailed",
                border="red",
                color=col,
                alpha=.7,
                back="lines",
                main=ifelse(is.null(tau),"Conditional Mean",paste("Conditional Quantile (tau = ",format(tau),")",sep="")))

        grid3d(c("x", "y+", "z"))

        play3d(spin3d(axis=c(0,0,1), rpm=5), duration=15)

      }

    }

    if(plot.behavior!="plot") {
      console <- printClear(console)
      console <- printPop(console)
      if(!persp.rgl) par(mfrow=c(1,1))
      return(mg)
    }

  }

  ## deriv

  if(deriv > 0) {

    ## Information required to compute predictions

    basis <- object$basis
    prune <- object$prune
    prune.index <- object$prune.index

    xz <- object$xz
    y <- object$y

    ## Divide into factors and numeric

    if(!object$kernel) {
      xztmp <- splitFrame(xz)
    } else {
      xztmp <- splitFrame(xz,factor.to.numeric=TRUE)
    }
    x <- xztmp$x
    z <- xztmp$z
    is.ordered.z <- xztmp$is.ordered.z
    rm(xztmp)

    ## Get degree vector and lambda vector

    complexity <- object$complexity
    knots <- object$knots
    K <- object$K
    degree <- object$degree
    segments <- object$segments
    include <- object$include
    lambda <- object$lambda
    is.ordered.z <- object$is.ordered.z
    tau <- object$tau
    weights <- object$weights

    ## End information required to compute predictions

    if(deriv > 0) {

      rg <- list()
      m <- 0
      i.numeric <- 0

      for(i in 1:NCOL(object$xz)) {

        if(!is.factor(object$xz[,i])) {
          i.numeric <- i.numeric + 1
          newdata <- matrix(NA,nrow=num.eval,ncol=NCOL(object$xz))
          neval <- num.eval
          m <- m + 1
        } else {
          newdata <- matrix(NA,nrow=length(levels(object$xz[,i])),ncol=NCOL(object$xz))
          neval <- length(levels(object$xz[,i]))
        }

        newdata <- data.frame(newdata)
        newdata.base <- data.frame(newdata)

        if(!is.factor(object$xz[,i])) {
          xlim <- trim.quantiles(object$xz[,i],xtrim)
          newdata[,i] <- seq(xlim[1],xlim[2],length=neval)
        } else {
          newdata[,i] <- factor(levels(object$xz[,i]),levels=levels(object$xz[,i]),ordered=is.ordered(object$xz[,i]))
          newdata.base[,i] <- factor(rep(levels(object$xz[,i])[1],neval),levels=levels(object$xz[,i]),ordered=is.ordered(object$xz[,i]))
        }

        for(j in (1:NCOL(object$xz))[-i]) {
          if(!is.factor(object$xz[,j])) {
            newdata[,j] <- rep(uocquantile(object$xz[,j],prob=xq[j]),neval)
            newdata.base[,j] <- rep(uocquantile(object$xz[,j],prob=xq[j]),neval)
          } else {
            newdata[,j] <- factor(rep(uocquantile(object$xz[,j],prob=xq[j]),neval),levels=levels(object$xz[,j]),ordered=is.ordered(object$xz[,j]))
            newdata.base[,j] <- factor(rep(uocquantile(object$xz[,j],prob=xq[j]),neval),levels=levels(object$xz[,j]),ordered=is.ordered(object$xz[,j]))
          }
        }

        newdata <- data.frame(newdata)
        names(newdata) <- names(object$xz)
        newdata.base <- data.frame(newdata.base)
        names(newdata.base) <- names(object$xz)

        if(!object$kernel) {
          xztmp <- splitFrame(data.frame(newdata))
        } else {
          xztmp <- splitFrame(data.frame(newdata),factor.to.numeric=TRUE)
        }
        xeval <- xztmp$x
        zeval <- xztmp$z
        is.ordered.z <- xztmp$is.ordered.z
        rm(xztmp)

        if(!object$kernel) {
          xztmp <- splitFrame(data.frame(newdata.base))
        } else {
          xztmp <- splitFrame(data.frame(newdata.base),factor.to.numeric=TRUE)
        }
        xeval.base <- xztmp$x
        zeval.base <- xztmp$z
        is.ordered.z <- xztmp$is.ordered.z
        rm(xztmp)

        ## Compute the predicted values.

        if(!object$kernel) {

          if(!is.factor(newdata[,i])) {
            if(deriv <= degree[i.numeric]) {
              tmp <- deriv.factor.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=K,
                                         I=include,
                                         xeval=xeval,
                                         zeval=zeval,
                                         knots=knots,
                                         basis=basis,
                                         deriv.index=m,
                                         deriv=deriv,
                                         prune.index=prune.index,
                                         tau=tau,
                                         weights=weights)
            } else {
              tmp <- matrix(0,nrow(newdata),3)
            }
            deriv.est <- tmp[,1]
            deriv.lwr <- tmp[,2]
            deriv.upr <- tmp[,3]
            rm(tmp)
          } else {
            zpred <- predict.factor.spline(x=x,
                                           y=y,
                                           z=z,
                                           K=K,
                                           I=include,
                                           xeval=xeval,
                                           zeval=zeval,
                                           knots=knots,
                                           basis=basis,
                                           prune=prune,
                                           prune.index=prune.index,
                                           tau=tau,
                                           weights=weights)$fitted.values

            zpred.base <- predict.factor.spline(x=x,
                                                y=y,
                                                z=z,
                                                K=K,
                                                I=include,
                                                xeval=xeval.base,
                                                zeval=zeval.base,
                                                knots=knots,
                                                basis=basis,
                                                prune=prune,
                                                prune.index=prune.index,
                                                tau=tau,
                                                weights=weights)$fitted.values

            deriv.est <- zpred[,1]-zpred.base[,1]
            deriv.lwr <- deriv.est - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
            deriv.upr <- deriv.est + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

          }

        } else {

          if(!is.factor(newdata[,i])) {
            if(deriv <= degree[i.numeric]) {
              tmp <- deriv.kernel.spline(x=x,
                                         y=y,
                                         z=z,
                                         K=K,
                                         lambda=lambda,
                                         is.ordered.z=is.ordered.z,
                                         xeval=xeval,
                                         zeval=zeval,
                                         knots=knots,
                                         basis=basis,
                                         deriv.index=m,
                                         deriv=deriv,
                                         tau=tau,
                                         weights=weights)

            } else {
              tmp <- matrix(0,nrow(newdata),3)
            }
            deriv.est <- tmp[,1]
            deriv.lwr <- tmp[,2]
            deriv.upr <- tmp[,3]
            rm(tmp)
          } else {
            zpred <- predict.kernel.spline(x=x,
                                           y=y,
                                           z=z,
                                           K=K,
                                           lambda=lambda,
                                           is.ordered.z=is.ordered.z,
                                           xeval=xeval,
                                           zeval=zeval,
                                           knots=knots,
                                           basis=basis,
                                           tau=tau,
                                           weights=weights)$fitted.values

            zpred.base <- predict.kernel.spline(x=x,
                                                y=y,
                                                z=z,
                                                K=K,
                                                lambda=lambda,
                                                is.ordered.z=is.ordered.z,
                                                xeval=xeval.base,
                                                zeval=zeval.base,
                                                knots=knots,
                                                basis=basis,
                                                tau=tau,
                                                weights=weights)$fitted.values

            deriv.est <- zpred[,1]-zpred.base[,1]
            deriv.lwr <- deriv.est - qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)
            deriv.upr <- deriv.est + qnorm(0.975)*sqrt(zpred[,4]^2+zpred.base[,4]^2)

          }

        }

        if(!ci) {

          rg[[i]] <- data.frame(newdata[,i],deriv.est)
          names(rg[[i]]) <- c(names(newdata)[i],"deriv")

        } else {

          rg[[i]] <- data.frame(newdata[,i],
                                deriv.est,
                                deriv.lwr,
                                deriv.upr)
          names(rg[[i]]) <- c(names(newdata)[i],"deriv","lwr","upr")

        }

        console <- printClear(console)
        console <- printPop(console)

      }

    }

    if(deriv > 0) {

      if(common.scale) {
        min.rg <- Inf
        max.rg <- -Inf
        for(i in 1:length(rg)) {
          min.rg <- min(min.rg,min(rg[[i]][,-1]))
          max.rg <- max(max.rg,max(rg[[i]][,-1]))
        }
        ylim <- c(min.rg,max.rg)
      } else {
        ylim <- NULL
      }

      if(plot.behavior!="data") {

        if(!is.null(object$num.z)||(object$num.x>1)) par(mfrow=dim.plot(NCOL(object$xz)))

        for(i in 1:NCOL(object$xz)) {

          if(!ci) {
            plot(rg[[i]][,1],rg[[i]][,2],
                 xlab=names(newdata)[i],
                 ylab=ifelse(!is.factor(newdata[,i]), paste("Order", deriv,"Derivative"), "Difference in Levels"),
                 ylim=ylim,
                 type="l",
                 ...)

          } else {
            if(!common.scale) ylim <- c(min(rg[[i]][,-1]),max(rg[[i]][,-1]))
            plot(rg[[i]][,1],rg[[i]][,2],
                 xlab=names(newdata)[i],
                 ylab=ifelse(!is.factor(newdata[,i]), paste("Order", deriv,"Derivative"), "Difference in Levels"),
                 ylim=ylim,
                 type="l",
                 ...)
            ## Need to overlay for proper plotting of factor errorbars
            par(new=TRUE)
            plot(rg[[i]][,1],rg[[i]][,3],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2,
                 ...)
            par(new=TRUE)
            plot(rg[[i]][,1],rg[[i]][,4],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2,
                 ...)
          }

        }

      }

      if(plot.behavior!="plot") {
        console <- printClear(console)
        console <- printPop(console)
        par(mfrow=c(1,1))
        return(rg)
      }

    }

  }

  ## Reset par to 1,1 (can be modified above)

  console <- printClear(console)
  console <- printPop(console)

}

crs.sigtest <- function(object,...) {

  ## This function for the asymptotic significance test can be
  ## airlifted in trivially... trace of the smoother matrix etc. will
  ## deliver correct F stat etc. Left for future 9/1/11 since we have
  ## crssigtest function...

  if(object$kernel) stop(" sigtest is currently available only when kernel=FALSE")

  sg <- list()

  ## Conduct the significance test in order variable by variable

  j.num.x <- 1
  j.num.z <- 1

  for(i in 1:NCOL(object$xz)) {

    if(!is.factor(object$xz[,i])) {
      degree <- object$degree
      degree[j.num.x] <- 0
      model.res <- crs(object$formula,cv="none",degree=degree,include=object$include,basis=object$basis,prune=object$prune,data=eval(object$call$data))
      sg[[i]] <- anova(model.res$model.lm,object$model.lm)
      j.num.x <- j.num.x + 1
    } else {
      include <- object$include
      include[j.num.z] <- 0
      model.res <- crs(object$formula,cv="none",degree=object$degree,include=include,basis=object$basis,prune=object$prune,data=eval(object$call$data))
      sg[[i]] <- anova(model.res$model.lm,object$model.lm)
      j.num.z <- j.num.z + 1
    }

    cat(paste("Predictor ", format(names(object$xz)[i]), ": Df = ", sg[[i]]$Df[2], ", F = ", format(sg[[i]]$F[2],digits=4), ", Pr(>F) = ", format(sg[[i]][[6]][2],digits=4), "\n", sep=""))

  }

}
