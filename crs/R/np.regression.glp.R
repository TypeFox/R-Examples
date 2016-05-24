## This file contains code for multivariate generalized local
## polynomial kernel regression with mixed datatypes. It relies on
## functions in the np package and on snomadr which currently resides
## in the crs package (June 29 2011).

## Note that the approach taken here is computationally efficient and
## relies on expressing the local polynomial method in slightly
## different form purely for computational simplicity. Both approaches
## are identical though.

scale.robust <- function(y){
 if(any(dim(as.matrix(y)) == 0))
      return(0)
  sd.vec <- apply(as.matrix(y),2,sd)
  IQR.vec <- apply(as.matrix(y),2,IQR)/(qnorm(.25,lower.tail=F)*2)
  mad.vec <- apply(as.matrix(y),2,mad)
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) max(x))
  if(any(a<=0)) warning(paste("variable ",which(a<=0)," appears to be constant",sep=""))
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) min(x[x>0]))  
  return(a)
}

mypoly <- function(x,
                   ex=NULL,
                   degree,
                   gradient.compute = FALSE,
                   r=0,
                   Bernstein = TRUE) {

  if(missing(x)) stop(" Error: x required")
  if(missing(degree)) stop(" Error: degree required")
  if(degree < 1) stop(" Error: degree must be a positive integer")
  if(!is.logical(Bernstein)) stop(" Error: Bernstein must be logical")

  if(!Bernstein) {

    ## Raw polynomials and their derivatives

    if(!is.null(ex)) x <- ex

    if(gradient.compute) {
      Z <- NULL
      for(i in 1:degree) {
        if((i-r) >= 0) {
          tmp <- (factorial(i)/factorial(i-r))*x^max(0,i-r)
        } else {
          tmp <- rep(0,length(x))
        }
        Z <- cbind(Z,tmp)
      }
    } else {
      Z <- outer(x,1L:degree,"^")
    }

  } else {
    ## Bernstein polynomials and their derivatives (i.e. Bezier curves
    ## i.e. B-splines with no interior knots)
    if(is.null(ex)) {
      if(gradient.compute) {
        Z <- gsl.bs(x,degree=degree,deriv=r)
      } else {
        Z <- gsl.bs(x,degree=degree)
      }
    } else {
      if(gradient.compute) {
        Z <- predict(gsl.bs(x,degree=degree,deriv=r),newx=ex)
      } else {
        Z <- predict(gsl.bs(x,degree=degree),newx=ex)
      }
    }
  }

  return(as.matrix(Z))

}

## W.glp is a modified version of the polym() function (stats). The
## function accepts a vector of degrees and provides a generalized
## polynomial with varying polynomial order.

W.glp <- function(xdat = NULL,
                  exdat = NULL,
                  degree = NULL,
                  gradient.vec = NULL,
                  Bernstein = TRUE) {

  if(is.null(xdat)) stop(" Error: You must provide data")
  if(is.null(degree) || any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  xdat <- as.data.frame(xdat)

  xdat.col.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})
  k <- ncol(as.data.frame(xdat[,xdat.col.numeric]))

  xdat.numeric <- NULL
  if(k > 0) {
    xdat.numeric <- as.data.frame(xdat[,xdat.col.numeric])
    if(!is.null(exdat)) {
      exdat.numeric <- as.data.frame(exdat[,xdat.col.numeric])
    } else {
      exdat.numeric <- NULL
    }
  }

  if(!is.null(gradient.vec) && (length(gradient.vec) != k)) stop(paste(" Error: gradient vector and number of numeric predictors must be conformable\n",sep=""))
  if(!is.null(gradient.vec) && any(gradient.vec < 0)) stop(paste(" Error: gradient vector must contain non-negative integers\n",sep=""))

  if(!is.null(gradient.vec)) {
    gradient.compute <- TRUE
  } else {
    gradient.compute <- FALSE
    gradient.vec <- rep(NA,k)
  }

  if(length(degree) != k) stop(" Error: degree vector and number of numeric predictors incompatible")

  if(all(degree == 0) || (k == 0)) {

    ## Local constant OR no continuous variables

    if(is.null(exdat)) {
      return(matrix(1,nrow=nrow(xdat.numeric),ncol=1))
    } else {
      return(matrix(1,nrow=nrow(exdat.numeric),ncol=1))
    }

  } else {

    degree.list <- list()
    for(i in 1:k) degree.list[[i]] <- 0:degree[i]
    z <- do.call("expand.grid", degree.list, k)
    s <- rowSums(z)
    ind <- (s > 0) & (s <= max(degree))
    z <- z[ind, ,drop=FALSE]
    if(!all(degree==max(degree))) {
      for(j in 1:length(degree)) {
        d <- degree[j]
        if((d < max(degree)) & (d > 0)) {
          s <- rowSums(z)
          d <- (s > d) & (z[,j,drop=FALSE]==matrix(d,nrow(z),1,byrow=TRUE))
          z <- z[!d, ]
        }
      }
    }
    if(is.null(exdat)) {
      res <- rep.int(1,nrow(xdat.numeric))
    } else {
      res <- rep.int(1,nrow(exdat.numeric))
    }
    res.deriv <- 1
    if(degree[1] > 0) {
      res <- cbind(1, mypoly(x=xdat.numeric[,1],
                             ex=exdat.numeric[,1],
                             degree=degree[1],
                             gradient.compute=gradient.compute,
                             r=gradient.vec[1],
                             Bernstein=Bernstein))[, 1 + z[, 1]]

      if(gradient.compute && gradient.vec[1] != 0) res.deriv <- cbind(1,matrix(NA,1,degree[1]))[, 1 + z[, 1],drop=FALSE]
      if(gradient.compute && gradient.vec[1] == 0) res.deriv <- cbind(1,matrix(0,1,degree[1]))[, 1 + z[, 1],drop=FALSE]
    }
    if(k > 1) for (i in 2:k) if(degree[i] > 0) {
      res <- res * cbind(1, mypoly(x=xdat.numeric[,i],
                                   ex=exdat.numeric[,i],
                                   degree=degree[i],
                                   gradient.compute=gradient.compute,
                                   r=gradient.vec[i],
                                   Bernstein=Bernstein))[, 1 + z[, i]]
      if(gradient.compute && gradient.vec[i] != 0) res.deriv <- res.deriv * cbind(1,matrix(NA,1,degree[i]))[, 1 + z[, i],drop=FALSE]
      if(gradient.compute && gradient.vec[i] == 0) res.deriv <- res.deriv *cbind(1,matrix(0,1,degree[i]))[, 1 + z[, i],drop=FALSE]
    }

    if(is.null(exdat)) {
      res <- matrix(res,nrow=NROW(xdat))
    } else {
      res <- matrix(res,nrow=NROW(exdat))
    }
    if(gradient.compute) res.deriv <- matrix(res.deriv,nrow=1)
    colnames(res) <- apply(z, 1L, function(x) paste(x, collapse = "."))
    if(gradient.compute) colnames(res.deriv) <- apply(z, 1L, function(x) paste(x, collapse = "."))

    if(gradient.compute) {
      res[,!is.na(as.numeric(res.deriv))] <- 0
      return(cbind(0,res))
    } else {
      return(cbind(1,res))
    }

  }

}

## This function determines the maximum value for k for
## k-nearest-neighbor-based estimation. This is necessary as data can
## contain repeated values or contain whole numbers both of which
## reduce the maximum value of k but this cannot be determined without
## actually computing all distances. This function will be numerically
## intensive for large datasets, however, it ensures that search for
## instance would take place over the largest permissible range. It
## gets called one time for each numeric predictor. It is redundant
## for draws from a continuous distribution.

## Since the distances are symmetric there are naturally savings to be
## achieved (currently we take a brute force approach).

## Note this can suffer from comparison of floating point numbers and
## needs to be improved (re-scaling may not be kind).

knn.max <- function(x) {

  k.max <- length(x)-1
  non.unique <- length(unique(x)) != length(x)
  if(non.unique) x <- unique(x)
  for(i in 1:length(x)) {
    x.diff.unique <- sort(unique(abs(x[i]-x[-i])))
    if(length(x.diff.unique) < k.max)
      if(length(x.diff.unique)-1 < k.max) k.max <- length(x.diff.unique)-1
  }

  return(k.max)

}

## This function will check whether the polynomial for a given
## predictor is ill conditioned, return TRUE if it is but also pass
## (as an attribute) a vector containing the maximum well-conditioned
## polynomial degree for each numeric predictor in xdat

check.max.degree <- function(xdat=NULL,degree=NULL,issue.warning=FALSE,Bernstein=TRUE) {

  if(is.null(xdat)) stop(" xdat must be provided")
  if(is.null(degree)) stop(" degree vector must be provided")

  xdat <- as.data.frame(xdat)

  ill.conditioned <- FALSE

  xdat.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})
  numeric.index <- which(xdat.numeric==TRUE)
  num.numeric <- sum(sapply(1:NCOL(xdat),function(i){is.numeric(xdat[,i])})==TRUE)
  d <- numeric(num.numeric)

  if(num.numeric > 0) {

    for(i in 1:num.numeric) {
      if(degree[i]>0) {
        ## First check if degree results in a well-conditioned basis
        d[i] <- degree[i]
        X <- mypoly(x=xdat[,numeric.index[i]],
                    ex=NULL,
                    degree=degree[i],
                    Bernstein=Bernstein)
        ## If the max degree is ill-conditioned then ascend
        ## (descending may require computation of large matrices that
        ## are discarded)
        if(!is.fullrank(X)) {
          for(j in 1:degree[i]) {
            d[i] <- j
            X <- mypoly(x=xdat[,numeric.index[i]],
                        ex=NULL,
                        degree=d[i],
                        Bernstein=Bernstein)
            if(!is.fullrank(X)) {
              d[i] <- j-1
              break()
            }
          }
        }
        if(d[i] < degree[i]) {
          if(issue.warning) warning(paste("\r Predictor ",i," polynomial is ill-conditioned beyond degree ",d[i],": see note in ?npglpreg",sep=""))
          ill.conditioned <- TRUE
        }
      }
    }

  }

  attr(ill.conditioned, "degree.max.vec") <- d
  return(ill.conditioned)

}

npglpreg <- function(...) UseMethod("npglpreg")

npglpreg.default <- function(tydat=NULL,
                             txdat=NULL,
                             eydat=NULL,
                             exdat=NULL,
                             bws=NULL,
                             degree=NULL,
                             leave.one.out=FALSE,
                             ckertype=c("gaussian", "epanechnikov","uniform","truncated gaussian"),
                             ckerorder=2,
                             ukertype=c("liracine","aitchisonaitken"),
                             okertype=c("liracine","wangvanryzin"),
                             bwtype = c("fixed","generalized_nn","adaptive_nn","auto"),
                             gradient.vec=NULL,
                             gradient.categorical=FALSE,
                             cv.shrink=TRUE,
                             cv.maxPenalty=sqrt(.Machine$double.xmax),
                             cv.warning=FALSE,
                             Bernstein=TRUE,
                             mpi=FALSE,
                             ...) {

  ckertype <- match.arg(ckertype)
  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)
  if(!any(ckerorder==c(2,4,6,8))) stop("ckerorder must be 2, 4, 6, or 8")

  ### Nov. 23, issue is that if you don't shrink for estimation you
  ### can barf, so we only set shrink=FALSE when conducting CV. This
  ### occurs because glpcvNOMAD uses the leave-one-out estimator, here
  ### we use full sample, so it is possible for the leave-one-out to
  ### have fullrank but the full sample to fail. Thus we need to use
  ### ridging for estimation hence override the value here.

  est <- glpregEst(tydat=tydat,
                   txdat=txdat,
                   eydat=eydat,
                   exdat=exdat,
                   bws=bws,
                   degree=degree,
                   leave.one.out=leave.one.out,
                   ckertype=ckertype,
                   ckerorder=ckerorder,
                   ukertype=ukertype,
                   okertype=okertype,
                   bwtype=bwtype,
                   gradient.vec=gradient.vec,
                   cv.shrink=TRUE, ## Override
                   cv.warning=cv.warning,
                   Bernstein=Bernstein,
                   ...)

  est$cv.shrink <- cv.shrink

  ## Gradients for categorical predictors (here we just compute them
  ## all though gradient for continuous predictors given above
  ## requires specific variables to be selected)

  est$gradient.categorical.mat <- NULL

  if(gradient.categorical) {

    num.numeric <- est$num.numeric
    num.categorical <- est$num.categorical

    if(num.categorical > 0) {

      num.eval <- ifelse(is.null(exdat),nrow(txdat),nrow(exdat))
      gradient.categorical.mat <- matrix(NA,nrow=num.eval,ncol=num.categorical)

      for(i in 1:num.categorical) {

        if(is.null(exdat)) {
          exdat.base <- txdat
        } else {
          exdat.base <- exdat
        }

        categorical.index <- est$categorical.index
        eval.base <- levels(txdat[,categorical.index[i]])[1]
        eval.levels <- levels(txdat[,categorical.index[i]])

        if(is.ordered(txdat[,categorical.index[i]])) {
          exdat.base[,categorical.index[i]] <- ordered(rep(eval.base,num.eval),levels=eval.levels)
        } else {
          exdat.base[,categorical.index[i]] <- factor(rep(eval.base,num.eval),levels=eval.levels)
        }

        est.base <- glpregEst(tydat=tydat,
                              txdat=txdat,
                              eydat=eydat,
                              exdat=exdat.base,
                              bws=bws,
                              degree=degree,
                              leave.one.out=leave.one.out,
                              ckertype=ckertype,
                              ckerorder=ckerorder,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,
                              gradient.vec=NULL,
                              cv.shrink=TRUE, ## Override
                              cv.warning=cv.warning,
                              Bernstein=Bernstein,
                              ...)

        gradient.categorical.mat[,i] <- est$fitted.values - est.base$fitted.values

      }

      est$gradient.categorical.mat <- gradient.categorical.mat

    }

  }

  ## Add results to estimated object.

  if(!is.null(eydat)) {
    est$r.squared <- RSQfunc(eydat,est$fitted.values)
    est$residuals <- eydat - est$fitted.values
  } else if(is.null(eydat)&&is.null(exdat)) {
    est$r.squared <- RSQfunc(tydat,est$fitted.values)
    est$residuals <- tydat - est$fitted.values
  } else {
    est$r.squared <- NULL
    est$residuals <- NULL
  }

  est$call <- match.call()

  ## Return object of type npglpreg

  class(est) <- "npglpreg"

  return(est)

}

## Basic print method.

print.npglpreg <- function(x,
                           ...) {
  cat("Call:\n")
  print(x$call)

}

summary.npglpreg <- function(object,
                             ...) {

  cat("Call:\n")
  print(object$call)
  cat("\nGeneralized Local Polynomial Kernel Regression\n",sep="")

  ## Summarize continuous predictors

  if(object$Bernstein) {
    cat("\nPolynomial type: Bernstein")
  } else {
    cat("\nPolynomial type: raw")
  }

  if(object$cv.shrink)
    cat("\nUsing (local) Seifert & Gasser shrinkage for cross-validation")

  if(object$num.numeric == 1){
    cat(paste("\nThere is ",format(object$num.numeric), " continuous predictor",sep=""),sep="")
  } else if(object$num.numeric > 1) {
    cat(paste("\nThere are ",format(object$num.numeric), " continuous predictors",sep=""),sep="")
  }

  cat(paste("\nBandwidth type: ", object$bwtype, sep=""))

  if(object$num.numeric >= 1) {
    cat(paste("\nContinuous kernel type: ", object$ckertype, sep=""))
    cat(paste("\nContinuous kernel order: ", object$ckerorder, sep=""))
    for(j in 1:object$num.numeric) {
      if(object$bwtype=="fixed") {
        cat(paste("\nBandwidth for ",format(object$xnames[object$numeric.index][j]),": ",format(object$bws[object$numeric.index][j]),sep=""),sep="")
        if(!is.null(object$bws.sf))
          cat(paste(" (scale factor = ", format(object$bws.sf[object$numeric.index][j]),")",sep=""),sep="")
      } else {
        cat(paste("\nKth nearest neighbor for ",format(object$xnames[object$numeric.index][j]),": ",format(object$bws[object$numeric.index][j]),sep=""),sep="")
      }
    }
  }

  for(j in 1:object$num.numeric)
    cat(paste("\nDegree for ",format(object$xnames[object$numeric.index][j]),": ",format(object$degree[j]),sep=""),sep="")

  ## Summarize categorical predictors

  if(object$num.categorical==1) {
    cat(paste("\nThere is ",format(object$num.categorical), " categorical predictor",sep=""),sep="")
  } else if(object$num.categorical > 1) {
    cat(paste("\nThere are ",format(object$num.categorical), " categorical predictors",sep=""),sep="")
  }

  if(object$num.categorical >= 1) {
    cat(paste("\nUnordered kernel type: ", object$ukertype, sep=""))
    cat(paste("\nOrdered kernel type: ", object$okertype, sep=""))
    for(j in 1:(object$num.categorical))
      cat(paste("\nBandwidth for ",format(object$xnames[object$categorical.index][j]),": ",format(object$bws[object$categorical.index][j]),sep=""),sep="")
  }

  ## Summary statistics

  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4), sep=""))
  if(!is.null(object$fv)) {
    cat(paste("\nCross-validation score: ", format(object$fv,digits=8), sep=""))
    cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
  }
  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))
  cat("\n\n")

}

## Method for predicting given a new data frame.

predict.npglpreg <- function(object,
                             newdata=NULL,
                             gradient.vec=NULL,
                             ...) {

  if(is.null(newdata)) {

    ## If no new data provided, return sample fit.
    fitted.values <- fitted(object)
    gradient <- object$gradient
    gradient.categorical.mat <- object$gradient.categorical.mat

  } else{

    ## Get training data from object (xz and y) and parse into factors
    ## and numeric.

    degree <- object$degree
    bws <- object$bws
    bwtype <- object$bwtype
    ckertype <- object$ckertype
    ckerorder <- object$ckerorder
    ukertype <- object$ukertype
    okertype <- object$okertype
    Bernstein <- object$Bernstein
    mpi <- object$mpi
    if(is.null(gradient.vec)) {
      gradient.vec <- object$gradient.vec
    }

    txdat <- object$x
    tydat <- object$y

    tt <- terms(object)
    has.ey <- succeedWithResponse(tt, newdata)
    if (has.ey) {
      eydat <- model.response(model.frame(tt,newdata))
    } else {
      eydat <- NULL
    }
    exdat <- model.frame(delete.response(tt),newdata,xlev=object$xlevels)

    ## Return the predicted values.

    est <- npglpreg.default(tydat=tydat,
                            txdat=txdat,
                            exdat=exdat,
                            eydat=eydat,
                            bws=bws,
                            degree=degree,
                            ckertype=ckertype,
                            ckerorder=ckerorder,
                            ukertype=ukertype,
                            okertype=okertype,
                            bwtype=bwtype,
                            gradient.vec=gradient.vec,
                            Bernstein=Bernstein,
                            mpi=mpi,
                            ...)

    fitted.values <- est$fitted.values
    gradient <- est$gradient
    gradient.categorical.mat <- est$gradient.categorical.mat

  }

  attr(fitted.values, "gradient") <- gradient
  attr(fitted.values, "gradient.categorical.mat") <- gradient.categorical.mat

  return(fitted.values)

}

## Note that to exploit NOMAD and scaling this function and the
## functions it calls operate on the scaling factors and not the raw
## bandwidths to avoid issues of scale associated with variable
## measurement. In effect we need to put the parameters on which NOMAD
## operates on the same `scale'. I achieve this by letting the
## bandwidth scale factors lie in bandwidth.min to bandwidth.max then
## rescale by this factor. This is the computational breakthrough I
## have been missing and it took many frustrating days but would not
## have occurred without some inspiration. But so far it is working
## brilliantly and much better than before. All other functions
## operate as one would expect.

npglpreg.formula <- function(formula,
                             data=list(),
                             tydat=NULL,
                             txdat=NULL,
                             eydat=NULL,
                             exdat=NULL,
                             bws=NULL,
                             degree=NULL,
                             leave.one.out=FALSE,
                             ckertype=c("gaussian", "epanechnikov","uniform","truncated gaussian"),
                             ckerorder=2,
                             ukertype=c("liracine","aitchisonaitken"),
                             okertype=c("liracine","wangvanryzin"),
                             bwtype = c("fixed","generalized_nn","adaptive_nn","auto"),
                             cv=c("degree-bandwidth","bandwidth","none"),
                             cv.func=c("cv.ls","cv.aic"),
                             nmulti=5,
                             random.seed=42,
                             degree.max=10,
                             degree.min=0,
                             bandwidth.max=.Machine$double.xmax,
                             bandwidth.min=sqrt(.Machine$double.eps),
                             bandwidth.min.numeric=1.0e-02,
                             bandwidth.switch=1.0e+06,
                             bandwidth.scale.categorical=1.0e+04,
                             max.bb.eval=10000,
                             min.epsilon=.Machine$double.eps,
                             initial.mesh.size.real=1,
                             initial.mesh.size.integer=1,
                             min.mesh.size.real=sqrt(.Machine$double.eps),
                             min.mesh.size.integer=sqrt(.Machine$double.eps),
                             min.poll.size.real=sqrt(.Machine$double.eps),
                             min.poll.size.integer=sqrt(.Machine$double.eps),
                             opts=list(),
                             restart.from.min=FALSE,
                             gradient.vec=NULL,
                             gradient.categorical=FALSE,
                             cv.shrink=TRUE,
                             cv.maxPenalty=sqrt(.Machine$double.xmax),
                             cv.warning=FALSE,
                             Bernstein=TRUE,
                             mpi=FALSE,
                             ...) {

  ## Basic error trapping...

#  if(cv.func=="cv.aic" && cv.shrink==TRUE) {
#      warning("cv.shrink and cv.aic currently incompatible, cv.shrink set to FALSE")
#      cv.shrink <- FALSE
#  }

  if(!is.logical(mpi)) stop(" Error: mpi must be logical (TRUE/FALSE)")
  if(!is.logical(Bernstein)) stop(" Error: Bernstein must be logical (TRUE/FALSE)")
  if(!is.logical(gradient.categorical)) stop(" Error: gradient.categorical must be logical (TRUE/FALSE)")
  if(!is.logical(cv.warning)) stop(" Error: cv.warning must be logical (TRUE/FALSE)")
  if(!is.logical(leave.one.out)) stop(" Error: leave.one.out must be logical (TRUE/FALSE)")
  if(degree.max > 100) stop(paste(" degree.max (",degree.max,") exceeds reasonable value (",100,")",sep=""))
  if(degree.max < 1) stop(paste(" degree.max (",degree.max,") must be a positive integer",sep=""))
  if(degree.max < degree.min) stop(" degree.max must be greater than or equal to degree.min")
  if(as.numeric(initial.mesh.size.real) <= 0) stop(" initial.mesh.size.real must be positive")
  if(as.numeric(initial.mesh.size.integer) <= 0) stop(" initial.mesh.size.integer must be positive")
  if(as.numeric(min.mesh.size.integer) <= 0) stop(" min.mesh.size.integer must be positive")
  if(as.numeric(min.mesh.size.real) <= 0) stop(" min.mesh.size.real must be positive")
  if(as.numeric(min.poll.size.integer) <= 0) stop(" min.poll.size.integer must be positive")
  if(as.numeric(min.poll.size.real) <= 0) stop(" min.poll.size.real must be positive")
  if(as.numeric(bandwidth.scale.categorical) <= 0) stop(" bandwidth.scale.categorical must be positive")
  if(as.numeric(bandwidth.min) <= 0) stop(" bandwidth.min must be positive")
  if(as.numeric(bandwidth.min.numeric) <= 0) stop(" bandwidth.min.numeric must be positive")
  if(as.numeric(bandwidth.switch) <= 0) stop(" bandwidth.switch must be positive")
  if(as.numeric(bandwidth.max) <= 0) stop(" bandwidth.max must be positive")
  if(as.numeric(bandwidth.max) <= as.numeric(bandwidth.min)) stop(" bandwidth.max must exceed bandwidth.min")
  if(as.numeric(min.epsilon) <= 0) stop(" min.epsilon must be positive")
  if(as.numeric(min.epsilon) >= as.numeric(min.mesh.size.real)) stop(" min.epsilon must be less than min.mesh.size.real")
  if(as.numeric(min.epsilon) >= as.numeric(min.mesh.size.integer)) stop(" min.epsilon must be less than min.mesh.size.integer")
  if(as.numeric(min.epsilon) >= as.numeric(min.poll.size.real)) stop(" min.epsilon must be less than min.poll.size.real")
  if(as.numeric(min.epsilon) >= as.numeric(min.poll.size.integer)) stop(" min.epsilon must be less than min.poll.size.integer")
  if(as.numeric(max.bb.eval) <= 0) stop(" max.bb.eval must be positive")
#  if(!mpi) {
#    if(!require(np)) stop(" Error: you must install the np package to use this function")
#  } else {
#    if(!require(npRmpi)) stop(" Error: you must install the npRmpi package to use this function")
#  }
  ## Set DISPLAY_DEGREE to 0 if crs.messages=FALSE and DISPLAY_DEGREE
  ## is not provided

  if(!options('crs.messages')$crs.messages && is.null(opts[["DISPLAY_DEGREE"]])) opts$"DISPLAY_DEGREE"=0

  ckertype <- match.arg(ckertype)
  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)
  if(!any(ckerorder==c(2,4,6,8))) stop("ckerorder must be 2, 4, 6, or 8")

  cv <- match.arg(cv)
  cv.func <- match.arg(cv.func)

  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  tydat <- model.response(mf)
  txdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

  fv <- NULL
  ptm <- system.time("")

  bws.sf <- NULL

  if(cv=="none"&&bwtype=="auto") stop(" Error: you cannot use bwtype==\"auto\" without running cross-validation")

  if(cv!="none"&&bwtype!="auto") warning(paste(" bwtype is ", bwtype, ": you could consider bwtype=\"auto\"",sep=""),immediate.=TRUE)

  if(cv!="none") {
    if(bwtype!="auto") {
      ptm <- ptm + system.time(model.cv <-glpcvNOMAD(ydat=tydat,
                                                     xdat=txdat,
                                                     cv=cv,
                                                     degree=degree,
                                                     bandwidth=bws,
                                                     bwmethod=cv.func,
                                                     ckertype=ckertype,
                                                     ckerorder=ckerorder,
                                                     ukertype=ukertype,
                                                     okertype=okertype,
                                                     bwtype=bwtype,
                                                     nmulti=nmulti,
                                                     random.seed=random.seed,
                                                     degree.max=degree.max,
                                                     degree.min=degree.min,
                                                     bandwidth.max=bandwidth.max,
                                                     bandwidth.min=bandwidth.min,
                                                     bandwidth.min.numeric=bandwidth.min.numeric,
                                                     bandwidth.switch=bandwidth.switch,
                                                     bandwidth.scale.categorical=bandwidth.scale.categorical,
                                                     max.bb.eval=max.bb.eval,
                                                     min.epsilon=min.epsilon,
                                                     initial.mesh.size.real=initial.mesh.size.real,
                                                     initial.mesh.size.integer=initial.mesh.size.integer,
                                                     min.mesh.size.real=min.mesh.size.real,
                                                     min.mesh.size.integer=min.mesh.size.integer,
                                                     min.poll.size.real=min.poll.size.real,
                                                     min.poll.size.integer=min.poll.size.integer,
                                                     opts=opts,
                                                     restart.from.min=restart.from.min,
                                                     cv.shrink=cv.shrink,
                                                     cv.maxPenalty=cv.maxPenalty,
                                                     cv.warning=cv.warning,
                                                     Bernstein=Bernstein,
                                                     mpi=mpi,
                                                     ...))
    } else {
      ptm <- ptm + system.time(model.cv <-glpcvNOMAD(ydat=tydat,
                                                     xdat=txdat,
                                                     cv=cv,
                                                     degree=degree,
                                                     bandwidth=bws,
                                                     bwmethod=cv.func,
                                                     ckertype=ckertype,
                                                     ckerorder=ckerorder,
                                                     ukertype=ukertype,
                                                     okertype=okertype,
                                                     bwtype="fixed",
                                                     nmulti=nmulti,
                                                     random.seed=random.seed,
                                                     degree.max=degree.max,
                                                     degree.min=degree.min,
                                                     bandwidth.max=bandwidth.max,
                                                     bandwidth.min=bandwidth.min,
                                                     bandwidth.min.numeric=bandwidth.min.numeric,
                                                     bandwidth.switch=bandwidth.switch,
                                                     bandwidth.scale.categorical=bandwidth.scale.categorical,
                                                     max.bb.eval=max.bb.eval,
                                                     min.epsilon=min.epsilon,
                                                     initial.mesh.size.real=initial.mesh.size.real,
                                                     initial.mesh.size.integer=initial.mesh.size.integer,
                                                     min.mesh.size.real=min.mesh.size.real,
                                                     min.mesh.size.integer=min.mesh.size.integer,
                                                     min.poll.size.real=min.poll.size.real,
                                                     min.poll.size.integer=min.poll.size.integer,
                                                     opts=opts,
                                                     restart.from.min=restart.from.min,
                                                     cv.shrink=cv.shrink,
                                                     cv.maxPenalty=cv.maxPenalty,
                                                     cv.warning=cv.warning,
                                                     Bernstein=Bernstein,
                                                     mpi=mpi,
                                                     ...))

      bwtype <- "fixed"
      fv <- model.cv$fv

      ptm <- ptm + system.time(model <-glpcvNOMAD(ydat=tydat,
                                                  xdat=txdat,
                                                  cv=cv,
                                                  degree=degree,
                                                  bandwidth=bws,
                                                  bwmethod=cv.func,
                                                  ckertype=ckertype,
                                                  ckerorder=ckerorder,
                                                  ukertype=ukertype,
                                                  okertype=okertype,
                                                  bwtype="generalized_nn",
                                                  nmulti=nmulti,
                                                  random.seed=random.seed,
                                                  degree.max=degree.max,
                                                  degree.min=degree.min,
                                                  bandwidth.max=bandwidth.max,
                                                  bandwidth.min=bandwidth.min,
                                                  bandwidth.min.numeric=bandwidth.min.numeric,
                                                  bandwidth.switch=bandwidth.switch,
                                                  bandwidth.scale.categorical=bandwidth.scale.categorical,
                                                  max.bb.eval=max.bb.eval,
                                                  min.epsilon=min.epsilon,
                                                  initial.mesh.size.real=initial.mesh.size.real,
                                                  initial.mesh.size.integer=initial.mesh.size.integer,
                                                  min.mesh.size.real=min.mesh.size.real,
                                                  min.mesh.size.integer=min.mesh.size.integer,
                                                  min.poll.size.real=min.poll.size.real,
                                                  min.poll.size.integer=min.poll.size.integer,
                                                  opts=opts,
                                                  restart.from.min=restart.from.min,
                                                  cv.shrink=cv.shrink,
                                                  cv.maxPenalty=cv.maxPenalty,
                                                  cv.warning=cv.warning,
                                                  Bernstein=Bernstein,
                                                  mpi=mpi,
                                                  ...))

      if(model$fv < model.cv$fv) {
        model.cv <- model
        bwtype <- "generalized_nn"
      }

      ptm <- ptm + system.time(model <-glpcvNOMAD(ydat=tydat,
                                                  xdat=txdat,
                                                  cv=cv,
                                                  degree=degree,
                                                  bandwidth=bws,
                                                  bwmethod=cv.func,
                                                  ckertype=ckertype,
                                                  ckerorder=ckerorder,
                                                  ukertype=ukertype,
                                                  okertype=okertype,
                                                  bwtype="adaptive_nn",
                                                  nmulti=nmulti,
                                                  random.seed=random.seed,
                                                  degree.max=degree.max,
                                                  degree.min=degree.min,
                                                  bandwidth.max=bandwidth.max,
                                                  bandwidth.min=bandwidth.min,
                                                  bandwidth.min.numeric=bandwidth.min.numeric,
                                                  bandwidth.switch=bandwidth.switch,
                                                  bandwidth.scale.categorical=bandwidth.scale.categorical,
                                                  max.bb.eval=max.bb.eval,
                                                  min.epsilon=min.epsilon,
                                                  initial.mesh.size.real=initial.mesh.size.real,
                                                  initial.mesh.size.integer=initial.mesh.size.integer,
                                                  min.mesh.size.real=min.mesh.size.real,
                                                  min.mesh.size.integer=min.mesh.size.integer,
                                                  min.poll.size.real=min.poll.size.real,
                                                  min.poll.size.integer=min.poll.size.integer,
                                                  opts=opts,
                                                  restart.from.min=restart.from.min,
                                                  cv.shrink=cv.shrink,
                                                  cv.maxPenalty=cv.maxPenalty,
                                                  cv.warning=cv.warning,
                                                  Bernstein=Bernstein,
                                                  mpi=mpi,
                                                  ...))
      if(model$fv < model.cv$fv) {
        model.cv <- model
        bwtype <- "adaptive_nn"
      }

    }

    degree <- model.cv$degree
    bws <- model.cv$bws
    bws.sf <- model.cv$bws.sf
    fv <- model.cv$fv
    if(isTRUE(all.equal(fv,cv.maxPenalty))) stop(" Search failed: restart with larger nmulti or smaller degree.max (or degree if provided)")
  }

  if(!is.null(degree))    {
    ill.conditioned <- check.max.degree(txdat,degree,issue.warning=TRUE,Bernstein=Bernstein)
    degree.max.vec <- attr(ill.conditioned, "degree.max.vec")
    degree <- ifelse(degree > degree.max.vec, degree.max.vec, degree)
  }

  ptm <- ptm + system.time(est <- npglpreg.default(tydat=tydat,
                                                   txdat=txdat,
                                                   eydat=eydat,
                                                   exdat=exdat,
                                                   bws=bws,
                                                   degree=degree,
                                                   leave.one.out=leave.one.out,
                                                   ckertype=ckertype,
                                                   ckerorder=ckerorder,
                                                   ukertype=ukertype,
                                                   okertype=okertype,
                                                   bwtype=bwtype,
                                                   gradient.vec=gradient.vec,
                                                   gradient.categorical=gradient.categorical,
                                                   cv.shrink=TRUE, ## Must override cv.shrink here
                                                   cv.maxPenalty=cv.maxPenalty,
                                                   cv.warning=cv.warning,
                                                   Bernstein=Bernstein,
                                                   mpi=mpi,
                                                   ...))

  est$call <- match.call()
  est$formula <- formula
  est$terms <- mt
  est$xlevels <- .getXlevels(mt, mf)
  est$x <- txdat
  est$y <- tydat
  est$fv <- fv
  est$nmulti <- nmulti
  est$ptm <- ptm
  est$bws.sf <- bws.sf

  return(est)

}

glpregEst <- function(tydat=NULL,
                      txdat=NULL,
                      eydat=NULL,
                      exdat=NULL,
                      bws=NULL,
                      degree=NULL,
                      leave.one.out=FALSE,
                      ckertype=c("gaussian", "epanechnikov","uniform","truncated gaussian"),
                      ckerorder=2,
                      ukertype=c("liracine","aitchisonaitken"),
                      okertype=c("liracine","wangvanryzin"),
                      bwtype=c("fixed","generalized_nn","adaptive_nn"),
                      gradient.vec=NULL,
                      cv.shrink=TRUE,
                      cv.maxPenalty=sqrt(.Machine$double.xmax),
                      cv.warning=FALSE,
                      Bernstein=TRUE,
                      ...) {

  ckertype <- match.arg(ckertype)
  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)
  if(!any(ckerorder==c(2,4,6,8))) stop("ckerorder must be 2, 4, 6, or 8")

  if(is.null(tydat)) stop(" Error: You must provide y data")
  if(is.null(txdat)) stop(" Error: You must provide X data")
  if(is.null(bws)) stop(" Error: You must provide a bandwidth object")
  if(is.null(degree) || any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  miss.ex = is.null(exdat)
  miss.ey = is.null(eydat)

  if (miss.ex){
    exdat <- txdat
  }

  txdat <- as.data.frame(txdat)
  exdat <- as.data.frame(exdat)

  n.train <- nrow(txdat)
  n.eval <- nrow(exdat)

  xdat.numeric <- sapply(1:ncol(txdat),function(i){is.numeric(txdat[,i])})
  categorical.index <- which(xdat.numeric==FALSE)
  numeric.index <- which(xdat.numeric==TRUE)
  num.numeric <- sum(sapply(1:NCOL(txdat),function(i){is.numeric(txdat[,i])})==TRUE)
  num.categorical <- NCOL(txdat)-num.numeric

  if(num.numeric == 0) stop("generalized local polynomial regression requires at least one numeric predictor")

  ## Test for invalid knn values

  ## Below is the worst case scenario where
  ## all(txdat[,numeric.index[i]]%%1==0). Otherwise, we would have to
  ## compute all distances and then take the smallest integer for which
  ## the distance is defined.

  if(bwtype!="fixed" && !is.null(bws) && num.numeric > 0) {
    for(i in 1:num.numeric) {
      k.max <- knn.max(txdat[,numeric.index[i]])
      if(bws[numeric.index[i]] > k.max) stop(paste("Error: invalid knn provided for predictor ",numeric.index[i],": max is ",k.max,sep=""))
    }
  }

  ## Check whether it appears that training and evaluation data are
  ## conformable

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  if(ncol(txdat)!=ncol(exdat))
    stop(" Error: training and evaluation data have unequal number of columns\n")

  if(all(degree == 0)) {

    ## Local constant using only one call to npksum

    if(leave.one.out == TRUE) {

      ## exdat not supported with leave.one.out, but this is only used
      ## for cross-validation hence no exdat

      tww <- npksum(txdat = txdat,
                    weights = as.matrix(data.frame(1,tydat)),
                    tydat = rep(1,length(tydat)),
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ckertype = ckertype,
                    ckerorder=ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

    } else {

      tww <- npksum(txdat = txdat,
                    exdat = exdat,
                    weights = as.matrix(data.frame(1,tydat)),
                    tydat = rep(1,length(tydat)),
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ckertype = ckertype,
                    ckerorder = ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

    }

    ## Note that as bandwidth approaches zero the local constant
    ## estimator undersmooths and approaches each sample realization,
    ## so use the convention that when the sum of the kernel weights
    ## equals 0, return y. This is unique to this code.

    mhat <- tww[2,]/NZD(tww[1,])

    return(list(fitted.values = mhat,
                gradient = NULL,
                coef.mat = NULL,
                bwtype = bwtype,
                ckertype = ckertype,
                ckerorder=ckerorder,
                ukertype = ukertype,
                okertype = okertype,
                degree = degree,
                bws = bws,
                nobs = n.train,
                num.numeric = num.numeric,
                num.categorical = num.categorical,
                xnames = names(txdat),
                categorical.index = categorical.index,
                numeric.index = numeric.index,
                gradient.vec = gradient.vec,
                cv.shrink = cv.shrink,
                Bernstein = Bernstein))

  } else {

    ## Test for negative degrees of freedom

    if(dim.bs(basis="glp",kernel=TRUE,degree=degree,segments=rep(1,length(degree)))>n.train-1)
      stop(" Ill-conditioned polynomial basis encountered: modify polynomial order")

    W <- W.glp(xdat=txdat,
               degree=degree,
               Bernstein=Bernstein)

    ## Test for ill-conditioned polynomial basis

    if(!is.fullrank(W))
      stop(" Ill-conditioned polynomial basis encountered: modify polynomial order")

    if(miss.ex) {
      W.eval <- W

      if(!is.null(gradient.vec)) {
        W.eval.deriv <- W.glp(xdat=txdat,
                              degree=degree,
                              gradient.vec=gradient.vec,
                              Bernstein=Bernstein)
      }
    } else {
      W.eval <- W.glp(xdat=txdat,
                      exdat=exdat,
                      degree=degree,
                      Bernstein=Bernstein)

      if(!is.null(gradient.vec)) {
        W.eval.deriv <- W.glp(xdat=txdat,
                              exdat=exdat,
                              degree=degree,
                              gradient.vec=gradient.vec,
                              Bernstein=Bernstein)
      }
    }

    ## Local polynomial via smooth coefficient formulation and one
    ## call to npksum

    if(leave.one.out == TRUE) {

      ## exdat not supported with leave.one.out, but this is only used
      ## for cross-validation hence no exdat

      tww <- npksum(txdat = txdat,
                    tydat = as.matrix(cbind(tydat,W)),
                    weights = W,
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ckertype = ckertype,
                    ckerorder=ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

    } else {

      tww <- npksum(txdat = txdat,
                    exdat = exdat,
                    tydat = as.matrix(cbind(tydat,W)),
                    weights = W,
                    bws = bws,
                    bandwidth.divide = TRUE,
                    leave.one.out = leave.one.out,
                    ckertype = ckertype,
                    ckerorder=ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

    }

    tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n.eval))[1,,]
    tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n.eval))[-1,,]

    if(!cv.shrink) {
      ## We can choose to use ridging or simply check for less than
      ## full column rank. If we test for full column rank ridging
      ## ought never occur but this is quite strong as if one
      ## observation only is sparse and inversion is problematic, for
      ## the rest of the sample this may not be the case.
      for(i in 1:n.eval) {
        if(!is.fullrank(tww[,,i])) {
          if(cv.warning) console <- printPush(paste("\rWarning: is.fullrank required for inversion at obs. ", i," failed      ",sep=""),console = console)
          return(cv.maxPenalty)
        }
      }
    }

    ## This traps the case with one evaluation point where we need to
    ## keep the extra dimension.

    if(!is.matrix(tyw)) {
      tyw <- matrix(tyw)
      tww <- array(tww,dim = c(ncol(W),ncol(W),n.eval))
    }

    coef.mat <- matrix(cv.maxPenalty,ncol(W),n.eval)
    epsilon <- 1.0/n.eval
    ridge <- double(n.eval)
    ridge.lc <- double(n.eval)    
    doridge <- !logical(n.eval)

    nc <- ncol(tww[,,1])

    ridger <- function(i) {
      doridge[i] <<- FALSE
      ridge.lc[i] <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
      tryCatch(chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc))))%*%tyw[,i],
               error = function(e){
                 ridge[i] <<- ridge[i]+epsilon
                 doridge[i] <<- TRUE
                 if(cv.warning) console <- printPush(paste("\rWarning: ridging required for inversion at obs. ", i, ", ridge = ",formatC(ridge[i],digits=4,format="f"),"        ",sep=""),console = console)
                 return(rep(cv.maxPenalty,nc))
               })
    }

    ## Test for singularity of the generalized local polynomial
    ## estimator, cv.shrink the mean towards the local constant mean.

    while(any(doridge)){
      iloo <- (1:n.eval)[doridge]
      coef.mat[,iloo] <- sapply(iloo, ridger)
    }

    ## Shrinking towards the local constant mean is accomplished via
    ## ridge.lc[i] which is ridge[i] times the local constant
    ## estimator

    mhat <- sapply(1:n.eval, function(i) {
      (1-ridge[i]) * W.eval[i,, drop = FALSE] %*% coef.mat[,i] + ridge.lc[i]
    })

    ## Ought to have a correction here for derivative when shrinking
    ## towards the local constant - XXX

    if(!is.null(gradient.vec)) {
      gradient <- sapply(1:n.eval, function(i) {
        W.eval.deriv[i,-1, drop = FALSE] %*% coef.mat[-1,i]
      })
    } else {
      gradient <- NULL
    }

    return(list(fitted.values = mhat,
                gradient = gradient,
                coef.mat = t(coef.mat[-1,]),
                bwtype = bwtype,
                ckertype = ckertype,
                ckerorder=ckerorder,
                ukertype = ukertype,
                okertype = okertype,
                degree = degree,
                bws = bws,
                nobs = n.train,
                num.numeric = num.numeric,
                num.categorical = num.categorical,
                xnames = names(txdat),
                categorical.index = categorical.index,
                numeric.index = numeric.index,
                gradient.vec = gradient.vec,
                cv.shrink = cv.shrink,
                Bernstein = Bernstein))

  }

}

## Create the function wrappers to be fed to the snomadr solver for
## leave-one-out cross-validation and Hurvich, Simonoff, and Tsai's
## AIC_c approach

minimand.cv.ls <- function(bws=NULL,
                           ydat=NULL,
                           xdat=NULL,
                           degree=NULL,
                           W=NULL,
                           ckertype=c("gaussian", "epanechnikov","uniform","truncated gaussian"),
                           ckerorder=2,
                           ukertype=c("liracine","aitchisonaitken"),
                           okertype=c("liracine","wangvanryzin"),
                           bwtype=c("fixed","generalized_nn","adaptive_nn"),
                           cv.shrink=TRUE,
                           cv.maxPenalty=sqrt(.Machine$double.xmax),
                           cv.warning=FALSE,
                           bandwidth.scale.categorical=NULL,
                           ...) {

  ckertype <- match.arg(ckertype)
  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)
  if(!any(ckerorder==c(2,4,6,8))) stop("ckerorder must be 2, 4, 6, or 8")

  if(is.null(ydat)) stop(" Error: You must provide y data")
  if(is.null(xdat)) stop(" Error: You must provide X data")
  if(is.null(W)) stop(" Error: You must provide a weighting matrix W")
  if(is.null(bws)) stop(" Error: You must provide a bandwidth object")
  if(is.null(degree) || any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  xdat <- as.data.frame(xdat)

  n <- length(ydat)

  ## Manually conduct bandwidth scaling so that NOMAD is operating on
  ## a scale free level.

  num.bw <- ncol(xdat)
  xdat.numeric <- sapply(1:num.bw,function(i){is.numeric(xdat[,i])})
  num.numeric <- ncol(as.data.frame(xdat[,xdat.numeric]))

  for(i in 1:num.bw) {
    if(xdat.numeric[i]==TRUE && bwtype=="fixed") {
      bws[i] <- bws[i]*scale.robust(xdat[,i])*length(ydat)^{-1/(num.numeric+2*ckerorder)}
    }
    if(xdat.numeric[i]!=TRUE) {
      bws[i] <- bws[i]/bandwidth.scale.categorical
    }
  }

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  if(!is.null(W)) {
    ## Check for positive degrees of freedom
    if(ncol(W) >= nrow(W)-1) {
      if(cv.warning) console <- printPush(paste("\rWarning: negative degrees of freedom                           ",sep=""),console = console)
      return(cv.maxPenalty)
    }
    ## Check for full column rank
    if(!is.fullrank(W)) {
      if(cv.warning) console <- printPush(paste("\rWarning: negative degrees of freedom                           ",sep=""),console = console)
      return(cv.maxPenalty)
    }
  }

  if(any(bws<=0)) {

    return(cv.maxPenalty)

  } else {

    if(all(degree == 0)) {

      ## Local constant via one call to npksum

      tww <- npksum(txdat = xdat,
                    weights = as.matrix(data.frame(1,ydat)),
                    tydat = rep(1,n),
                    bws = bws,
                    leave.one.out = TRUE,
                    bandwidth.divide = TRUE,
                    ckertype = ckertype,
                    ckerorder=ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

      mean.loo <- tww[2,]/NZD(tww[1,])

      if (!any(mean.loo == cv.maxPenalty)){
        fv <- mean((ydat-mean.loo)^2)
      } else {
        fv <- cv.maxPenalty
      }

      fv <- ifelse(is.finite(fv),fv,cv.maxPenalty)

      console <- printPush("\r                                                                         ",console = console)
      console <- printPush(paste("\rfv = ",format(fv)," ",sep=""),console = console)

      return(fv)

    } else {

      ## Generalized local polynomial via smooth coefficient
      ## formulation and one call to npksum

      tww <- npksum(txdat = xdat,
                    tydat = as.matrix(cbind(ydat,W)),
                    weights = W,
                    bws = bws,
                    leave.one.out = TRUE,
                    bandwidth.divide = TRUE,
                    ckertype = ckertype,
                    ckerorder = ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

      tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[1,,]
      tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[-1,,]

      if(!cv.shrink) {
        ## We can choose to use ridging or simply check for less than
        ## full column rank. If we test for full column rank ridging
        ## ought never occur but this is quite strong as if one
        ## observation only is sparse and inversion is problematic, for
        ## the rest of the sample this may not be the case.
        for(i in 1:n) {
          if(!is.fullrank(tww[,,i])) {
            if(cv.warning) console <- printPush(paste("\rWarning: is.fullrank required for inversion at obs. ", i," failed      ",sep=""),console = console)
            return(cv.maxPenalty)
          }
        }
      }

      mean.loo <- rep(cv.maxPenalty,n)
      epsilon <- 1.0/n
      ridge <- double(n)
      ridge.lc <- double(n)      
      doridge <- !logical(n)

      nc <- ncol(tww[,,1])

      ## Test for singularity of the generalized local polynomial
      ## estimator, shrink the mean towards the local constant mean.

      ridger <- function(i) {
        doridge[i] <<- FALSE
        ridge.lc[i] <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
        W[i,, drop = FALSE] %*% tryCatch(chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc))))%*%tyw[,i],
                error = function(e){
                  ridge[i] <<- ridge[i]+epsilon
                  doridge[i] <<- TRUE
                  if(cv.warning) console <- printPush(paste("\rWarning: ridging required for inversion at obs. ", i, ", ridge = ",formatC(ridge[i],digits=4,format="f"),"        ",sep=""),console = console)
                  return(rep(cv.maxPenalty,nc))
                })
      }

      ## Shrinking towards the local constant mean is accomplished via
      ## ridge.lc[i] which is ridge[i] times the local constant
      ## estimator

      while(any(doridge)){
        iloo <- (1:n)[doridge]
        mean.loo[iloo] <- (1-ridge[i])*sapply(iloo, ridger) + ridge.lc[i]
      }

      if (!is.na(any(mean.loo == cv.maxPenalty)) && !any(mean.loo == cv.maxPenalty)){
        fv <- mean((ydat-mean.loo)^2)
      } else {
        fv <- cv.maxPenalty
      }

      fv <- ifelse(is.finite(fv),fv,cv.maxPenalty)

      console <- printPush("\r                                                                         ",console = console)
      console <- printPush(paste("\rfv = ",format(fv)," ",sep=""),console = console)

      return(fv)

    }

  }

}

minimand.cv.aic <- function(bws=NULL,
                            ydat=NULL,
                            xdat=NULL,
                            degree=NULL,
                            W=NULL,
                            ckertype=c("gaussian", "epanechnikov","uniform","truncated gaussian"),
                            ckerorder=2,
                            ukertype=c("liracine","aitchisonaitken"),
                            okertype=c("liracine","wangvanryzin"),
                            bwtype = c("fixed","generalized_nn","adaptive_nn"),
                            cv.shrink=TRUE,
                            cv.maxPenalty=sqrt(.Machine$double.xmax),
                            cv.warning=FALSE,
                            bandwidth.scale.categorical=NULL,
                            ...) {

  ckertype <- match.arg(ckertype)
  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)
  if(!any(ckerorder==c(2,4,6,8))) stop("ckerorder must be 2, 4, 6, or 8")

  if(is.null(ydat)) stop(" Error: You must provide y data")
  if(is.null(xdat)) stop(" Error: You must provide X data")
  if(!all(degree==0)) if(is.null(W)) stop(" Error: You must provide a weighting matrix W")
  if(is.null(bws)) stop(" Error: You must provide a bandwidth object")
  if(is.null(degree) || any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  xdat <- as.data.frame(xdat)

  n <- length(ydat)

  if(!is.null(W)) {
    ## Check for positive degrees of freedom
    if(ncol(W) >= nrow(W)-1) {
      if(cv.warning) console <- printPush(paste("\rWarning: negative degrees of freedom                           ",sep=""),console = console)
      return(cv.maxPenalty)
    }
    ## Check for full column rank
    if(!is.fullrank(W)) {
      if(cv.warning) console <- printPush(paste("\rWarning: is.fullrank required for inversion at obs. ", i," failed      ",sep=""),console = console)
      return(cv.maxPenalty)
    }
  }

  ## Manually conduct bandwidth scaling so that NOMAD is operating on
  ## a scale free level.

  num.bw <- ncol(xdat)
  xdat.numeric <- sapply(1:num.bw,function(i){is.numeric(xdat[,i])})
  num.numeric <- ncol(as.data.frame(xdat[,xdat.numeric]))

  for(i in 1:num.bw) {
    if(xdat.numeric[i]==TRUE && bwtype=="fixed") {
      bws[i] <- bws[i]*scale.robust(xdat[,i])*length(ydat)^{-1/(num.numeric+2*ckerorder)}
    }
    if(xdat.numeric[i]!=TRUE) {
      bws[i] <- bws[i]/bandwidth.scale.categorical
    }
  }

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)

  if(any(bws<=0)) {

    return(cv.maxPenalty)

  } else {

    ## This computes the kernel function when i=j (i.e., K(0))

    kernel.i.eq.j <- npksum(txdat = xdat[1,],
                            weights = as.matrix(data.frame(1,ydat)[1,]),
                            tydat = 1,
                            bws = bws,
                            bandwidth.divide = TRUE,
                            ckertype = ckertype,
                            ckerorder=ckerorder,
                            ukertype = ukertype,
                            okertype = okertype,
                            bwtype = bwtype,
                            ...)$ksum[1,1]

    if(all(degree == 0)) {

      ## Local constant via one call to npksum

      tww <- npksum(txdat = xdat,
                    weights = as.matrix(data.frame(1,ydat)),
                    tydat = rep(1,n),
                    bws = bws,
                    bandwidth.divide = TRUE,
                    ckertype = ckertype,
                    ckerorder=ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

      ghat <- tww[2,]/NZD(tww[1,])

      trH <- kernel.i.eq.j*sum(1/NZD(tww[1,]))

      aic.penalty <- (1+trH/n)/(1-(trH+2)/n)

      if (!any(ghat == cv.maxPenalty) && (aic.penalty > 0)){
        fv <- log(mean((ydat-ghat)^2)) + aic.penalty
      } else {
        fv <- cv.maxPenalty
      }

      return(ifelse(is.finite(fv),fv,cv.maxPenalty))

    } else {

      ## Generalized local polynomial via smooth coefficient
      ## formulation and one call to npksum

      tww <- npksum(txdat = xdat,
                    tydat = as.matrix(cbind(ydat,W)),
                    weights = W,
                    bws = bws,
                    bandwidth.divide = TRUE,
                    ckertype = ckertype,
                    ckerorder=ckerorder,
                    ukertype = ukertype,
                    okertype = okertype,
                    bwtype = bwtype,
                    ...)$ksum

      tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[1,,]
      tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[-1,,]

      if(!cv.shrink) {
        ## We can choose to use ridging or simply check for less than
        ## full column rank. If we test for full column rank ridging
        ## ought never occur but this is quite strong as if one
        ## observation only is sparse and inversion is problematic, for
        ## the rest of the sample this may not be the case.
        for(i in 1:n) {
          if(!is.fullrank(tww[,,i])) {
            if(cv.warning) console <- printPush(paste("\rWarning: is.fullrank required for inversion at obs. ", i," failed      ",sep=""),console = console)
            return(cv.maxPenalty)
          }
        }
      }

      ghat <- rep(cv.maxPenalty,n)
      epsilon <- 1.0/n
      ridge <- double(n)
      ridge.lc <- double(n)      
      doridge <- !logical(n)

      nc <- ncol(tww[,,1])

      ## Test for singularity of the generalized local polynomial
      ## estimator, shrink the mean towards the local constant mean.

      ridger <- function(i) {
        doridge[i] <<- FALSE
        ridge.lc[i] <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
        W[i,, drop = FALSE] %*% tryCatch(chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc))))%*%tyw[,i],
                error = function(e){
                  ridge[i] <<- ridge[i]+epsilon
                  doridge[i] <<- TRUE
                  if(cv.warning) console <- printPush(paste("\rWarning: ridging required for inversion at obs. ", i, ", ridge = ",formatC(ridge[i],digits=4,format="f"),"        ",sep=""),console = console)
                  return(rep(cv.maxPenalty,nc))
                })
      }

      ## Shrinking towards the local constant mean is accomplished via
      ## ridge.lc[i] which is ridge[i] times the local constant
      ## estimator

      while(any(doridge)){
        ii <- (1:n)[doridge]
        ghat[ii] <- (1-ridge[i])*sapply(ii, ridger) + ridge.lc[i]
      }

      trH <- kernel.i.eq.j*sum(sapply(1:n,function(i){
        (1-ridge[i]) * W[i,, drop = FALSE] %*% chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc)))) %*% t(W[i,, drop = FALSE]) + ridge[i]/NZD(tww[,,i][1,1])
      }))

      aic.penalty <- (1+trH/n)/(1-(trH+2)/n)

      if (!any(ghat == cv.maxPenalty) && (aic.penalty > 0)){
        fv <- log(mean((ydat-ghat)^2)) + aic.penalty
      } else {
        fv <- cv.maxPenalty
      }

      fv <- ifelse(is.finite(fv),fv,cv.maxPenalty)

      console <- printPush("\r                                                                         ",console = console)
      console <- printPush(paste("\rfv = ",format(fv)," ",sep=""),console = console)

      return(fv)

    }

  }

}

## Note that for the function below which uses NOMAD for optimization
## we cross-validate on scale factors throughout (i.e. bandwidths for
## continuous predictors are scaled according to
## scale.robust()*length(ydat)^{-1/(num.numeric+2*ckerorder)}). We adjust
## the upper bounds for the categorical variables accordingly (i.e. 1
## and (c-1)/c). This is done so that the grid search takes place on a
## _somewhat_ common scale and also allows us to sidestep the issue
## with upper bounds where we were previously using relative tolerance
## (I am grateful for Sebastien Le Digabel for providing hints that
## allowed me to overcome this issue). This appears to resolve a
## longstanding issue caused by the use of relative tolerance and
## bounds and appears to put npglpreg() on equal footing with npreg()
## in the np package. But npglpreg() may in fact scale better with n
## while it admits generalized polynomials that npreg() lacks.

glpcvNOMAD <- function(ydat=NULL,
                       xdat=NULL,
                       degree=NULL,
                       bandwidth=NULL,
                       bwmethod=c("cv.ls","cv.aic"),
                       ckertype=c("gaussian", "epanechnikov","uniform","truncated gaussian"),
                       ckerorder=2,
                       ukertype=c("liracine","aitchisonaitken"),
                       okertype=c("liracine","wangvanryzin"),
                       bwtype = c("fixed","generalized_nn","adaptive_nn"),
                       cv=c("degree-bandwidth", "bandwidth"),
                       nmulti=NULL,
                       random.seed=42,
                       degree.max=10,
                       degree.min=0,
                       bandwidth.max=.Machine$double.xmax,
                       bandwidth.min=sqrt(.Machine$double.eps),
                       bandwidth.min.numeric=1.0e-02,
                       bandwidth.switch=1.0e+06,
                       bandwidth.scale.categorical=1.0e+04,
                       max.bb.eval=10000,
                       min.epsilon=.Machine$double.eps,
                       initial.mesh.size.real=1,
                       initial.mesh.size.integer=1,
                       min.mesh.size.real=sqrt(.Machine$double.eps),
                       min.mesh.size.integer=sqrt(.Machine$double.eps),
                       min.poll.size.real=sqrt(.Machine$double.eps),
                       min.poll.size.integer=sqrt(.Machine$double.eps),
                       cv.shrink=TRUE,
                       cv.maxPenalty=sqrt(.Machine$double.xmax),
                       cv.warning=FALSE,
                       Bernstein=TRUE,
                       mpi=FALSE,
                       restart.from.min=FALSE,
                       opts=list(),
                       ...) {

  ## Functions which, if defined outside, cause R warnings currently
  ## (Jan 19 2014, reported to Zhenghua)

  ## Jan 19 2014, optional args error reported to zhenghua... want to
  ## pass ... but having ... in eval.lscv and in snomadr causes an error
  ## message thrown by checkFunctionArguments
  
  ## In snomadr() we write wrappers around user-defined functions to
  ## pass additional arguments
  ##  eval.f.wrapper <- function(x){ eval.f(x,...) }
  
  eval.lscv <- function(input,
                        params){
  
    ydat <- params$ydat
    xdat <- params$xdat
    xdat.numeric <- params$xdat.numeric
    num.bw <- params$num.bw
    num.numeric <- params$num.numeric
    cv.maxPenalty <- params$cv.maxPenalty
    degree <- params$degree
    cv <- params$cv
    ckertype <- params$ckertype
    ckerorder <- params$ckerorder
    ukertype <- params$ukertype
    okertype <- params$okertype
    bwtype <- params$bwtype
    cv.shrink <- params$cv.shrink
    cv.warning <- params$cv.warning
    Bernstein <- params$Bernstein
    bw.switch <- params$bw.switch
    bandwidth.scale.categorical=params$bandwidth.scale.categorical
  
    bw.gamma <- input[1:num.bw]
    if(cv=="degree-bandwidth")
      degree <- round(input[(num.bw+1):(num.bw+num.numeric)])
  
    ## Test for negative degrees of freedom
  
    if(dim.bs(basis="glp",kernel=TRUE,degree=degree,segments=rep(1,length(degree)))>length(ydat)-1)
      return(cv.maxPenalty)
  
    W <- W.glp(xdat=xdat,
               degree=degree,
               Bernstein=Bernstein)
  
    console <- newLineConsole()
    console <- printClear(console)
    console <- printPop(console)
  
    if(all(bw.gamma < bw.switch)) {
      ## No bandwidths hit their upper bounds
      lscv <- minimand.cv.ls(bws=bw.gamma,
                             ydat=ydat,
                             xdat=xdat,
                             degree=degree,
                             W=W,
                             ckertype=ckertype,
                             ckerorder=ckerorder,
                             ukertype=ukertype,
                             okertype=okertype,
                             bwtype=bwtype,
                             cv.shrink=cv.shrink,
                             cv.maxPenalty=cv.maxPenalty,
                             cv.warning=cv.warning,
                             bandwidth.scale.categorical=bandwidth.scale.categorical,
                             ...)
    } else if(all(bw.gamma >= bw.switch)) {
      ## All bandwidths hit their upper bounds
      if(all(degree==0)) {
        lscv <- mean((ydat-mean(ydat))^2)
      } else {
        model <- lm(ydat~W-1)
        htt <- hatvalues(model)
        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
        lscv <- mean((residuals(model)/(1-htt))^2)
      }
    } else {
      ## Some bandwidths hit their upper bounds
      degree.sub <- degree[bw.gamma[xdat.numeric]<bw.switch[xdat.numeric]]
      lscv <- minimand.cv.ls(bws=bw.gamma[bw.gamma<bw.switch],
                             ydat=ydat,
                             xdat=xdat[,bw.gamma<bw.switch,drop=FALSE],
                             degree=if(length(degree.sub)==0){0}else{degree.sub},
                             W=W,
                             ckertype=ckertype,
                             ckerorder=ckerorder,
                             ukertype=ukertype,
                             okertype=okertype,
                             bwtype=bwtype,
                             cv.shrink=cv.shrink,
                             cv.maxPenalty=cv.maxPenalty,
                             cv.warning=cv.warning,
                             bandwidth.scale.categorical=bandwidth.scale.categorical,
                             ...)
    }
    console <- printPush("\r                                                                         ",console = console)
    console <- printPush(paste("\rfv = ",format(lscv)," ",sep=""),console = console)
    return(lscv)
  }
  
  ## Jan 19 2014, optional args error reported to zhenghua... want to
  ## pass ... but having ... in eval.aicc and in snomadr causes an error
  ## message thrown by checkFunctionArguments
  
  ## In snomadr() we write wrappers around user-defined functions to
  ## pass additional arguments
  ##  eval.f.wrapper <- function(x){ eval.f(x,...) }
  
  eval.aicc <- function(input,
                        params){
  
    ydat <- params$ydat
    xdat <- params$xdat
    xdat.numeric <- params$xdat.numeric
    num.bw <- params$num.bw
    num.numeric <- params$num.numeric
    cv.maxPenalty <- params$cv.maxPenalty
    degree <- params$degree
    cv <- params$cv
    ckertype <- params$ckertype
    ckerorder <- params$ckerorder
    ukertype <- params$ukertype
    okertype <- params$okertype
    bwtype <- params$bwtype
    cv.shrink <- params$cv.shrink
    cv.warning <- params$cv.warning
    Bernstein <- params$Bernstein
    bw.switch <- params$bw.switch
    bandwidth.scale.categorical=params$bandwidth.scale.categorical
  
    bw.gamma <- input[1:num.bw]
    if(cv=="degree-bandwidth")
      degree <- round(input[(num.bw+1):(num.bw+num.numeric)])
  
    ## Test for negative degrees of freedom
  
    if(dim.bs(basis="glp",kernel=TRUE,degree=degree,segments=rep(1,length(degree)))>length(ydat)-1)
      return(cv.maxPenalty)
  
    W <- W.glp(xdat=xdat,
               degree=degree,
               Bernstein=Bernstein)
  
    console <- newLineConsole()
    console <- printClear(console)
    console <- printPop(console)
  
    if(all(bw.gamma < bw.switch)) {
      ## No bandwidths hit their upper bounds
      aicc <- minimand.cv.aic(bws=bw.gamma,
                              ydat=ydat,
                              xdat=xdat,
                              degree=degree,
                              W=W,
                              ckertype=ckertype,
                              ckerorder=ckerorder,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,
                              cv.shrink=cv.shrink,
                              cv.maxPenalty=cv.maxPenalty,
                              cv.warning=cv.warning,
                              bandwidth.scale.categorical=bandwidth.scale.categorical,
                              ...)
    } else if(all(bw.gamma >= bw.switch)) {
      ## All bandwidths hit their upper bounds
      if(all(degree==0)) {
        aicc <- mean((ydat-mean(ydat))^2)
      } else {
        model <- lm(ydat~W-1)
        htt <- hatvalues(model)
        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
        aicc <- mean((residuals(model)/(1-htt))^2)
      }
    } else {
      ## Some bandwidths hit their upper bounds
      degree.sub <- degree[bw.gamma[xdat.numeric]<bw.switch[xdat.numeric]]
      aicc <- minimand.cv.aic(bws=bw.gamma[bw.gamma<bw.switch],
                              ydat=ydat,
                              xdat=xdat[,bw.gamma<bw.switch,drop=FALSE],
                              degree=if(length(degree.sub)==0){0}else{degree.sub},
                              W=W,
                              ckertype=ckertype,
                              ckerorder=ckerorder,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,
                              cv.shrink=cv.shrink,
                              cv.maxPenalty=cv.maxPenalty,
                              cv.warning=cv.warning,
                              bandwidth.scale.categorical=bandwidth.scale.categorical,
                              ...)
    }
    console <- printPush("\r                                                                         ",console = console)
    console <- printPush(paste("\rfv = ",format(aicc)," ",sep=""),console = console)
    return(aicc)
  }

  ## Save the seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  ckertype <- match.arg(ckertype)
  ukertype <- match.arg(ukertype)
  okertype <- match.arg(okertype)
  bwtype <- match.arg(bwtype)
  if(!any(ckerorder==c(2,4,6,8))) stop("ckerorder must be 2, 4, 6, or 8")

  bwmethod <- match.arg(bwmethod)
  cv <- match.arg(cv)

  if(is.null(ydat)) stop(" Error: You must provide y data")
  if(is.null(xdat)) stop(" Error: You must provide X data")
  if(!is.null(nmulti) && nmulti < 1) stop(paste(" Error: nmulti must be a positive integer (minimum 1)\nnmulti is (", nmulti, ")\n",sep=""))

  if(!is.null(degree) && any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

  if(degree.min < 0 ) stop(" Error: degree.min must be a non-negative integer")
  if(degree.max < degree.min) stop(" Error: degree.max must exceed degree.min")

  if(bandwidth.min < 0) stop(" Error: bandwidth.min must be non-negative")
  if(bandwidth.max < bandwidth.min) stop(" Error: bandwidth.max must exceed bandwidth.min")

  if(cv=="degree-bandwidth") {
    if(!is.null(degree) && any(degree>degree.max)) stop(" Error: degree supplied but exceeds degree.max")
    if(!is.null(degree) && any(degree<degree.min)) stop(" Error: degree supplied but less than degree.min")
  }

  ## For nearest neighbour bandwidths override default bandwidth.min
  ## and bandwidth.max and use sample size information.

  num.bw <- ncol(xdat)
  num.obs <- nrow(xdat)

  if(bwtype!="fixed") {
    bandwidth.min <- 1
    bandwidth.max <- num.obs-1
  }

  xdat <- as.data.frame(xdat)

  if(!is.null(bandwidth) && (length(bandwidth) != num.bw)) stop(" Error: bandwidth supplied but length not compatible with X data")

  if(is.null(nmulti)) nmulti <- min(5,num.bw)

  ## Determine which predictors are categorical and which are
  ## numeric... we care about unordered categorical kernels if the
  ## Aitchison & Aitken kernel is used since its bandwidth bounds are
  ## [0,(c-1)/c] and not 0/1 as are the rest of the unordered and
  ## ordered kernel bandwidth bounds.

  xdat.numeric <- sapply(1:num.bw,function(i){is.numeric(xdat[,i])})
  num.numeric <- ncol(as.data.frame(xdat[,xdat.numeric]))
  numeric.index <- which(xdat.numeric==TRUE)

  if(num.numeric == 0) stop("generalized local polynomial regression requires at least one numeric predictor")

  if(!is.null(degree) && length(degree) != num.numeric) stop(paste(" Error: degree vector supplied has ", length(degree), " elements but there exist ", num.numeric," numeric.predictors",sep=""))

  xdat.unordered <- sapply(1:num.bw,function(i){is.factor(xdat[,i])&!is.ordered(xdat[,i])})
  num.unordered <- ncol(as.data.frame(xdat[,xdat.unordered]))

  if(cv=="degree-bandwidth") {
    bbin <- c(rep(0, num.bw), rep(1, num.numeric))
    lb <- c(rep(bandwidth.min, num.bw), rep(degree.min, num.numeric))
    ub <- c(rep(bandwidth.max, num.bw), rep(degree.max, num.numeric))
  } else {
    bbin <- c(rep(0, num.bw))
    lb <- c(rep(bandwidth.min, num.bw))
    ub <- c(rep(bandwidth.max, num.bw))
  }

  ## This is to provide flexibility wrt ridging... very small
  ## bandwidths will wreak havoc on local polynomial regression. But
  ## sensible lower bounds exist. They could be data-driven (minimum
  ## unique distance between observations, say) or not. Regardless,
  ## using the default for both categorical and numeric bandwidths is
  ## likely not sensible.

  if(!is.null(bandwidth.min.numeric)) {
    for(i in 1:num.numeric) {
      lb[numeric.index[i]] <- bandwidth.min.numeric
    }
  }
  ## The input `bandwidth.switch' is the number of (scaled) standard
  ## deviations that will trigger the move to the global categorical
  ## kernel weighted polynomial fit for numeric predictors (now that
  ## we are using absolute initial meshes the upper bound for search
  ## is .Machine$double.xmax and the algorithm will step quickly into
  ## this region so this ought to assure k(0) holds). ub will be used
  ## for the categorical predictors. This helps speed up kd-guided
  ## search in the presence of large bandwidths (which otherwise would
  ## descend to the root of the tree).

  bw.switch <- c(rep(bandwidth.switch, num.bw))

  if(bwtype!="fixed" && num.numeric > 0) {
    for(i in 1:num.numeric) {
      k.max <- knn.max(xdat[,numeric.index[i]])
      if(ub[numeric.index[i]] > k.max) ub[numeric.index[i]] <- k.max
      if(bw.switch[numeric.index[i]] > k.max) bw.switch[numeric.index[i]] <- k.max
    }
  }

  ## Set parameters for NOMAD per variables, and the upper bounds,
  ## `trigger' for switching etc.

  INITIAL.MESH.SIZE <- list()
  MIN.MESH.SIZE <- list()
  MIN.POLL.SIZE <- list()

  for(i in 1:num.bw) {
    INITIAL.MESH.SIZE[[i]] <- initial.mesh.size.real
    MIN.MESH.SIZE[[i]] <- min.mesh.size.real
    MIN.POLL.SIZE[[i]] <- min.poll.size.real
    ## Need to do integer search for numeric predictors when bwtype is
    ## a nearest-neighbour, so set bbin appropriately.
    if(xdat.numeric[i]==TRUE && bwtype!="fixed") {
      bbin[i] <- 1
    }
    if(!xdat.numeric[i]) {
      lb[i] <- lb[i]*bandwidth.scale.categorical
      ub[i] <- 1*bandwidth.scale.categorical
      bw.switch[i] <- ub[i]
      INITIAL.MESH.SIZE[[i]] <- initial.mesh.size.integer
      MIN.MESH.SIZE[[i]] <- min.mesh.size.integer
      MIN.POLL.SIZE[[i]] <- min.poll.size.integer
    }
    ## Check for unordered and Aitchison/Aitken kernel
    if(xdat.unordered[i]==TRUE && ukertype=="aitchisonaitken") {
      c.num <- length(unique(xdat[,i]))
      ub[i] <- (c.num-1)/c.num*bandwidth.scale.categorical
      bw.switch[i] <- ub[i]
    }
  }

  ## Use degree for initial values if provided, and check whether the
  ## polynomial basis is well-conditioned or not. If it is, use the
  ## maximum well-conditioned basis.

  ill.conditioned <- check.max.degree(xdat,rep(degree.max,num.numeric),issue.warning=TRUE,Bernstein=Bernstein)
  degree.max.vec <- attr(ill.conditioned, "degree.max.vec")

  if(cv == "degree-bandwidth") {
    ub[(num.bw+1):(num.bw+num.numeric)] <- ifelse(ub[(num.bw+1):(num.bw+num.numeric)] > degree.max.vec, degree.max.vec, ub[(num.bw+1):(num.bw+num.numeric)])
    for(i in (num.bw+1):(num.bw+num.numeric)) {
      INITIAL.MESH.SIZE[[i]] <- initial.mesh.size.integer
      MIN.MESH.SIZE[[i]] <- min.mesh.size.integer
      MIN.POLL.SIZE[[i]] <- min.poll.size.integer
    }
  }

  ## Assign the NOMAD parameters to opts which is passed to snomadr()

  opts$"EPSILON" <- min.epsilon
  opts$"MAX_BB_EVAL" <- max.bb.eval
  opts$"INITIAL_MESH_SIZE" <- INITIAL.MESH.SIZE
  opts$"MIN_MESH_SIZE" <- MIN.MESH.SIZE
  opts$"MIN_POLL_SIZE" <- MIN.POLL.SIZE

  if(is.null(degree)) {
    if(cv == "degree-bandwidth") {
      degree <- numeric(num.numeric)
      for(i in 1:num.numeric) {
        ## We adopt the convention to always start from a polynomial
        ## or degree 1 for the first search attempt
        degree[i] <- 1
      }
    }
    else {
      stop(paste(" Error: degree must be given when optimizing only bandwidth"))
    }
  }

  ## Here we attempt to do some `smart' branching during search. If
  ## the bandwidth for a categorical predictor hits its upper bound it
  ## is `smoothed out', its kernel becomes a constant function, and
  ## the categorical predictor and does not influence the fit,
  ## delete-one or otherwise. If the bandwidth for any continuous
  ## predictor is `large' (> bandwidth.max), the variable is smoothed
  ## out if its polynomial degree is zero otherwise it is the global
  ## polynomial fit (W) in that dimension If all bandwidths hit their
  ## upper bound/are large, we get the global polynomial OLS fit. Note
  ## we have both the ls.cv and aic.cv methods.

  ## Generate the params fed to the snomadr solver

  params <- list()
  params$xdat <- xdat
  params$ydat <- ydat
  params$xdat.numeric <- xdat.numeric
  params$num.bw <- num.bw
  params$num.numeric <- num.numeric
  params$cv.maxPenalty <- cv.maxPenalty
  params$cv <- cv
  params$degree <- degree
  params$ckertype <- ckertype
  params$ckerorder <- ckerorder
  params$ukertype <- ukertype
  params$okertype <- okertype
  params$bwtype <- bwtype
  params$cv.shrink <- cv.shrink
  params$cv.warning <- cv.warning
  params$Bernstein <- Bernstein
  params$bw.switch <- bw.switch
  params$bandwidth.scale.categorical=bandwidth.scale.categorical

  ## Multistarting

  fv.vec <- numeric(nmulti)

  ## Whether or not to display the information in snomadr

  print.output <- FALSE
  console <- newLineConsole()
  if(!is.null(opts$DISPLAY_DEGREE)){
    if(opts$DISPLAY_DEGREE>0){
      print.output <-TRUE
      console <- printPush("\rCalling NOMAD (Nonsmooth Optimization by Mesh Adaptive Direct Search)\n",console = console)
    }
  } else {
    print.output <-TRUE
    console <- printPush("\rCalling NOMAD (Nonsmooth Optimization by Mesh Adaptive Direct Search)\n",console = console)
  }

  ## Use bandwidth for initial values if provided

  if(is.null(bandwidth)) {
    init.search.vals <- numeric()
    for(i in 1:num.bw) {
      if(xdat.numeric[i]==TRUE && bwtype=="fixed") {
        init.search.vals[i] <- runif(1,lb[i],1.5)
      }
      if(xdat.numeric[i]==TRUE && bwtype!="fixed") {
        init.search.vals[i] <- round(runif(1,lb[i],sqrt(ub[i])))
      }
      if(xdat.numeric[i]!=TRUE) {
        init.search.vals[i] <- runif(1,lb[i],ub[i])
      }
    }
  } else {
    ## If bandwidths are provided, need to convert those for the
    ## numeric variables into scale factors (level on which snomadr is
    ## optimizing)
    init.search.vals <- bandwidth
    for(i in 1:num.bw) {
      if(xdat.numeric[i]==TRUE && bwtype=="fixed") {
        init.search.vals[i] <- bandwidth[i]/(scale.robust(xdat[,i])*length(ydat)^{-1/(num.numeric+2*ckerorder)})
      }
      if(xdat.numeric[i]!=TRUE) {
        init.search.vals[i] <- bandwidth[i]*bandwidth.scale.categorical
      }
    }
  }

  ## Generate all initial points for the multiple restarting

  x0.pts <- matrix(numeric(1), nmulti, length(bbin))

  for(iMulti in 1:nmulti) {
    if(iMulti != 1) {
      init.search.vals <- numeric()
      for(i in 1:num.bw) {
        if(xdat.numeric[i]==TRUE && bwtype=="fixed") {
          init.search.vals[i] <- runif(1,lb[i],1.5)
        }
        if(xdat.numeric[i]==TRUE && bwtype!="fixed") {
          init.search.vals[i] <- round(runif(1,lb[i],sqrt(ub[i])))
        }
        if(xdat.numeric[i]!=TRUE) {
          init.search.vals[i] <- runif(1,lb[i],ub[i])
        }
      }
    }

    if(cv == "degree-bandwidth" && iMulti != 1) {
      for(i in 1:num.numeric) {
        degree[i] <- sample(degree.min:ub[(num.bw+1):(num.bw+num.numeric)][i], 1)
      }
    }

    if(cv =="degree-bandwidth"){
      x0.pts[iMulti, ] <- c(init.search.vals, degree)
    }
    else {
      x0.pts[iMulti, ] <- c(init.search.vals)
    }

  }

  if(bwmethod == "cv.ls" ) {
      solution<-snomadr(eval.f=eval.lscv,
                        n=length(bbin),
                        x0=as.numeric(x0.pts),
                        bbin=bbin,
                        bbout=0,
                        lb=lb,
                        ub=ub,
                        nmulti=ifelse(nmulti==1,0,nmulti),
                        random.seed=random.seed,
                        opts=opts,
                        print.output=print.output,
                        params=params,
                        ...);

      if(restart.from.min) solution<-snomadr(eval.f=eval.lscv,
                                             n=length(bbin),
                                             x0=solution$solution,
                                             bbin=bbin,
                                             bbout=0,
                                             lb=lb,
                                             ub=ub,
                                             nmulti=1,
                                             random.seed=random.seed,
                                             opts=opts,
                                             print.output=print.output,
                                             params=params,
                                             ...);

  } else {
      solution<-snomadr(eval.f=eval.aicc,
                        n=length(bbin),
                        x0=as.numeric(x0.pts),
                        bbin=bbin,
                        bbout=0,
                        lb=lb,
                        ub=ub,
                        nmulti=ifelse(nmulti==1,0,nmulti),
                        random.seed=random.seed,
                        opts=opts,
                        print.output=print.output,
                        params=params,
                        ...);

      if(restart.from.min) solution<-snomadr(eval.f=eval.aicc,
                                             n=length(bbin),
                                             x0=solution$solution,
                                             bbin=bbin,
                                             bbout=0,
                                             lb=lb,
                                             ub=ub,
                                             nmulti=1,
                                             random.seed=random.seed,
                                             opts=opts,
                                             print.output=print.output,
                                             params=params,
                                             ...);
  }

  fv.vec[1] <- solution$objective

  bw.opt.sf <- solution$solution[1:num.bw]

  ## We optimize at the level of scaling factors (multiples of
  ## bandwidths) so we need to back out the unscaled (raw) bandwidths.

  bw.opt <- bw.opt.sf

  if(bwtype=="fixed") {
    for(i in 1:num.numeric) {
      sd.xdat <- scale.robust(xdat[,numeric.index[i]])
      bw.opt[numeric.index[i]] <- bw.opt[numeric.index[i]]*sd.xdat*length(ydat)^{-1/(num.numeric+2*ckerorder)}
    }
  }

  for(i in 1:num.bw) {
    if(!xdat.numeric[i]) {
      bw.opt[i] <- bw.opt[i]/bandwidth.scale.categorical
    }
  }

  if(cv == "degree-bandwidth") {
    degree.opt <- solution$solution[(num.bw+1):(num.bw+num.numeric)]
    for(i in num.numeric){
      if(isTRUE(all.equal(degree.opt[i],ub[num.bw+i]))) warning(paste(" Optimal degree for numeric predictor ",i," equals search upper bound (", ub[num.bw+i],")",sep=""))
    }
  } else {
    degree.opt <- degree
  }

  if(!is.null(opts$MAX_BB_EVAL)){
    if(nmulti>0) {if(nmulti*opts$MAX_BB_EVAL <= solution$bbe) warning(paste(" MAX_BB_EVAL reached in NOMAD: perhaps use a larger value...", sep=""))}
    if(nmulti==0) {if(opts$MAX_BB_EVAL <= solution$bbe) warning(paste(" MAX_BB_EVAL reached in NOMAD: perhaps use a larger value...", sep="")) }
  }

  fv <- solution$objective

  best <- NULL
  numimp <- 0

  for(i in num.bw) {
    if(isTRUE(all.equal(bw.opt[i],lb[i]))) warning(paste(" Optimal bandwidth for predictor ",i," equals search lower bound (", formatC(lb[i],digits=3,format="g"),"): rerun with smaller bandwidth.min",sep=""))
  }

  console <- printPush("\r                        ",console = console)
  console <- printClear(console)
  console <- printPop(console)

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  if(isTRUE(all.equal(fv,cv.maxPenalty))) stop(" Search failed: restart with larger nmulti or smaller degree.max")

  return(list(bws=bw.opt,
              bws.sf=bw.opt.sf,
              fv=fv,
              numimp=numimp,
              best=best,
              fv.vec=fv.vec,
              degree=degree.opt,
              bwtype=bwtype,
              ckertype=ckertype,
              ckerorder=ckerorder,
              ukertype=ukertype,
              okertype=okertype))

}

compute.bootstrap.errors <- function(tydat,
                                     txdat,
                                     exdat,
                                     eydat,
                                     bws,
                                     degree,
                                     ckertype,
                                     ckerorder,
                                     ukertype,
                                     okertype,
                                     bwtype,boot.object=c("fitted","gradient","gradient.categorical"),
                                     plot.errors.boot.num=99,
                                     plot.errors.type=c("quantiles","standard"),
                                     plot.errors.quantiles=c(.025,.975),
                                     alpha=0.05,
                                     gradient.vec=NULL,
                                     gradient.categorical=FALSE,
                                     gradient.categorical.index=NULL,
                                     Bernstein=TRUE,
                                     ...){

  plot.errors.type <- match.arg(plot.errors.type)
  boot.object <- match.arg(boot.object)

  if(missing(exdat)) {
    neval <- nrow(txdat)
    exdat <- txdat
  } else {
    neval <- nrow(exdat)
  }

  boot.err <- matrix(data = NA, nrow = neval, ncol = 2)

  ## Conduct simple iid residual bootstrap

  boot.func.mean <- function(model.fitted,indices){
    est.boot <- npglpreg(tydat=model.fitted+(tydat-model.fitted)[indices],
                         txdat=txdat,
                         exdat=exdat,
                         bws=bws,
                         degree=degree,
                         ckertype=ckertype,
                         ckerorder=ckerorder,
                         ukertype=ukertype,
                         okertype=okertype,
                         bwtype=bwtype,
                         gradient.vec=gradient.vec,
                         gradient.categorical=gradient.categorical,
                         gradient.categorical.index=gradient.categorical.index,
                         Bernstein=Bernstein,
                         ...)
    if(boot.object=="fitted") {
      return(est.boot$fitted.values)
    } else if(boot.object=="gradient") {
      return(est.boot$gradient)
    } else if(boot.object=="gradient.categorical") {
      return(est.boot$gradient.categorical.mat[,gradient.categorical.index])
    }
  }

  ## Fitted values for the training data required

  est <- npglpreg(tydat=tydat,
                  txdat=txdat,
                  bws=bws,
                  degree=degree,
                  ckertype=ckertype,
                  ckerorder=ckerorder,
                  ukertype=ukertype,
                  okertype=okertype,
                  bwtype=bwtype,
                  Bernstein=Bernstein,
                  ...)

  model.fitted <- est$fitted.values

  boot.out <- boot(data = model.fitted,
                   statistic = boot.func.mean,
                   R = plot.errors.boot.num)

  if (plot.errors.type == "standard") {
    boot.err[,1:2] <- qnorm(1-alpha/2)*sqrt(diag(cov(boot.out$t)))
    boot.err[,1] <- boot.out$t0 - boot.err[,1]
    boot.err[,2] <- boot.out$t0 + boot.err[,2]
  }
  else if (plot.errors.type == "quantiles") {
    boot.err[,1:2] <- t(sapply(as.data.frame(boot.out$t),
                               function (y) {
                                 quantile(y,probs = plot.errors.quantiles)
                               }))
  }

  return(cbind(boot.out$t0,boot.err))

}

plot.npglpreg <- function(x,
                          mean=TRUE,
                          deriv=0,
                          ci=FALSE,
                          num.eval=100,
                          common.scale=TRUE,
                          xtrim = 0.0,
                          xq = 0.5,
                          plot.behavior = c("plot","plot-data","data"),
                          plot.errors.boot.num=99,
                          plot.errors.type=c("quantiles","standard"),
                          plot.errors.quantiles=c(.025,.975),
                          persp.rgl=FALSE,
                          ...) {

  plot.behavior <- match.arg(plot.behavior)
  plot.errors.type <- match.arg(plot.errors.type)

  object <- x

  ## Needed for correctly obtaining predictions

  degree <- object$degree
  bws <- object$bws
  bwtype <- object$bwtype
  ckertype <- object$ckertype
  ckerorder <- object$ckerorder
  ukertype <- object$ukertype
  okertype <- object$okertype
  Bernstein <- object$Bernstein

  txdat <- object$x
  tydat <- object$y

  xq <- double(ncol(txdat)) + xq

  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)
  console <- printPush("\rWorking...",console = console)

  ## Mean

  if(mean==TRUE && deriv==0) {

    if(!persp.rgl) {

      mg <- list()

      for(i in 1:NCOL(object$x)) {

        if(!is.factor(object$x[,i])) {
          exdat <- matrix(NA,nrow=num.eval,ncol=NCOL(object$x))
          neval <- num.eval
        } else {
          neval <- length(unique(object$x[,i]))
          exdat <- matrix(NA,nrow=neval,ncol=NCOL(object$x))
        }

        exdat <- data.frame(exdat)

        if(!is.factor(object$x[,i])) {
          xlim <- trim.quantiles(object$x[,i],xtrim)
          exdat[,i] <- seq(xlim[1],xlim[2],length=neval)
        } else {
          exdat[,i] <- sort(unique(object$x[,i]))
        }

        for(j in (1:NCOL(object$x))[-i]) {
          exdat[,j] <- rep(uocquantile(object$x[,j],prob=xq[j]),neval)
        }

        names(exdat) <- object$xnames

        est <- npglpreg.default(tydat=tydat,
                                txdat=txdat,
                                exdat=exdat,
                                bws=bws,
                                degree=degree,
                                ckertype=ckertype,
                                ckerorder=ckerorder,
                                ukertype=ukertype,
                                okertype=okertype,
                                bwtype=bwtype,
                                Bernstein=Bernstein,
                                ...)

        fitted.values <- est$fitted.values

        if(!ci) {

          mg[[i]] <- data.frame(exdat[,i],fitted.values)
          names(mg[[i]]) <- c(names(exdat)[i],"mean")

        } else {

          console <- printClear(console)
          console <- printPop(console)
          console <- printPush(paste("\rConducting ",plot.errors.boot.num," bootstrap resamples for predictor ",i,"...",sep=""),console = console)

          ci.out <- compute.bootstrap.errors(tydat=tydat,
                                             txdat=txdat,
                                             exdat=exdat,
                                             bws=bws,
                                             degree=degree,
                                             ckertype=ckertype,
                                             ckerorder=ckerorder,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,
                                             boot.object="fitted",
                                             plot.errors.boot.num=plot.errors.boot.num,
                                             plot.errors.type=plot.errors.type,
                                             plot.errors.quantiles=plot.errors.quantiles)

          mg[[i]] <- data.frame(exdat[,i],ci.out)
          names(mg[[i]]) <- c(names(exdat)[i],"mean","lwr","upr")

        }

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

        if(!is.null(object$num.categorical)||(object$num.numeric>1)) par(mfrow=dim.plot(NCOL(object$x)))

        for(i in 1:NCOL(object$x)) {

          if(!ci) {

            plot(mg[[i]][,1],mg[[i]][,2],
                 xlab=names(exdat)[i],
                 ylab="Conditional Mean",
                 ylim=ylim,
                 type="l")

          } else {
            if(!common.scale) ylim <- c(min(mg[[i]][,-1]),max(mg[[i]][,-1]))
            plot(mg[[i]][,1],mg[[i]][,2],
                 xlab=names(exdat)[i],
                 ylab="Conditional Mean",
                 ylim=ylim,
                 type="l")
            ## Need to overlay for proper plotting of factor errorbars
            par(new=TRUE)
            plot(mg[[i]][,1],mg[[i]][,3],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2)
            par(new=TRUE)
            plot(mg[[i]][,1],mg[[i]][,4],
                 xlab="",
                 ylab="",
                 ylim=ylim,
                 type="l",
                 axes=FALSE,
                 lty=2,
                 col=2)
          }

        }

      }

    } else {

#      if(!require(rgl)) stop(" Error: you must first install the rgl package")

      if(object$num.categorical != 0) stop(" Error: persp3d is for continuous predictors only")
      if(object$num.numeric != 2) stop(" Error: persp3d is for cases involving two continuous predictors only")

      newdata <- matrix(NA,nrow=num.eval,ncol=2)
      newdata <- data.frame(newdata)

      xlim <- trim.quantiles(object$x[,1],xtrim)
      ylim <- trim.quantiles(object$x[,2],xtrim)

      x1.seq <- seq(xlim[1],xlim[2],length=num.eval)
      x2.seq <- seq(ylim[1],ylim[2],length=num.eval)
      x.grid <- expand.grid(x1.seq,x2.seq)
      newdata <- data.frame(x.grid[,1],x.grid[,2])
      names(newdata) <- names(object$x)

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
                xlab=names(object$x)[1],ylab=names(object$x)[2],zlab="Y",
                ticktype="detailed",
                border="red",
                color=col,
                alpha=.7,
                back="lines",
                main="Conditional Mean")

        grid3d(c("x", "y+", "z"))

        play3d(spin3d(axis=c(0,0,1), rpm=5), duration=15)

      }

    }

    if(plot.behavior!="plot") {
      console <- printClear(console)
      console <- printPop(console)
      return(mg)
    }

  }

  ## deriv

  if(deriv > 0) {

    rg <- list()

    i.numeric <- 1
    i.categorical <- 1

    for(i in 1:NCOL(object$x)) {

      gradient.vec <- NULL

      if(!is.factor(object$x[,i])) {
        newdata <- matrix(NA,nrow=num.eval,ncol=NCOL(object$x))
        neval <- num.eval
        gradient.vec <- rep(0,object$num.numeric)
        gradient.vec[i.numeric] <- deriv
      } else {
        neval <- length(unique(object$x[,i]))
        newdata <- matrix(NA,nrow=neval,ncol=NCOL(object$x))
      }

      newdata <- data.frame(newdata)

      if(!is.factor(object$x[,i])) {
        xlim <- trim.quantiles(object$x[,i],xtrim)
        newdata[,i] <- seq(xlim[1],xlim[2],length=neval)
      } else {
        newdata[,i] <- sort(unique(object$x[,i]))
      }

      for(j in (1:NCOL(object$x))[-i]) {
        newdata[,j] <- rep(uocquantile(object$x[,j],prob=xq[j]),neval)
      }

      names(newdata) <- object$xnames

      est <- npglpreg.default(tydat=tydat,
                              txdat=txdat,
                              exdat=newdata,
                              bws=bws,
                              degree=degree,
                              ckertype=ckertype,
                              ckerorder=ckerorder,
                              ukertype=ukertype,
                              okertype=okertype,
                              bwtype=bwtype,
                              gradient.vec=gradient.vec,
                              gradient.categorical=TRUE,
                              Bernstein=Bernstein,
                              ...)

      if(!is.factor(object$x[,i])) {
        fitted.values <- est$gradient
      } else {
        fitted.values <- est$gradient.categorical.mat[,i.categorical]
      }

      if(!ci) {

        rg[[i]] <- data.frame(newdata[,i],fitted.values)
        names(rg[[i]]) <- c(names(newdata)[i],"deriv")

      } else {

        console <- printClear(console)
        console <- printPop(console)
        console <- printPush(paste("\rConducting ",plot.errors.boot.num," bootstrap resamples for predictor ",i,"...",sep=""),console = console)

        if(!is.factor(object$x[,i])) {
          ci.out <- compute.bootstrap.errors(tydat=tydat,
                                             txdat=txdat,
                                             exdat=newdata,
                                             bws=bws,
                                             degree=degree,
                                             ckertype=ckertype,
                                             ckerorder=ckerorder,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,
                                             boot.object="gradient",
                                             plot.errors.boot.num=plot.errors.boot.num,
                                             plot.errors.type=plot.errors.type,
                                             plot.errors.quantiles=plot.errors.quantiles,
                                             gradient.vec=gradient.vec)
        } else {
          ci.out <- compute.bootstrap.errors(tydat=tydat,
                                             txdat=txdat,
                                             exdat=newdata,
                                             bws=bws,
                                             degree=degree,
                                             ckertype=ckertype,
                                             ckerorder=ckerorder,
                                             ukertype=ukertype,
                                             okertype=okertype,
                                             bwtype=bwtype,
                                             boot.object="gradient.categorical",
                                             plot.errors.boot.num=plot.errors.boot.num,
                                             plot.errors.type=plot.errors.type,
                                             plot.errors.quantiles=plot.errors.quantiles,
                                             gradient.categorical=TRUE,
                                             gradient.categorical.index=i.categorical)
        }

        rg[[i]] <- data.frame(newdata[,i],ci.out)
        names(rg[[i]]) <- c(names(newdata)[i],"deriv","lwr","upr")

      }

      if(!is.factor(object$x[,i])) {
        i.numeric <- i.numeric + 1
      } else {
        i.categorical <- i.categorical + 1
      }

    }

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

      if(!is.null(object$num.categorical)||(object$num.numeric>1)) par(mfrow=dim.plot(NCOL(object$x)))

      for(i in 1:NCOL(object$x)) {

          if(!ci) {
          plot(rg[[i]][,1],rg[[i]][,2],
               xlab=names(newdata)[i],
               ylab=ifelse(!is.factor(newdata[,i]), paste("Order", deriv,"Gradient"), "Difference in Levels"),
               ylim=ylim,
               type="l")

        } else {
          if(!common.scale) ylim <- c(min(rg[[i]][,-1]),max(rg[[i]][,-1]))
          plot(rg[[i]][,1],rg[[i]][,2],
               xlab=names(newdata)[i],
               ylab=ifelse(!is.factor(newdata[,i]), paste("Order", deriv,"Gradient"), "Difference in Levels"),
               ylim=ylim,
               type="l")
          ## Need to overlay for proper plotting of factor errorbars
          par(new=TRUE)
          plot(rg[[i]][,1],rg[[i]][,3],
               xlab="",
               ylab="",
               ylim=ylim,
               type="l",
               axes=FALSE,
               lty=2,
               col=2)
          par(new=TRUE)
          plot(rg[[i]][,1],rg[[i]][,4],
               xlab="",
               ylab="",
               ylim=ylim,
               type="l",
               axes=FALSE,
               lty=2,
               col=2)
        }

      }

    }

    if(plot.behavior!="plot") {
      console <- printClear(console)
      console <- printPop(console)
      return(rg)
    }

  }

  console <- printClear(console)
  console <- printPop(console)

  ## Reset par to 1,1 (can be modified above)

  if(!persp.rgl) par(mfrow=c(1,1))

}
