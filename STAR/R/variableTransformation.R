##############################################################################
##############################################################################
##############################################################################
mkM2U <- function(df,
                  vN,
                  low,
                  high,
                  delta,
                  alpha=2,
                  ...) {
#######################################################################
### Function mkM2U
### Finds the transformation making the distribution of one variable of 
### a data frame uniform on its definition domain. The transformatin
### is a smooth version of the ecdf computed optionally on a subset of
### of the data frame.
### ----------------------------------------------------------
### Arguments:
###  df: a data frame.
###  vN: a character or an integer specifying the variable to transform.
###  low: a numeric specifying the time from which the ecdf should be
###       be evaluated. If missing it is set to the smallest time.
###  high: a numeric specifying the time up to which the ecdf should be
###        be evaluated. If missing it is set to the largest time.
###  delta: the bin width used to build the histogram of the variable
###         values which is later smoothed.
###  alpha: alpha argument of ssden see this gss function documentation
###  ...: additional arguments passed to ssden
### -------------------------------------------------------------------
### Value: a smooth version of the ecdf.
### The returned function has in addition 4 attributes:
### qFct: the smooth quantile function, that is if mySmooth is what is
###       returned by mkM2U, then attr(mySmooth,"qFct")(mySmooth(x))=x
### dFct: the smooth pdf of the selected variable
### range: a two elements vector with the range of the variable
### call: the matched call
#######################################################################
  
  ii <- df[[vN]]
  if (missing(low)) low <- min(df$time)
  if (missing(high)) high <- max(df$time)  
  if (missing(delta)) {
    iii <- sort(diff(sort(unique(ii))))
    delta <- min(iii[iii>=diff(range(iii))/1000])
    rm(iii)
  } ## End of conditional on missing(delta)
  
  iiD <- c(min(df[[vN]])-.Machine$double.eps,
           max(df[[vN]])+.Machine$double.eps)
  iiRef <- ii[low <= df$time & df$time <= high]
  rm(df,ii,low,high)

  iiB <- seq(iiD[1],iiD[2]+delta,delta)
  iiH <- as.data.frame(hist(iiRef,breaks=iiB,plot=FALSE)[c("mids","counts")])
  names(iiH) <- c("x","counts")
  riiD <- range(iiH$x)+delta*c(-1,1)/2
  ii.fit <- ssden(~x, data=iiH,
                  domain=data.frame(x=riiD),
                  weights=iiH$counts,alpha=alpha,
                  ...
                  )
  rm(iiH,iiB,iiRef,iiD)

  mn <- min(ii.fit$domain)
  mx <- max(ii.fit$domain)
  ## find out the true integral of the density function on its definition
  ## domain
  Z <- integrate(function(x) dssden(ii.fit,x),mn,mx,subdivisions=1000)$value
  ## define dFct returning a density actually summing to 1
  dFct <- function(x) {
    result <- numeric(length(x))
    good <- mn <= x & x <= mx
    if (any(good)) result[good] <- dssden(ii.fit,x[good])/Z
    result
  }

  ## create a lookup table of quantiles and corresponding probabilities
  Q <- seq(mn,mx,len=101)
  P <- numeric(101)
  for (i in 2:101) P[i] <- P[i-1] + integrate(dFct,Q[i-1],Q[i])$value
  
  result <- function(q,...) {
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q <= mn] <- 0
    p[q >= mx] <- 1
    kk <- (1:length(q))[q > mn & q < mx]
    p.dup <- 0
    for (i in kk) {
        if (q.dup[i]) {
            p[i] <- p.dup
            next
        }
        idx <- findInterval(q[i],Q)
        lwr <- Q[idx]
        if (q[i]-lwr < 1e-5) {
          p[i] <- P[idx]+(P[idx+1]-P[idx])/(Q[idx+1]-lwr)*(q[i]-lwr)
        } else {
          p[i] <- integrate(dFct,
                            lower = lwr, upper = q[i],
                            ...)$value + P[idx]
          p.dup <- p[i]
        } ## End of conditional on q[i]-lwr < 1e-5 
    } ## End of loop on i
    p[order.q]
  }

  qFct <- function(p,...) {
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    kk <- (1:length(p))[p > 0 & p < 1]
    for (i in kk) {
      if (p.dup[i]) {
        q[i] <- q[i-1]
      } else {
        idx <- findInterval(p[i],P)
        if (identical(p[i],P[idx])) {
          q[i] <- Q[idx]
        } else {
          q[i] <- uniroot(function(x) result(x) - p[i],Q[idx:(idx+1)],...)$root
        }
      } ## End of conditional on p.dup[i]
    } ## End of loop on i
    q
  }
  
  attr(result,"qFct") <- qFct 
  attr(result,"dFct") <- dFct
  attr(result,"range") <- riiD
  attr(result,"call") <- match.call()
  result
}
##############################################################################
##############################################################################
##############################################################################
mkAR <- function(df,
                 low,
                 high,
                 max.order,
                 selfName="lN.1",
                 ...
                 ) {
#######################################################################
### Function mkAR
### Generates a data frame with variables suitable for an AR like model.
### These variables are the former inter spike intervals. The variables
### are moreover transformed with mkM2U inorder to have a uniform
### distribution on their definition domain.
### ----------------------------------------------------------
### Arguments:
###  df: a data frame.
###  low: a numeric specifying the time from which the ecdf should be
###       be evaluated. If missing it is set to the smallest time.
###  high: a numeric specifying the time up to which the ecdf should be
###        be evaluated. If missing it is set to the largest time.  
###  max.order: the maximal order of the process.
###  selfName: the name of the variable of df containing the elapsed time
###            since the last spike of the neuron whose discharge is
###            modeled.
###  ... : additional arguments passed to mkM2U
### -------------------------------------------------------------------
### Value: a data frame.
### In addition to the variables of df the returned data frame contains
### a variable "est" with the transformed elapsed time since the last
### spike of the neuron and "i1t", "i2t",...,"i max.order t", the
### transformed previous inter spike intervals.
### The returned data frame has also four attributes:
###  fmla: a formula suitable for a first argument of, say, gssanova
###  tfL: the function returned by mkM2U transforming the elasped time
###       since the last spike of the neuron.
###  tfI: the function returned by mkM2U transforming the first former
###       inter spike interval.
###  call: the matched call.
#######################################################################  
  if (missing(low)) low <- min(df$time) - .Machine$double.eps
  if (missing(high)) high <- max(df$time) + .Machine$double.eps
  vNames <- paste("i",1:max.order,sep="")
  vNames2 <- paste("i",1:max.order,"t",sep="")
  for (i in 1:max.order) df[[vNames[i]]] <- isi(df,lag=i)
  df <- df[complete.cases(df),]
  m2uL <- mkM2U(df,selfName,low,high,...)
  m2uI <- mkM2U(df,"i1",low,high,...)
  df <- within(df, est <- m2uL(lN.1))
  for (i in 1:max.order) df[[vNames2[i]]] <- m2uI(df[[vNames[i]]])
  for (i in 1:max.order) df[[vNames[i]]] <- NULL
  fmla <- as.formula(paste("event ~ est+", paste(vNames2, collapse= "+")))
  attr(df,"fmla") <- fmla
  attr(df,"m2uL") <- m2uL
  attr(df,"m2uI") <- m2uI
  attr(df,"call") <- match.call()
  df
}
##############################################################################
##############################################################################
##############################################################################
changeScale <- function(obj,
                        xFct,
                        yFct
                        ) {
#######################################################################
### Function changeScale
### Designed to transform results of quickPredict obtained on interaction
### terms from the transformed scale (on which the variables are approxi-
### mately uniformly distributed) onto the "native", linear scale.
### ----------------------------------------------------------
### Arguments:
###  obj: a quickPredict object.
###  xFct: a function to be applied on the "xx" element of obj.
###  yFct: a function to be applied on the "yy" element of obj.
### -------------------------------------------------------------------
### Value: a quickPredict object.
#######################################################################  
  ## Check that obj inherits of "quickPredict"
  if (!inherits(obj,"quickPredict")) stop("obj should be a quickPredict object.")

  if (missing(xFct)) xFct <- function(x) x
  if (missing(yFct)) yFct <- function(y) y

  est.mean <- obj$est.mean
  if (!is.null(obj$est.sd)) {
    est.sd <- obj$est.sd
  } else {
    est.sd <- NULL
  }

  theCall <- match.call()
  theCall[["obj"]] <- obj$call
  xx <- xFct(obj$xx)
  n <- length(xx)
  equal2min <- sum(abs(xx-min(xx)) <= .Machine$double.eps)
  equal2max <- sum(abs(xx-max(xx)) <= .Machine$double.eps)
  goodX <- !logical(n)
  if (equal2min > 1) goodX[1:(equal2min-1)] <- FALSE
  if (equal2max > 1) goodX[(n-equal2max+2):n] <- FALSE
  est.mean <- est.mean[goodX,]
  if (!is.null(est.sd)) est.sd <- est.sd[goodX,]
  xx <- xx[goodX]
  goodX <- !duplicated(xx)
  xx <- xx[goodX]
  est.mean <- est.mean[goodX,]
  if (!is.null(est.sd)) est.sd <- est.sd[goodX,]
  
  yy <- yFct(obj$yy)
  equal2min <- sum(abs(yy-min(yy)) <= .Machine$double.eps)
  equal2max <- sum(abs(yy-max(yy)) <= .Machine$double.eps)
  goodY <- !logical(n)
  if (equal2min > 1) goodY[1:(equal2min-1)] <- FALSE
  if (equal2max > 1) goodY[(n-equal2max+2):n] <- FALSE
  est.mean <- est.mean[,goodY]
  if (!is.null(est.sd)) est.sd <- est.sd[,goodY]
  yy <- yy[goodY]
  goodY <- !duplicated(yy)
  yy <- yy[goodY]
  est.mean <- est.mean[,goodY]
  if (!is.null(est.sd)) est.sd <- est.sd[,goodY]

  result <- list(xx=xx,
                 yy=yy,
                 call=theCall,
                 include=obj$include,
                 est.mean=est.mean,
                 est.sd=est.sd)
  class(result) <- "quickPredict"
  result
  
}
