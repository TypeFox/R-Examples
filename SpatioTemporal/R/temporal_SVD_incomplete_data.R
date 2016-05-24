#######################################################################
## FUNCTIONS THAT COMPUTE SVDs FOR TEMPORAL TRENDS WITH MISSING DATA ##
#######################################################################
##Functions in this file:
## SVDmiss             - EX:ok
## SVDsmooth           - EX:ok
## SVDsmoothCV         - EX:with SVDsmooth
## calcSmoothTrends    - EX:ok
## plot.SVDcv          - EX:with SVDsmooth
## boxplot.SVDcv       - EX:with SVDsmooth
## print.SVDcv         - EX:with SVDsmooth
## updateTrend.STdata  - EX:ok
## updateTrend.STmodel - EX:with updateTrend.STdata
## updateTrend         - EX:with updateTrend.STdata
## updateSTdataTrend   - EX:with updateTrend.STdata
## internalUpdateTrendFnc    - internal function for updateTrend
## internalUpdateTrendSmooth - internal function for updateTrend
## internalCreateTrendFnc    - internal function that creates a trend function

##' Function that completes a data matrix using iterative svd as described in
##' Fuentes et. al. (2006). The function iterates between computing the svd for
##' the matrix and replacing the missing values by linear regression of the
##' columns onto the first \code{ncomp} svd components. As initial replacement
##' for the missing values regression on the column averages are used. The
##' function \emph{will fail} if entire rows and/or columns are missing from the
##' data matrix.
##'
##' @title Missing Data SVD
##' @param X Data matrix, with missing values marked by \code{NA}.
##' @param niter Maximum number of iterations to run before exiting,
##'   \code{Inf} will run until the \code{conv.reldiff} criteria is met.
##' @param ncomp Number of SVD components to use in the reconstruction (>0).
##' @param conv.reldiff Assume the iterative procedure has converged when
##'   the relative difference between two consecutive iterations is less than
##'   \code{conv.reldiff}.
##' @return A list with the following components:
##'   \item{Xfill}{The completed data matrix with missing values replaced by
##'                fitting the data to the \code{ncomp} most important svd
##'                components}
##'   \item{svd}{The result of svd on the completed data matrix, i.e.
##'              \code{svd(Xfill)}}
##'   \item{status}{A vector of status variables: \code{diff}, the
##'                 absolute difference between the two last iterations;
##'                 \code{rel.diff}, the relative difference;
##'                 \code{n.iter}, the number of iterations; and
##'                 \code{max.iter}, the requested maximum number of
##'                 iterations.}
##' 
##' @references
##' M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
##'  Analyze Space-Time Processes in Statistical methods for spatio-temporal
##'  systems (B. Finkenstädt, L. Held, V. Isham eds.) 77-150
##' 
##' @author Paul D. Sampson and Johan Lindström
##'
##' @example Rd_examples/Ex_SVDmiss.R
##' @family SVD for missing data
##' @family data matrix
##' @export
SVDmiss <- function(X, niter=25, ncomp=min(4,dim(X)[2]), conv.reldiff=0.001)
{
  ##ensure at least one iteration
  niter <- max(niter,1)
  ##and ensure sane number of components
  if( ncomp<1 ){
    stop("ncomp should be >0, is: ", ncomp)
  }
  
  ##First find missing values
  Ina <- is.na(X)
  if( all(!Ina) ){
    ##if X has no missing data this is simple
    svd0 <- svd(X)
    XF <- X
    i <- diff <- reldiff <- 0
  }else{
    ##X has missing data, use iterative method
    ##Iterative svd calculation with missing data.
    ##Initial first element of U matrix is average curve.
    u1 <- rowMeans(X, na.rm = TRUE)
    XM <- matrix(1, nrow(X), ncol(X))
    XM[Ina] <- 0
    XZ <- X
    # v1 is proportional to X'u1/(u1'u1), but with calculations
    # filling in zeros for missing values in the sums.
    XZ[Ina] <- 0.
    # Fill in missing values for initial complete SVD calculation.
    # Then iterate  using only complete data.
    v1 <- diag(t(XZ) %*% (XM * u1))/diag(t(XM * u1) %*% (XM * u1))
    XF <- X
    XF[Ina] <- (matrix(u1, ncol = 1) %*% matrix(v1, nrow = 1))[Ina]
    if( any(is.na(XF)) )
      stop("Unable to complete matrix, too much missing data")
    reldiff <- conv.reldiff+1
    i <- 0
    while(i<niter && reldiff>conv.reldiff){
      svd0 <- svd(XF)
      Xnew <- X
      Xnew[Ina] <- (svd0$u[, 1:ncomp] %*%
                    diag(svd0$d[1:ncomp],nrow=length(svd0$d[1:ncomp])) %*%
                    t(svd0$v[,1:ncomp]))[Ina]
      diff <- max(abs(Xnew - XF))
      reldiff <- diff/max(abs(XF[Ina]))
      XF <- Xnew
      i <- i+1
    }
  }#if( all(!is.na(X)) ) ... else ...
  final.diff <- c(diff,reldiff,i,niter)
  names(final.diff) <- c("diff","rel.diff","n.iter","max.iter")
  return(list(svd=svd0, Xfill=XF, status=final.diff))
}##function SVDmiss

##' Function that computes smooth functions for a data matrix with missing
##' values, as described in Fuentes et. al. (2006), or does cross validation to
##' determine a suitable number of basis functions. The function uses
##' \code{\link{SVDmiss}} to complete the matrix and then computes smooth basis
##' functions by applying \code{\link[stats:smooth.spline]{smooth.spline}} to the
##' SVD of the completed data matrix.
##' 
##' \code{SVDsmoothCV} uses leave-one-column-out cross-validation; holding one column
##' out from \code{X}, calling \code{SVDsmooth}, and then regressing
##' the held out column on the resulting smooth functions. Cross-validation
##' statistics computed for each of these regressions include MSE, R-squared,
##' AIC and BIC. The weighted average (weighted by number of observations in the
##' colum) is then reported as CV-statistics.
##' 
##' @title Smooth Basis Functions for Data Matrix with Missing Values
##' 
##' @param X Data matrix, with missing values marked by \code{NA} (use
##'   \code{\link{createDataMatrix}}). Rows and/or columns 
##'   that are completely missing will be dropped (with a message), for the rows the
##'   smooths will be interpolated using
##'   \code{\link[stats:predict.smooth.spline]{predict.smooth.spline}}.
##' @param n.basis Number of smooth basis functions to compute, will be passed
##'   as \code{ncomp} to \code{\link{SVDmiss}}; for \code{SVDsmoothCV} a
##'   vector with the different number of basis functions to evaluate (including 0).
##' @param date.ind Vector giving the observation time of each row in
##'   \code{X}, used as \code{x} in \cr
##'   \code{\link[stats:smooth.spline]{smooth.spline}} when computing the smooth
##'   basis functions. If missing \code{\link{convertCharToDate}} is used to
##'   coerce the \code{rownames(X)}.
##' @param scale If \code{TRUE}, will use \code{\link[base:scale]{scale}} to scale
##'   \code{X} before calling \code{\link{SVDmiss}}.
##' @param niter,conv.reldiff Controls convergence, passed to \code{\link{SVDmiss}}.
##' @param df,spar The desired degrees of freedom/smoothing parameter for the
##'   spline, \cr see \code{\link[stats:smooth.spline]{smooth.spline}}
##' @param fnc If \code{TRUE} return a function instead of the trend-matrix, see
##'   Value below.
##' @param ... Additional parameters passed to \code{SVDsmooth}; i.e. \code{date.ind},
##'   \code{scale}, \code{niter}, \code{conv.reldiff}, \code{df}, \code{spar},
##'   and/or \code{fnc}.
##' 
##' @return Depends on the function:
##'   \item{SVDsmooth}{A matrix (if \code{fnc==FALSE}) where each column is a
##'                    smooth basis function based on the SVD of the completed
##'                    data matrix. The left most column contains the smooth of
##'                    the most important SVD. If \code{fnc==TRUE} a function
##'                    that will create the data matrix if called as
##'                    \code{fnc(date.ind)}, \code{fnc(1:dim(X)[1])}, or
##'                    \code{fnc(convertCharToDate( rownames(X) ))}.}
##'   \item{SVDsmoothCV}{A list of class \code{SVDcv} with components:
##'     \describe{
##'       \item{CV.stat,CV.sd}{\code{data.frame}s with mean and standard
##'             deviation of the CV statistics for each of the number of basis
##'             functions evaluated.} 
##'       \item{MSE.all,R2.all,AIC.all,BIC.all}{\code{data.frame}s with the
##'             individual MSE, R2, AIC, and BIC values for each column in the
##'             data matrix and for each number of basis functions evaluated.}
##'       \item{smoothSVD}{A list with \code{length(n.basis)} components. If
##'                        \code{fnc==FALSE} each component contains an array
##'                        where \code{smoothSVD[[j]][,,i]} is the result of
##'                        \code{SVDsmooth} applied to \code{X[,-i]} with
##'                        \code{n.basis[j]} smooth functions; if
##'                        \code{fnc==FALSE} each component contains a list
##'                        of functions as \code{smoothSVD[[j]][[i]]}.}
##'     }
##'   }
##' 
##' @references
##' M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
##'  Analyze Space-Time Processes in Statistical methods for spatio-temporal
##'  systems (B. Finkenstädt, L. Held, V. Isham eds.) 77-150
##' 
##' @author Paul D. Sampson and Johan Lindström
##'
##' @example Rd_examples/Ex_SVDsmooth.R
##' @family SVD for missing data
##' @family data matrix
##' @family SVDcv methods
##' @export
SVDsmooth <- function(X, n.basis=min(2,dim(X)[2]), date.ind=NULL, scale=TRUE,
                      niter=100, conv.reldiff=0.001, df=NULL, spar=NULL,
                      fnc=FALSE){

  ##date.ind contains NA/missing, use rownames of X
  if( missing(date.ind) || is.null(date.ind) || any(is.na(date.ind)) ){
    date.ind <- convertCharToDate( rownames(X) )
  }
  if( is.null(date.ind) ){
    date.ind <- 1:dim(X)[1]
  }
  ##check length od date.ind
  if( length(date.ind)!=dim(X)[1] ){
    stop("length(date.ind)!=dim(X)[1]")
  }
  ##check number of columns
  if( max(n.basis)>dim(X)[2] ){
    stop("Number of basis functions cannot exceed dim(X)[2]")
  }
  ##scale data to zero mean and unit variance
  if( scale ){
    X <- scale(X)
  }
  ##drop missing columns and rows
  Icol <- colSums(is.na(X))!=dim(X)[1]
  Irow <- rowSums(is.na(X))!=dim(X)[2]
  X.reduced <- X[Irow,Icol]
  if( any(!Icol) ){
    message( paste("Dropping column(s):", paste(which(!Icol), collapse=", ")) )
  }
  if( any(!Irow) ){
    message( paste("Smooth interpolated at", sum(!Irow), "time(s).") )
  }
  ##missing data SVD
  X.svd <- SVDmiss(X.reduced, niter=niter, ncomp=n.basis,
                   conv.reldiff=conv.reldiff)

  ##create function that compute trend matrix
  ##first compute splines for each basis-set
  spline <- lapply(1:n.basis,
                   function(j){
                     if( is.null(df) && is.null(spar) ){
                       return(smooth.spline(date.ind[Irow], X.svd$svd$u[,j]))
                     }else{
                       return(smooth.spline(date.ind[Irow], X.svd$svd$u[,j],
                                            df=df, spar=spar))
                     }
                   })
  ##second compute mean and scaling for each spline
  scale.spline <- lapply(spline, function(x){ c(mean(x$y), sd(x$y)) })
  ##define the function
  trend.fnc <- function(x=date.ind){
    X.comps <- matrix(NA, length(x), length(spline))
    for(i in 1:length(spline)){
      ##scale components to unit variance and zero mean
      X.comps[,i] <- scale(predict(spline[[i]], as.double(x))$y,
                           center=scale.spline[[i]][1],
                           scale=scale.spline[[i]][2])
      ##Ensure that components have alternating sign
      X.comps[,i] <- (-1)^i*X.comps[,i]*sign(X.comps[1,i])
    }
    ##add names
    rownames(X.comps) <- as.character(x)
    colnames(X.comps) <- paste("V", 1:length(spline), sep="")
    return( X.comps )
  }
  ##clean up the environment of the function (reduces overhead)
  tmp <- ls(environment(trend.fnc))
  tmp <- tmp[ !(tmp %in% c("date.ind", "spline", "scale.spline")) ]
  rm(tmp, pos=environment(trend.fnc))

  ##should we return the function or the matrix?
  if( fnc ){
    return( trend.fnc )
  }
  return( trend.fnc() )
}##function SVDsmooth

##' @rdname SVDsmooth
##' @export
SVDsmoothCV <- function(X, n.basis, ...){
  ##check number of columns
  if( max(n.basis)>(dim(X)[2]-1) )
    stop("Number of basis functions cannot exceed dim(X)[2]-1.")
   
  ##matrices that holds the individual CV-statistics
  CV.all <- array(NA, c(length(n.basis), dim(X)[2], 4))
  dimnames(CV.all) <- list(paste("n.basis", as.character(n.basis), sep="."),
                           colnames(X), c("MSE", "R2", "AIC", "BIC"))
  N <- matrix(colSums(!is.na(X)), 1, dim(X)[2])
  colnames(N) <- colnames(X)

  ##loop over the number of basis functions requested
  trend <- vector("list", length(n.basis))
  trend.fnc <- vector("list", length(n.basis))
  names(trend) <- names(trend.fnc) <- dimnames(CV.all)[[1]]
  for(i in 1:length(n.basis)){
    ##for each basis functions do leave one out CV
    trend[[i]] <- array(NA, c(dim(X)[1], n.basis[i], dim(X)[2]))
    trend.fnc[[i]] <- vector("list", dim(X)[2])
    for(j in 1:dim(X)[2]){
      if( n.basis[i]!=0 ){
        if(i==1 && j==1){
          tmp.trend <- SVDsmooth(X[,-j,drop=FALSE], n.basis[i], ...)
        }else{
          ##suppress message given above to avoid noisy CV output
          suppressMessages( tmp.trend <- SVDsmooth(X[,-j,drop=FALSE],
                                                   n.basis[i], ...))
        }
        ##did we return function or matrix?
        if( is.function(tmp.trend) ){
          trend.fnc[[i]][[j]] <- tmp.trend
          trend[[i]][,,j] <- tmp.trend()
        }else{
          trend[[i]][,,j] <- tmp.trend
        }
      }
      
      if( N[j] > (n.basis[i]+1) ){
        ##compute CV-error
        if( n.basis[i]!=0 ){
          tmp.lm <- lm(X[,j]~trend[[i]][,,j])
        }else{
          tmp.lm <- lm(X[,j]~1)
        }
        ##compute CV-statistics
        CV.all[i,j,"MSE"] <- mean(tmp.lm$residuals^2)
        CV.all[i,j,"R2"] <- summary(tmp.lm)$r.squared
        CV.all[i,j,"AIC"] <- extractAIC(tmp.lm)[2]
        CV.all[i,j,"BIC"] <- extractAIC(tmp.lm, k=log(N[j]))[2]
      }
    }
    #add names to the trend functions
    if( n.basis[i]!=0 ){
      dimnames(trend[[i]]) <- list(rownames(X), paste("V", 1:n.basis[i], sep=""),
                                   colnames(X))
    }else{
      dimnames(trend[[i]]) <- list(rownames(X), NULL, colnames(X))
    }
    names(trend.fnc[[i]]) <- colnames(X)
  }##for(i in 1:length(n.basis))

  weights <- N/sum(N)
  CV.stat <- apply(CV.all, c(1,3), weighted.mean, w=weights, na.rm=TRUE)
  CV.sd <- CV.stat
  for(i in 1:dim(CV.sd)[1]){
    for(j in 1:dim(CV.sd)[2]){
      CV.sd[i,j] <- weighted.mean((CV.all[i,,j] - CV.stat[i,j])^2,
                                  w=weights, na.rm=TRUE)
    }
  }
  CV.sd <- sqrt(CV.sd)

  ##return a class object
  out <- list(CV.stat=as.data.frame(CV.stat),
              CV.sd=as.data.frame(CV.sd),
              MSE.all=as.data.frame(t(CV.all[,,"MSE"])),
              R2.all=as.data.frame(t(CV.all[,,"R2"])),
              AIC.all=as.data.frame(t(CV.all[,,"AIC"])),
              BIC.all=as.data.frame(t(CV.all[,,"BIC"])),
              smoothSVD=switch(is.function(tmp.trend)+1, trend, trend.fnc))
  class(out) <- "SVDcv"
  return( out )
}##function SVDsmoothCV


##' A front end function for calling \code{\link{SVDsmooth}} (and
##' \code{\link{SVDsmoothCV}}), with either a \code{STdata} object
##' or vectors containing observations, dates and locations.
##' 
##' The function uses \code{\link{createDataMatrix}} to create a
##' data matrix which is passed to \code{\link{SVDsmooth}} (and
##' \code{\link{SVDsmoothCV}}). The output can be used as \cr
##' \code{STdata$trend = calcSmoothTrends(...)$trend}, or \cr
##' \code{STdata$trend = calcSmoothTrends(...)$trend.cv[[i]]}.
##' However, it is recommended to use \code{\link{updateTrend.STdata}}.
##'
##' @title Smooth Basis Functions for a STdata Object
##' 
##' @param STdata A \code{STdata}/\code{STmodel} data structure containing
##'   observations, see \code{\link{mesa.data.raw}}. Use either this or the \code{obs},
##'   \code{date}, and \code{ID} inputs.
##' @param obs A vector of observations.
##' @param date A vector of observation times.
##' @param ID A vector of observation locations.
##' @param subset A subset of locations to extract the data matrix for. A warning
##'   is given for each name not found in \code{ID}.
##' @param extra.dates Additional dates for which smooth trends should be
##'   computed (any duplicates will be removed).
##' @param n.basis Number of basis functions to compute, see
##'   \code{\link{SVDsmooth}}.
##' @param cv Also compute smooth functions using leave one out
##'   cross-validation, \cr see \code{\link{SVDsmoothCV}}.
##' @param ... Additional parameters passed to \code{\link{SVDsmooth}}
##'   and \code{\link{SVDsmoothCV}}; except \code{fnc}, which is always
##'   \code{TRUE}.
##' 
##' @return Returns a list with
##'   \item{trend}{A data.frame containing the smooth trends and the dates.
##'                This can be used as the \code{trend} in \code{STdata$trend}.}
##'   \item{trend.cv}{If \code{cv==TRUE} a list of data.frames; each one containing
##'                   the smooth trend obtained when leaving one site out.
##'                   Similar to \cr \code{SVDsmoothCV(data)$smoothSVD[[1]]}).}
##'   \item{trend.fnc,trend.fnc.cv}{Functions that produce the content of the above
##'        data.frames, see \code{\link{SVDsmooth}}.}
##' 
##' @author Johan Lindström and Paul D. Sampson
##' 
##' @example Rd_examples/Ex_calcSmoothTrends.R
##' @family SVD for missing data
##' @family STdata
##' @export
calcSmoothTrends <- function(STdata=NULL, obs=STdata$obs$obs,
                             date=STdata$obs$date, ID=STdata$obs$ID,
                             subset=NULL, extra.dates=NULL, n.basis=2,
                             cv=FALSE, ...){
  ##add extra dates
  date <- c(date, extra.dates)
  ##and expand obs and ID
  obs <- c(obs, rep(NA,length(extra.dates)) )
  ID <- c(ID, rep(ID[1],length(extra.dates)) )
  ##create data matrix (this removes duplicated extra.dates)
  data <- createDataMatrix(obs=obs, date=date, ID=ID, subset=subset)

  ##internal function
  extractTrend <- function(x){
    x <- as.data.frame(x)
    x$date <- convertCharToDate( rownames(x) )
    rownames(x) <- NULL
    return(x)
  }
  ##now let's do SVD
  data.comps.fnc <- SVDsmooth(data, n.basis, fnc=TRUE, ...)
  data.comps <- extractTrend(data.comps.fnc())
  ##and cross-validation
  if(cv){
    ##Message of interpolation already displayed above...
    suppressMessages( svd.cv <- SVDsmoothCV(data, n.basis, fnc=TRUE, ...) )
    svd.fnc <- svd.cv$smoothSVD[[1]]
    svd.tmp <- vector("list", length(svd.fnc))
    for(i in 1:length(svd.fnc)){
      svd.tmp[[i]] <- extractTrend( svd.fnc[[i]]() )
    }
    names(svd.tmp) <- names(svd.fnc)
  }else{
    svd.tmp <- NULL
    svd.fnc <- NULL
  }
  return( list(trend=data.comps, trend.cv=svd.tmp, trend.fnc=data.comps.fnc,
               trend.fnc.cv=svd.fnc) )
}##function calcSmoothTrends


##########################
## S3-METHODS FOR SVDcv ##
##########################

##' \code{\link[graphics:plot]{plot}} and
##' \code{\link[graphics:boxplot]{boxplot}} methods for class \code{SVDcv}.
##' Plots summary statistics for the cross-validation. Plots include
##' RMSE, R2, BIC, and scatter plots of BIC for each column.
##'
##' @title Plot and Boxplot cross-validation statistics for \code{SVDcv} object
##' @param x \code{SVDcv} object to plot.
##' @param y Which CV-statistic to plot. For pairs \code{"all"} implies
##'   \code{"BIC"}.
##' @param pairs \code{TRUE}/\code{FALSE} plot cross-validation statistics,
##'   or scatter plot of individual BIC:s.
##' @param sd \code{TRUE}/\code{FALSE} add uncertainty to each CV-statistic.
##' @param ... Additional parameters passed to \code{\link[graphics:plot]{plot}} or
##'   \code{\link[graphics:plot]{pairs}}.
##' @return Nothing
##'
##' @examples
##'   ##See SVDsmooth example
##' 
##' @author Johan Lindström
##' 
##' @family SVDcv methods
##' @family SVD for missing data
##' @method plot SVDcv
##' @export
plot.SVDcv <- function(x, y=c("all","MSE","R2","AIC","BIC"),
                       pairs=FALSE, sd=FALSE, ...){
  stCheckClass(x, "SVDcv", "'x'")

  ##we have to use y, cast to resonable name
  y <- match.arg(y)
  ##basis
  n.basis <- sapply(x$smoothSVD, dim)[2,]
  
  if( !pairs ){
    ##plot summary of cross-validation statistics
    if( y=="all" ){
      par(mfrow=c(2,2),mar=c(4,4,.5,.5))
      y <- c("MSE", "R2", "AIC", "BIC")
    }
    for(i in y){
      if( sd ){
        ylim <- range(x$CV.stat[[i]] + x$CV.sd[[i]],
                      x$CV.stat[[i]] - x$CV.sd[[i]])
      }else{
        ylim <- range(x$CV.stat[[i]])
      }
      args <- internalPlotFixArgs(list(...),
                                  default=list(xlab="n.basis", ylab=i),
                                  add=list(x=n.basis, y=x$CV.stat[[i]],
                                    type="l", ylim=ylim))
      do.call(plot, args)
      if( sd ){
        lines(n.basis, x$CV.stat[[i]] + x$CV.sd[[i]], lty=3)
        lines(n.basis, x$CV.stat[[i]] - x$CV.sd[[i]], lty=3)
      }
    }
  }else{
    ##plot the CV.stats for each column
    if( y=="all" ){ y <- "BIC" }
    y <- paste(y, "all", sep=".")
    pairs(x[[y]], panel=function(x,y){points(x,y); abline(0,1)}, ...)
  }
  return(invisible())
}##function plot.SVDcv

##' @rdname plot.SVDcv
##' @importFrom graphics boxplot
##' @method boxplot SVDcv
##' @export
boxplot.SVDcv <- function(x, y=c("all","MSE","R2","AIC","BIC"), ...){
  stCheckClass(x, "SVDcv", "'x'")

  ##we have to use y, cast to resonable name
  y <- match.arg(y)
  ##basis
  n.basis <- sapply(x$smoothSVD, dim)[2,]
  
  ##plot summary of cross-validation statistics
  if( y=="all" ){
    par(mfrow=c(2,2),mar=c(4,4,.5,.5))
    y <- c("MSE", "R2", "AIC", "BIC")
  }
  for(i in y){
    args <- internalPlotFixArgs(list(...),
                                default=list(ylab=i),
                                add=list(at=n.basis,
                                  x=x[[paste(i, "all", sep=".")]]))
    do.call(boxplot, args)
  }
  return(invisible())
}##function boxplot.SVDcv

##' \code{\link[base:print]{print}} and \code{\link[base:summary]{summary}}
##' methods for class \code{SVDcv}, prints cross-validation statistics.
##'
##' @title Print details for \code{SVDcv} object
##' @param x \code{SVDcv} object to print information for.
##' @param ... ignored additional arguments.
##' @return Nothing
##' 
##' @examples
##'   ##See SVDsmooth example
##' 
##' @author Johan Lindström
##' 
##' @family SVDcv methods
##' @family SVD for missing data
##' @method print SVDcv
##' @export
print.SVDcv <- function(x, ...){
  stCheckClass(x, "SVDcv", "'x'")
  cat("Result of SVDsmoothCV, average of CV-statistics:\n")
  print(x$CV.stat)
  return(invisible())
}##function print.SVDcv 

##' @rdname print.SVDcv
##' @param object \code{SVDcv} object to compute summary for.
##' @method summary SVDcv
##' @export
summary.SVDcv <- function(object, ...){
  stCheckClass(object, "SVDcv", "'object'")
  
  cat("Individual MSE:s by column:\n")
  print( apply(object$MSE.all, 2, summary) )
  cat("\nIndividual R2:s by column:\n")
  print( apply(object$R2.all, 2, summary) )
  cat("\nIndividual AIC:s by column:\n")
  print( apply(object$AIC.all, 2, summary) )
  cat("\nIndividual BIC:s by column:\n")
  print( apply(object$BIC.all, 2, summary) )
  return(invisible())
}##function summary.SVDcv 


############################################
## S3-METHODS that update temporal trends ##
############################################

##' Updates/sets the temporal trend for \code{STdata} or \code{STmodel}
##' objects. It also checks that the spatio-temporal covariate exists for all
##' dates in the trend, mainly an issue if \code{extra.dates!=NULL} adds
##' additional times at which to do predictions.
##'
##' If \code{n.basis} is given this will use \code{\link{calcSmoothTrends}} to
##' compute smoothed SVDs of data for use as temporal trends. If \code{fnc} is
##' given, \code{n.basis} is ignored and \code{fnc} should be a function that,
##' given a vector of dates, returns an object that can be coerced to a
##' data.frame with \emph{numeric} temporal trends; recall that an intercept is
##' \strong{always} added.
##'
##' For a \code{STmodel} object the new trend \emph{must have no more} components
##' than the existing trend; if a function is given \code{colnames} of the new
##' trend \emph{must match} those of the existing trend.
##' In both cases the returned \code{STdata} or \code{STmodel} object will have
##' both a \code{$trend} and \code{$trend.fnc} field.
##'
##' Function \code{updateSTdataTrend} is deprecated and will be removed in
##' future versions of the package.
##' 
##' @title Update Trend in \code{STdata} or \code{STmodel} Object
##' 
##' @param object A \code{STdata} or \code{STmodel} object, see
##'   \code{\link{mesa.data.raw}}. 
##' @param n.basis number of basis functions for the temporal trend
##' @param extra.dates Additional dates for which smooth trends should be
##'   computed (otherwise only those in \code{object$obs$date} are used);
##'   \emph{only} for \code{STdata}.
##' @param fnc Function that defines the trend, see Details and Example.
##' @param ... Additional parameters passed to \code{\link{calcSmoothTrends}}.
##' 
##' @return Returns a modfied version of the input, with an added/altered
##'   trend.
##' 
##' @example Rd_examples/Ex_updateTrend.R
##' 
##' @author Johan Lindström
##' @family STdata functions
##' @family SVD for missing data
##' @method updateTrend STdata
##' @export
updateTrend.STdata <- function(object, n.basis=0, fnc=NULL,
                               extra.dates=NULL, ...){
  ##check class belonging
  stCheckClass(object, "STdata", name="object")

  ##is there an existing trend
  if( !is.null(object$trend) ){
    message("Replacing existing trend.")
  }
  
  ##dates
  if( length(object$obs$date)==0 ){
    ##special case since c(..., Date) casts to numeric
    ##(an issue when object$obs$date is empty)
    dates <- sort(unique(extra.dates))
  }else{
    dates <- sort(unique(c(object$obs$date,extra.dates)))
  }
  
  if( !is.null(fnc) ){
    if( !missing(n.basis) ){
      warning("fnc defined, ignoring n.basis")
    }
    ##call function and update the trend
    object <- internalUpdateTrendFnc(object, fnc, dates)
  }else{
    object <- internalUpdateTrendSmooth(object, n.basis, dates, ...)
  }
  ##ensure that we have all dates (mainly checks that added dates are also in ST-covariate.
  trend.dates <- object$trend$date
  ST.dates <- convertCharToDate( rownames(object$SpatioTemporal) )
  if( any(!(trend.dates %in% ST.dates)) && !is.null(ST.dates) ){
    stop( paste("Trend dates not found in rownames(object$SpatioTemporal):",
                paste(trend.dates[!(trend.dates %in% ST.dates)],
                      collapse=", ")) )
  }
  ##return modified object
  return(object)
}##function updateTrend.STdata

##' @rdname updateTrend.STdata
##' @family STmodel functions
##' @method updateTrend STmodel
##' @export
updateTrend.STmodel <- function(object, n.basis=0, fnc=NULL, ...){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")

  ##can only use exactly the existing dates
  dates <- object$trend$date
  ##keep the old trend for comparison
  trend.old <- object$trend
  F.old <- object$F

  if( !is.null(fnc) ){
    if( !missing(n.basis) ){
      warning("fnc defined, ignoring n.basis")
    }
    object <- internalUpdateTrendFnc(object, fnc, dates)
  }else{
    object <- internalUpdateTrendSmooth(object, n.basis, dates, ...)
  }
  ##ensure that we have exatcly the same dates as before
  if( !isTRUE(all.equal(dates, object$trend$date)) ){
    stop("Dates for the trend have changed (this should NOT happen)")
  }
  ##check that we dont expand the size of the trend
  if( dim(trend.old)[2] < dim(object$trend)[2] ){
    stop("New trend must have no more components than existing trend.")
  }

  ##and update the F
  object$F <- internalSTmodelCreateF(object)
  
  Ind <- match( colnames(object$F), colnames(F.old) )
  if( any(is.na(Ind)) ){
    stop("colnames of new trend not found in the old: ",
         paste(colnames(object$F)[is.na(Ind)], collapse=", "))
  }
  ##Lets pick out components relevant to the new trend (LUR and beta-fields)
  object$LUR <- object$LUR[Ind]
  object$LUR.all <- object$LUR.all[Ind]
  object$LUR.list <- object$LUR.list[Ind]
  object$cov.beta$covf <- object$cov.beta$covf[Ind]
  object$cov.beta$nugget <- object$cov.beta$nugget[Ind]
  
  ##return modified object
  return(object)
}##function updateTrend.STdata

##' @rdname updateTrend.STdata
##' @export
updateTrend <- function(object, n.basis=0, fnc=NULL, ...){
  UseMethod("updateTrend")
}

##' @rdname updateTrend.STdata
##' @export
updateSTdataTrend <- function(object, n.basis=0, extra.dates=NULL, fnc=NULL, ...){
  message("Deprecated, use updateTrend instead. NOTE change in parameter order!")
  
  if( missing(n.basis) ){
    return( updateTrend.STdata(object=object, fnc=fnc,
                               extra.dates=extra.dates, ...) )
  }else{
    return( updateTrend.STdata(object=object, n.basis=n.basis, fnc=fnc,
                               extra.dates=extra.dates, ...) )
  }
}##function updateSTdataTrend


########################################
## Internal functions for updateTrend ##
########################################
internalUpdateTrendFnc <- function(object, fnc, dates){
    ##call function
    object$trend <- as.data.frame(fnc(dates))
    object$trend.fnc <- fnc
    ##drop rownames
    rownames(object$trend) <- NULL
    ##sanity checks
    if( !is.data.frame(object$trend) || dim(object$trend)[1]!=length(dates) ){
      stop("fnc(x) should return a data.frame with dim(.)[1]=length(x).")
    }
    if( any(!sapply(object$trend, is.numeric)) ){
      stop("fnc(x) should return a data.frame with 'numeric' columns")
    }
    if( "date" %in% names(object$trend) ){
      stop("Non of the columns of fnc(x) may be named 'date'.")
    }
    ##add a column of dates
    object$trend$date <- dates
    return( object )
}##function internalUpdateTrendFnc

internalUpdateTrendSmooth <- function(object, n.basis, dates, ...){
    if(n.basis==0){
      ##no temporal trend, just add dates for prediction
      object$trend <- data.frame( date=dates )
      ##constant, just hide the trend.fnc
      object$trend.fnc <- NULL
    }else{
      tmp <- calcSmoothTrends(STdata=object, n.basis=n.basis,
                              extra.dates=dates, ...)
      object$trend <- tmp$trend
      object$trend.fnc <- tmp$trend.fnc
    }
    return( object )
}##function internalUpdateTrendFnc


##############################################
## Internal functions for creation of trend ##
##############################################
internalCreateTrendFnc <- function(trend){
  X.trend <- trend
  Y.trend <- X.trend[,colnames(X.trend)!="date",drop=FALSE]
  X.trend <- X.trend$date
  trend.fnc <- function(x){
    out <- matrix(NA, length(x), dim(Y.trend)[2])
    colnames(out) <- colnames(Y.trend)
    for(i in 1:dim(out)[2]){
      out[,i] <- spline(X.trend, Y.trend[,i], xout=x)$y
    }
    return( out )
  }
  return( trend.fnc )
}##function internalCreateTrendFnc
