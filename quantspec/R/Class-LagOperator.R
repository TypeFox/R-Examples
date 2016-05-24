#' @include generics.R
NULL

################################################################################
#' Interface Class to access different types of operators on time series.
#'
#' \code{LagOperator} is an S4 class that provides a common interface to
#' implementations of an operator \eqn{\Gamma(Y)}{Gamma(Y)} which is calculated on 
#' all pairs of observations \eqn{(Y_0,Y_k)}{(Y0,Yk)} with lag smaller than maxLag 
#'
#' Currently one implementation is available:
#'     (1) \code{\link{ClippedCov}}.
#' 
#' Currently, the implementation of this class allows only for the analysis of
#' univariate time series.  
#'
#' @name LagOperator-class
#' @aliases LagOperator
#'
#' @keywords S4-classes
#'
#' @slot values an array of dimension \code{c(maxLag,length(levels.1),length(levels.2))}
#' 			 containing the values of the operator.
#' @slot Y is the time series the operator shall be applied to
#' @slot maxLag maximum lag between two observations
#' @slot levels a vector of numerics that determines the levels of the operator
#' @slot isRankBased A flag that is \code{FALSE} if the determined \code{values}
#'                     are based on the original time series and \code{TRUE} if it
#'                     is based on the ranks.
#' @slot positions.boot An object of type \code{\link{BootPos}},
#'                      that is used to determine the block bootstrapped
#'                      replicates of the time series.
#' @slot B Number of bootstrap replications to perform.
#'
################################################################################


setClass(
  Class = "LagOperator",
  representation=representation(
        Y = "numeric",
        values = "array",
        maxLag = "numeric",
        levels.1 = "numeric",
        levels.2 = "numeric",
        isRankBased = "logical",
        positions.boot = "BootPos",
        B = "numeric"
    )
)

################################################################################
#' Get attribute \code{values} from a \code{LagOperator}.
#'
#' @name getValues-LagOperator 
#' @aliases getValues,LagOperator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagOperator} from which to get the \code{values}.
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values    
#' @return Returns the \code{values} attribute.
################################################################################
setMethod(f = "getValues",
          signature = "LagOperator",
          definition = function(object,levels.1,levels.2) {
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(levels.1)) {
              levels.1 <- object@levels.1
            }
            if (!hasArg(levels.2)) {
              levels.2 <- object@levels.2
            }
            # end: workaround
            
            
            # Select columns
        
            c.1.pos <- match(levels.1,object@levels.1)
            c.2.pos <- match(levels.2,object@levels.2)
            
            if(is.na(c.1.pos[1])){
              stop("no 'values' for 'levels.1' requested were found")
            }
            if(is.na(c.2.pos[1])){
              stop("no 'values' for 'levels.2' requested were found")
            }
            
            if(!(length(c.1.pos)==length(levels.1))){
              warning("not all requested 'levels.1' were found")
            }
            
            if(!(length(c.2.pos)==length(levels.2))){
              warning("not all requested 'levels.2' were found")
            }
            
            ln.1 = length(levels.1)
            ln.2 = length(levels.2)
            ln = dim(object@values)[1]
            return(array(object@values[,c.1.pos,c.2.pos,],dim = c(ln,ln.1,ln.2,getB(object)+1)))
          }
)

################################################################################
#' Get \code{maxLag} from a \code{\link{LagOperator}} object.
#'
#' @name getMaxLag-LagOperator
#' @aliases getMaxLag,LagOperator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagOperator} of which to get the \code{maxLag}
#' @return Returns the attribute \code{maxLag} that's a slot of \code{object}.
################################################################################
setMethod(f = "getMaxLag",
    signature = signature("LagOperator"),
    definition = function(object) {
      return(object@maxLag)
    }
)

################################################################################
#' Get attribute \code{levels} from a \code{LagOperator}.
#'
#' If the optional parameter \code{j} is supplied, then the \code{j}th vector of
#' levels will be returned, a list with all vectors otherwise.
#'
#' @name getLevels-LagOperator
#' @aliases getLevels,LagOperator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagOperator} from which to get the \code{levels}.
#' @param j Index pointing to a set of levels in the list; optional.
#'
#' @return Returns levels attribute, as a vector of real numbers.
################################################################################
setMethod(f = "getLevels",
  signature = "LagOperator",
  definition = function(object,j) {
    if (missing("j")) {
      levels = list(object@levels.1,object@levels.2)
      names(levels) = c("levels.1","levels.2")
      return(levels)
    } else {
      if (!(j==1 | j==2)) {
        stop("Index needs to be either 1 or 2.")
      } else {
        if(j==1){return(object@levels.1)}
        if(j==2){return(object@levels.2)}
      }
    }
  }
)


################################################################################
#' Get \code{isRankBased} from a \code{\link{LagOperator}} object
#'
#' @name getIsRankBased-LagOperator
#' @aliases getIsRankBased,LagOperator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagOperator} of which to get the \code{isRankBased}
#'
#' @return Returns the attribute \code{isRankBased} that's a slot of \code{object}.
################################################################################
setMethod(f = "getIsRankBased",
    signature = signature("LagOperator"),
    definition = function(object) {
      return(object@isRankBased)
    }
)

################################################################################
#' Get \code{B} from a \code{\link{LagOperator}} object.
#'
#' @name getB-LagOperator
#' @aliases getB,LagOperator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagOperator} of which to get the \code{B}
#' @return Returns the attribute \code{B} that's a slot of \code{object}.
################################################################################
setMethod(f = "getB",
    signature = signature("LagOperator"),
    definition = function(object) {
      return(object@B)
    }
)

################################################################################
#' Get associated \code{\link{BootPos}} from a
#' \code{\link{LagOperator}}.
#'
#' @name getBootPos-LagOperator
#' @aliases getBootPos,LagOperator-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{LagOperator} from which to get the
#'                \code{\link{BootPos}}.
#' @return Returns the \code{\link{BootPos}} object associated.
################################################################################
setMethod(f = "getBootPos",
    signature = "LagOperator",
    definition = function(object) {
      return(object@positions.boot)
    }
)


setMethod(f = "show",
    signature = "LagOperator",
    definition = function(object) {

      maxLag <- getMaxLag(object)
      K1 <- length(getLevels(object, 1))
      K2 <- length(getLevels(object, 2))
      B <- getB(object)
      
      values <- getValues(object, levels.1 = getLevels(object, 1), levels.2 = getLevels(object, 2))
      values <- array(values, dim = c(maxLag+1, K1, K2, B+1))
      
      cat(paste("\n",class(object)," (maxLag = ",maxLag,", K1 = ",K1,", K2 = ",K2,", B+1 = ",B+1,")\n", sep=""))

      if (K1 <= 10) {
        cat("Levels 1   : ", getLevels(object, 1),"\n")
      } else {
        cat("Levels 1   : ", getLevels(object, 1)[1:5],"..",object@levels[[1]][(K1-4):K1],"\n")
      }
      
      if (K2 <= 10) {
        cat("Levels 2   : ", getLevels(object, 2),"\n")
      } else {
        cat("Levels 2   : ", getLevels(object, 2)[1:5],"..",object@levels[[2]][(K2-4):K2],"\n")
      }
      
      cat("\nValues:\n")
      
      resultMatr <- matrix(nrow=maxLag+1, ncol=K1*K2)
      cn <- rep(0,K1*K2)
      for (k1 in 1:K1) {
        for (k2 in 1:K2) {
          resultMatr[,k1+(k2-1)*K1] <- values[,k1,k2, 1]
          cn[k1+(k2-1)*K1] <- paste(getLevels(object, 1)[k1],"/", getLevels(object, 2)[k2], sep="")
        }
      }
      nrowShow <- min(10,nrow(resultMatr))
      ncolShow <- min(5,ncol(resultMatr))
      
      res <- apply(resultMatr[1:nrowShow, 1:ncolShow, drop=FALSE], c(1,2),
          #function(x){complex(real=round(Re(x),3), imaginary=round(Im(x),3))})
          function(x){round(x,3)})
      rownames(res) <- (1:(maxLag+1))[1:nrowShow]
      colnames(res) <- cn[1:ncolShow]
      
      show(res)
    }
)

################################################################################
#' Plot the values of the \code{\link{LagOperator}}.
#'
#' Creates a \code{K} x \code{K} plot (where \code{K} is the length of the \code{levels} parameter)
#' showing the values of the \code{\link{LagOperator}}. The plots below the diagonal show the positive
#' Lags and the plots above display the negative ones. 
#'
#' @name plot-LagOperator
#' @aliases plot,LagOperator,ANY-method
#' @export
#'
#' @param x The \code{\link{LagOperator}} to plot.
#' @param maxLag maximum Lag that should be displayed. It defaults to the
#'               maximum number of Lags available but usually a smaller number
#'               yields a more informative result.
#' @param ratio quotient of width over height of the subplots; use this
#'               parameter to produce landscape or portrait shaped plots.
#' @param widthlab width for the labels (left and bottom); default is
#'                  \code{lcm(1)}, cf. \code{\link[graphics]{layout}}.
#' @param xlab label that will be shown on the bottom of the plots; can be
#'              an expression (for formulas), characters or \code{NULL} to
#'              force omission (to save space).
#' @param ylab label that will be shown on the left side of the plots;
#'              can be an expression (for formulas), characters or
#'              \code{NULL} to force omission (to save space).
#' @param levels a set of levels for which the values are to be plotted.
#'
################################################################################
setMethod(f = "plot",
    signature = signature(x = "LagOperator"),
    definition = function(x, levels=intersect(x@levels.1, x@levels.2), maxLag = maxLag, widthlab = lcm(1),ratio = 3/2, xlab = expression(omega/2*pi), ylab = NULL) {
      def.par <- par(no.readonly = TRUE) # save default, for resetting...
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(maxLag)) {
        maxLag = x@maxLag
      }
      if (!hasArg(levels)) {
        levels <- intersect(x@levels.1, x@levels.1)
      }
      if (!hasArg(widthlab)){
        widthlab = lcm(1)
      }
      if (!hasArg(xlab)) {
        xlab <- "Lag"
      }
      if (!hasArg(ylab)) {
        ylab <- NULL
      }
      if (!hasArg(ratio)) {
        ratio <- 3/2
      }
      # end: workaround
      
      if (length(levels) == 0) {
        stop("There has to be at least one level to plot.")
      }
      if (maxLag <= 0){
        stop("maxLag has to be a positive integer.")
      }
      if(maxLag > x@maxLag){
        maxLag = x@maxLag
        warning("maxLag too large, set to maximum")
      }
      
        K <- length(levels)
        values <- getValues(x, levels.1=levels, levels.2=levels)
      
        p <- K
        M <- matrix(1:p^2, ncol=p)
        M.heights <- rep(1,p)
        M.widths  <- rep(ratio,p)
        
        # Add places for tau labels
        M <- cbind((p^2+1):(p^2+p),M)
        M.widths <- c(widthlab,M.widths)
        M <- rbind(M,c(0,(p^2+p+1):(p^2+2*p)))
        M.heights <- c(M.heights, widthlab)
        
        i <- (p^2+2*p+1)
        # Add places for labels
        if (length(xlab)>0) {
          M.heights <- c(M.heights, widthlab)
          M <- rbind(M,c(rep(0,length(M.widths)-p),rep(i,p)))
          i <- i + 1
        }
        
        if (length(ylab)>0) {
          M <- cbind(c(rep(i,p),rep(0,length(M.heights)-p)),M)
          M.widths <- c(widthlab,M.widths)
        }
        
        nf <- layout(M, M.widths, M.heights, TRUE)
        
        par(mar=c(2,2,1,1))
        
      
      for(i in 1:K){
        for(j in 1:K){
          if (i <= j){
          plot(values[1:maxLag,i,j,1],type = "h",xlab = "",ylab = "")
          abline(h=0)
          }else{
            plot((-(maxLag-1):0),values[maxLag:1,i,j,1],type = "h",xlab = "",ylab = "")
            abline(h=0)
          }
        }
      }
      p = K
      par(mar=c(0,0,0,0))
      for (i in 1:p) {
        plot.new()
        text(0.5,0.5,substitute(paste(tau[1],"=",k),list(k=levels[i])), srt=90)
      }
      
      for (i in 1:p) {
        plot.new()
        text(0.5,0.5,substitute(paste(tau[2],"=",k),list(k=levels[i])))
      }
      if (length(xlab)>0) {
        plot.new()
        text(0.5, 0.5, xlab)
      }
      if (length(ylab)>0) {
        plot.new()
        text(0.5, 0.5, ylab, srt=90)
      }
      
   par(def.par)
  }
)