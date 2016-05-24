##' blockSeg fitting procedure
##'
##' Produce a blockwise estimation of a matrix.
##'
##' @param Y matrix of observations.
##' @param max.break a positive integer less than number of columns and number of rows.
##' By default, floor(min(ncol(Y),nrow(Y))/10+1).
##' @param max.var a positive integer less than number of columns times number of rows.
##' By default, ncol(Y)**2/2.
##' @param verbose logical. To display each step. By default TRUE.
##' @param Beta logical. To save each Beta associated at each lambda. By default FALSE (very heavy in memory space).
##'
##' @rdname blockSeg-proc
##'
##' @examples
##'  ## model parameters 
##' n <- 100 
##' K <- 5
##' mu <- suppressWarnings(matrix(rep(c(1,0),ceiling(K**2/2)), K,K))
##' Y <- rblockdata(n,mu,sigma=.5)$Y
##' res <- blockSeg(Y, 50)
##'
##' @export blockSeg
blockSeg <- function(Y, max.break=floor(min(ncol(Y),nrow(Y))/10+1), max.var = floor(ncol(Y)**2/2), verbose=TRUE, Beta=FALSE) {
  if (!(is.matrix(Y)||(class(Y)=="dgeMatrix"))){
    stop("Y must be the observations data (or a transformation)")
  }
  
  if (!is.numeric(max.break)){
    stop("max.break must be an integer between 1 and n")
  } else if ((max.break<=0)||(max.break>min(dim(Y)))||(length(max.break)!=1)||(floor(max.break)!=max.break)){
    stop("max.break must be an integer between 1 and n")
  }
  if (!is.numeric(max.var)){
    stop("max.var must be an integer between 1 and n1 times n2")
  } else if ((max.var<=0)||(max.var>(length(Y)))||(length(max.var)!=1)||(floor(max.var)!=max.var)){
    stop("max.var must be an integer between 1 and n1 times n2")
  }
  if (!is.logical(verbose)){
    stop("verbose must be logical : TRUE if you want the details of the procedure")
  }
  if (!is.logical(Beta)){
    stop("Beta must be logical : TRUE if you want to have the list of Beta")
  }
    
  out <- doLARS2D(
      R_Y         = as.matrix(Y),
      R_maxBreaks = max.break,
      R_maxVar    = max.var,
      R_verbose   = verbose,
      R_Beta      = Beta)
  
    return(new(Class = "blockSeg",
               Beta   = out$Beta  ,
               Lambda = out$Lambda,
               RowBreaks = out$RowBreaks,
               ColBreaks = out$ColBreaks,
               Actions   = out$Actions
               ))  
}
