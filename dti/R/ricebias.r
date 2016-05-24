dwiRiceBias <- function(object,  ...) cat("No Rice Bias correction defined for this class:",class(object),"\n")

setGeneric("dwiRiceBias", function(object,  ...) standardGeneric("dwiRiceBias"))

setMethod("dwiRiceBias", "dtiData", function(object, 
                                             sigma = NULL, 
                                             ncoils = 1) {
  
  ## the method requires specification of sigma
  if (is.null(sigma) || sigma < 1) {
    cat("Please provide a value for sigma ... returning the original object!\n")
    return(object)
  }
  
  ## get all previous calls that generated this object
  args <- object@call
  
  ## test wether Bias correction has already been performed
  corrected <- FALSE
  for (i in 1:length(args)) { ## for all previous calls
    if (length(grep("dwiRiceBias", args[i][[1]])) > 0) {
      corrected <- TRUE
      cat("Rice bias correction already performed by\n")
      print(args[i][[1]])
      cat("\n ... returning the original object!\n")
    }
  }
  
  ## if Bias correction has not yet been performed, do it now
  if (!corrected) {
    ## replace data
    object@si <- array(ricebiascorr(object@si, sigma, ncoils), dim(object@si))
    ## add call top the list of calls
    object@call <- c(args, sys.call(-1))
  }
  
  ## return the object
  invisible(object)
})




ricebiascorr <- function(x, s = 1, ncoils = 1){
  
  ## get noncentrality parameter (ncp),
  ## mean (mu), standard deviation (sd) and variance (s2) of ncChi-distr for 
  ## ncp from 0 to 50 with step size 0.002
  varstats <- sofmchi(ncoils, 50, .002)
  
  ## standardize data x with sigma s, such that xt is the expectation of a ncChi-distr
  xt <- x/s
  
  ## cut at minimum value
  xt <- pmax(varstats$minlev, xt)
  
  ## find the indices of xt in varstats$mu ...  
  ind <- 
    findInterval(xt, varstats$mu, rightmost.closed = FALSE, all.inside = FALSE)
  ## ... use the corresponding ncp and rescale with s
  varstats$ncp[ind] * s
  
  ## done
}

