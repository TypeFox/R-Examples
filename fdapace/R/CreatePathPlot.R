#' Create the sample path plot based on the results from FPCA().
#'
#' @param fpcaObj Returned object from FPCA().
#' @param subset A vector of indices or a logical vector for subsetting the
#' observations.
#' @param k The number of components to reconstruct the sample paths.
#' @param inputData A list of length 2 containing the sparse/dense
#' (unsupported yet) observations. \code{inputData} needs to contain two
#' fields: \code{t} for a list of time points and \code{y} for a list of
#' observations. Default to the `inputData` field within `fpcaObj`.
#' @param showObs Whether to plot the original observations for each subject.
#' @param derOptns A list of options to control derivation parameters; see `fitted.FPCA'. (default = NULL)
#' @param ... other arguments passed into matplot for plotting options
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan',
#'             verbose=TRUE))
#' CreatePathPlot(res, subset=1:5)
#' @export

CreatePathPlot = function(fpcaObj, subset, k=NULL, inputData=fpcaObj[['inputData']], 
                          showObs=!is.null(inputData),  derOptns = NULL, ...){
  
  n <- dim(fpcaObj[['xiEst']])[1]
  
  if (!is.null(derOptns[['p']]) && derOptns[['p']] >= 1 && missing(showObs)) {
  # makes no sense to show original observations with derivatives.
    showObs <- FALSE 
  }
  if (!is.null(inputData)) {
    if (!all(c('t', 'y') %in% names(inputData))) {
      stop('inputData does not contain the required fields `t` and `y`')
    }
  } 
  if (showObs) {
    if (is.null(inputData)) {
      stop('Cannot show the sparse observations due to unspecified input data')
    } else {
      if (length(inputData[['t']]) != n)
        stop('length of inputData mismatches that in fpcaObj')
    }
  }
  
  if (missing(subset)) {
    subset <- seq_len(n)
  }
  # browser()
  workGrid <- fpcaObj[['workGrid']]
  fit <- fitted(fpcaObj, k=k, derOptns = derOptns)[subset, , drop=FALSE]
  
  args1 <- list( xlab= 's', ylab= ' ')                    
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  
  # make a matrix with NAs for the sparse observations.
  maxN_i <- max(sapply(inputData[['t']][subset], length))
  obst <- sapply(inputData[['t']][subset], function(x) c(x, rep(NA, maxN_i - length(x))))
  obsy <- sapply(inputData[['y']][subset], function(x) c(x, rep(NA, maxN_i - length(x))))
  
  #matplot(obst, obsy, type='p',...)
  #args2 = list (x = obst, y = obsy, type='p' )
  
  if( showObs ){
    do.call(matplot, c(list(x=obst, y=obsy, type='p'), args1))
    do.call(matplot, c(list(x=workGrid, y=t(fit), type='l', add=TRUE)))
  } else {
      do.call(matplot, c(list(x=workGrid, y=t(fit), type='l'), args1))
  }
   
  
  invisible()
}
