#' Create the design plot of the functional data.
#'
#' This function will open a new device if not instructed otherwise.
#'
#' @param t a list of observed time points for functional data
#' @param obsGrid a vector of sorted observed time points
#' @param isColorPlot an option for colorful plot: 
#'                    TRUE: create color plot with color indicating counts
#'                    FALSE: create black and white plot with dots indicating observed time pairs
#' @param noDiagonal an option specifying plotting the diagonal design points:
#'                   TRUE:  remove diagonal time pairs
#'                   FALSE:  do not remove diagonal time pairs
#' @param ... Other arguments passed into \code{plot()}. 
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' CreateDesignPlot(sampWiener$tList, sort(unique(unlist(sampWiener$tList))))
#' @export

CreateDesignPlot = function(t, obsGrid = NULL, isColorPlot=TRUE, noDiagonal=TRUE, ...){
  
  if( class(t) != 'list'){
    stop("You do need to pass a list argument to 'CreateDesignPlot'!");
  }
  if( is.null(obsGrid)){
    obsGrid = sort(unique(unlist(t)))
  }
  
  args1 <- list( main="Design Plot", xlab= 'Observed time grid', ylab= 'Observed time grid')
  inargs <- list(...)
  args1[names(inargs)] <- inargs 
  
  
  # Check if we have very dense data (for visualization) on a regular grid
  if( (length(obsGrid) > 101) & all(sapply(t, function(u) identical(obsGrid, u)))){
    res = matrix(5, nrow = 101, ncol = 101)
    obsGrid = approx(x = seq(0,1,length.out = length(obsGrid)), y = obsGrid, 
                     xout = seq(0,1,length.out = 101))$y
  } else {
    res = DesignPlotCount(t, obsGrid, noDiagonal, isColorPlot)
  }
  
  oldpty <- par()[['pty']]
  par(pty="s")
  if(isColorPlot == TRUE){
    createColorPlot(res, obsGrid, args1)
  } else {
    createBlackPlot(res, obsGrid, args1)
  }
  par(pty=oldpty)
  
}

createBlackPlot = function(res, obsGrid, args1){
  
  if( is.null(args1$col)){
    args1$col = 'black'
  }
  if (is.null(args1$cex)){
    args1$cex = 0.33
  }
  if (is.null(args1$pch)){
    args1$pch = 19
  }
  
  u1 = as.vector(res)
  u2 = as.vector(t(res))
  t1 = rep(obsGrid, times = length(obsGrid) )
  t2 = rep(obsGrid, each = length(obsGrid)) 
  do.call( plot, c(args1, list( x = t1[u1 != 0], y = t2[u2 !=0] ) ) )  
  
}

createColorPlot = function(res, obsGrid, args1){
  
  res[res > 4] = 4;
  notZero <- res != 0
  nnres <- res[notZero]
  
  if ( is.null(args1$col) ){
    colVec <- c(`1`='black', `2`='blue', `3`='green', `4`='red')
    args1$col = colVec[nnres];
  } else {
    colVec = args1$col;
  }
  
  if ( is.null(args1$pch) ){
    pchVec <- rep(19, length(colVec))
    args1$pch = pchVec[nnres];
  } else {
    pchVec = args1$pch;
  }
  
  if ( is.null(args1$cex) ){
    cexVec <- seq(from=0.3, by=0.1, length.out=length(colVec))
    args1$cex <- cexVec[nnres]
  } else {
    cexVec <- args1$cex;
  }
  
  # pchVec <- rep(19, length(colVec))
  # names(pchVec) <- names(colVec)
  # if (!is.null(ddd[['pch']])){
  #   pchVec[] <- ddd[['pch']]
  # }
  
  # cexVec <- seq(from=0.3, by=0.1, length.out=length(colVec))
  # names(cexVec) <- names(colVec)
  # if (!is.null(ddd[['cex']])){
  #   cexVec[] <- ddd[['cex']]
  # }
  
  t1 = rep(obsGrid, times = length(obsGrid))
  t2 = rep(obsGrid, each = length(obsGrid)) 
  do.call( plot, c(args1, list( x = t1[notZero], y = t2[notZero]) ))
  
  if (!identical(unique(nnres), 1)){
    legend('right', c('1','2','3','4+'), pch = pchVec, col=colVec, pt.cex=1.5, title = 'Count',bg='white' )
  }
  
}


