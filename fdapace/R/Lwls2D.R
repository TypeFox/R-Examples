#' Two dimensional local linear kernel smoother.
#'
#' Two dimensional local weighted least squares smoother. Only local linear smoother for estimating the original curve is available (no higher order, no derivative). 
#' @param bw A scalar or a vector of length 2 specifying the bandwidth.
#' @param kern Kernel used: 'gauss', 'rect', 'gausvar', 'epan' (default), 'quar'.
#' @param xin An n by 2 dataframe or matrix of x-coordinate.
#' @param yin A vector of y-coordinate.
#' @param win A vector of weights on the observations. 
#' @param xout1 a p1-vector of first output coordinate grid. Defaults to the input gridpoints if left unspecified.
#' @param xout2 a p2-vector of second output coordinate grid. Defaults to the input gridpoints if left unspecified.
#' @param xout alternative to xout1 and xout2. A matrix of p by 2 specifying the output points (may be inefficient if the size of \code{xout} is small).
#' @param crosscov using function for cross-covariance estimation (Default: FALSE)
#' @param subset  a vector with the indeces of x-/y-/w-in to be used (Default: NULL)
#' @return a p1 by p2 matrix of fitted values if xout is not specified. Otherwise a vector of length p corresponding to the rows of xout. 
#' @export

# Uses Pantelis' cpp code.
Lwls2D <- function(bw, kern='epan', xin, yin, win=NULL, xout1=NULL, xout2=NULL, xout=NULL, subset=NULL, crosscov = FALSE) {
  userNumCores = NULL
  if (length(bw) == 1){
    bw <- c(bw, bw)
  }
  if (!is.matrix(xin) ||  (dim(xin)[2] != 2) ){
    stop('xin needs to be a n by 2 matrix')
  }
  # xin <- matrix(xin, ncol=2) # This causes unexcepted/wrong results.
  
  if (is.null(win)){
    win <- rep(1, nrow(xin))
  }  
  if (!is.null(subset)) {
    xin <- xin[subset, ]
    yin <- yin[subset]
    win <- win[subset]
  }
  
  if (!is.null(xout1) && !is.null(xout2) && !is.null(xout)) {
    stop('Either xout1/xout2 or xout should be specified, but not both.')
  }
  
  if (is.null(xout1)) 
    xout1 <- sort(unique(xin[, 1]))
  
  if (is.null(xout2)) 
    xout2 <- sort(unique(xin[, 2]))
  
  # For passing numerics into the cpp smoother.
  storage.mode(bw) <- 'numeric'
  storage.mode(xin) <- 'numeric'
  storage.mode(yin) <- 'numeric'
  storage.mode(win) <- 'numeric'
  storage.mode(xout1) <- 'numeric'
  storage.mode(xout2) <- 'numeric'
  if (!is.null(xout))
    storage.mode(xout) <- 'numeric' 
  
  if (!crosscov){
    if( is.null(userNumCores) ){
      ret <- Rmullwlsk(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
    } else {
      
      if ( length(xout1) < userNumCores){  
        warning('You have allocated more cores than grid-points to evaluate. Nice... \n (We will use only as many cores as grid-points)')
        userNumCores =  length(xout1)
      }       
      if ( length(xout1)*0.01 < userNumCores){   
        warning('Less than 10000 points to be processed per core. Probably you will get no speed-up.')
      }
      if ( !all.equal( xout1, xout2)){
        stop('Both xout1 and xout2 must be the same to use multicore.')
      }
      
      parSmooth2 <- function(bw, xin, yin, kernel_type, win, nc, xout){ 
        N = length(xout) 
        breakPoints =   sort( N-  round(sqrt((1:(nc-1))/nc) * N))
        m2 <- Matrix::Matrix(0, nrow = N, N, sparse = TRUE)
        i_indx = c( breakPoints,  N ) # 1 3 7
        j_indx = c( 1, breakPoints+1) # 1 2 4 
        for (ij in 1:nc){ 
          q = i_indx[ij]
          p = j_indx[ij]
          m2[(p:q) , (p:N)] =  Rmullwlsk(bw, xin, cxxn= yin, xgrid=xout[p:N], ygrid=xout[p:q], 
                                         kernel_type=kernel_type, win=rep(1,length(yin)), FALSE, FALSE)
        } 
        m3 = as.matrix(m2) + t(as.matrix(m2))
        diag(m3) = 0.5 * diag(m3)
        return(m3)
      }
      ret = parSmooth2(bw = bw, xin = t(xin), yin = yin, kernel_type = 'epan', win = win, nc = 3,  xout = xout1) 
    } 
  } else {
    ret <- RmullwlskCC(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
  }
  
  if (!is.null(xout)) {
    ret <- interp2lin(xout1, xout2, ret, xout[, 1], xout[, 2])
  }
  
  return(ret)
}
