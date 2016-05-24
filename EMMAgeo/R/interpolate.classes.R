#' Function to interpolate classes.
#' 
#' This function interpolates grain-size classes, either to higher or to lower
#' resolution.
#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param boundaries.in Numeric vector with class boundaries of the input data.
#' @param boundaries.out Numeric vector with class boundaries of the output
#' data.
#' @param method Logical scalar, interpolation method, one out of "linear"
#' (linear interpolation), "fmm" (cubic spline), "natural" (natural spline),
#' "periodic" (periodic spline). Default is \code{"natural"}.
#' @param fixed.start Logical scalar, specifying if the outer boundaries should
#' be set to the same values as in the original matrix, default is \code{TRUE}.
#' This may become necessary to avoid interpolation errors, see example.
#' @return Numeric matrix with interpolated class values.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{approx}}, \code{\link{spline}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data
#' data(X, envir = environment())
#' classes.in <- seq(from = 1, to = 10, length.out = ncol(X))
#'   
#' ## Example 1 - decrease the class numbers
#' ## define number of output classes
#' classes.out <- seq(1, 10, length.out = 20)
#' 
#' ## interpolate the data set
#' Y <- interpolate.classes(X = X, 
#'                          boundaries.in = classes.in, 
#'                          boundaries.out = classes.out,
#'                          method = "linear")
#' 
#' ## show original vs. interpolation for first 10 samples
#' plot(NA, xlim = c(1, 10), ylim = c(0, 11))
#' for(i in 1:10) {
#'   lines(classes.in, X[i,] * 20 + i)
#'   lines(classes.out, Y[i,] * 20 + i, col = 2)
#' }
#' 
#' ## Example 2 - increase the class numbers
#' ## define number of output classes
#' classes.out <- seq(1, 10, length.out = 200)
#' 
#' ## interpolate the data set
#' Y <- interpolate.classes(X = X, 
#'                          boundaries.in = classes.in, 
#'                          boundaries.out = classes.out)
#' 
#' ## show original vs. interpolation for first 10 samples
#' plot(NA, xlim = c(1, 10), ylim = c(0, 11))
#' for(i in 1:10) {
#'   lines(classes.in, X[i,] * 20 + i)
#'   lines(classes.out, Y[i,] * 20 + i, col = 2)
#' }
#' 
#' @export interpolate.classes
interpolate.classes <- function(
  X,
  boundaries.in,
  boundaries.out,
  method = "natural",
  fixed.start = TRUE
){

  ## check input data
  if(min(boundaries.in, na.rm = TRUE) > min(boundaries.out, na.rm = TRUE) |
       max(boundaries.in, na.rm = TRUE) < max(boundaries.out, na.rm = TRUE)) {
    stop("Output boundary range is wider than input boundary range!")
  }
  
  ## transform data structure of X
  if(is.matrix(X) == FALSE) {X <- t(as.matrix(X))}
  
  ## create cumulative class sums matrix
  Xcum <- t(apply(X, 1, cumsum))

  ## create output matrices
  Y <- matrix(nrow = nrow(X), ncol = length(boundaries.out))
  Ycum <- Y

  if(method == "linear") {
  ## linear interpolation
    for(i in 1:nrow(X)) {
      Ycum[i,] <- approx(x = boundaries.in,
                         y = Xcum[i,],
                         xout = boundaries.out,
                         rule = 1)$y
    }
  } else {
  ## spline interpolation
    for(i in 1:nrow(X)) {
      Ycum[i,] <- spline(x = boundaries.in, 
                         y = Xcum[i,], 
                         xout = boundaries.out, 
                         method = method)$y
    }
  }
  
  ## reverse cumulative sums to individual class sums
  Y[,1] <- Ycum[,1]
  Y[,2:ncol(Y)] <- Ycum[,2:ncol(Ycum)] - Ycum[,1:(ncol(Ycum) - 1)]
        
  ## rescale interpolated data to original limits
  Y.norm <- Y
  Y.scaled <- Y
  for(i in 1:nrow(X)){
    Y.norm[i,] <- (Y[i,] - min(Y[i,])) / max(Y[i,] - min(Y[i,]))
    Y.scaled[i,] <- Y.norm[i,] *  (max(X[i,]) - min(X[i,])) + min(X[i,])
  }

  ## optionally, assign original values to first class
  if(fixed.start == TRUE) {
    Y.scaled[,1] <- X[,1]
  }
  
  ## return result matrix
  return(Y.scaled)
}