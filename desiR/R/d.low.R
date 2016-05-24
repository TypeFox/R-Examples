#' @title Low values are desirable
#'
#' @description Maps a numeric variable to a 0-1 scale such that low values are
#' desirable.
#'
#' @details Values less than \code{cut1} will have a high desirability. Values
#' greater than \code{cut2} will have a low desirability. Values between
#' \code{cut1} and \code{cut2} will have intermediate values.
#'
#' @param x Vector of numeric or integer values.
#' @param cut1,cut2 Values of the original data that define where the
#' desirability function changes.
#' @param des.min,des.max Minimum and maximum desirability values. Defaults to
#' zero and one, respectively.
#' @param scale Controls how steeply the function increases or decreases.
#'
#' @return Numeric vector of desirability values.
#' @seealso \code{\link{d.high}},  \code{\link{d.4pl}}
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(1000, mean=100, sd =5) # generate data
#' d <- d.low(x, cut1=90, cut2=110, scale=1)
#' 
#' # plot data
#' hist(x, breaks=30)
#' # add line
#' des.line(x, "d.low", des.args=c(cut1=90, cut2=110, scale=1))
#'
#' hist(x, breaks=30)
#' des.line(x, "d.low", des.args=c(cut1=90, cut2=110, des.min=0.1,
#' des.max=0.95, scale=1.5))

d.low <- function(x, cut1, cut2,  des.min = 0, des.max = 1, scale = 1){
  
  if(cut1 >= cut2) stop("cut1 must be less than cut2\n")
  if(des.min < 0 | des.min > 1) stop("des.min must be between zero and one\n")
  if(des.max < 0 | des.max > 1) stop("des.max must be between zero and one\n")
  if(scale <= 0) stop("scale must be greater than zero\n")
  
  y <- rep(NA,length(x))
  y <- ((x - cut2)/(cut1 - cut2))^scale
  y[x < cut1] <- 1
  y[x > cut2] <- 0

  # rescale:  des.min to des.max
  y <- (y * (des.max - des.min))  + des.min

  return(y)
}
