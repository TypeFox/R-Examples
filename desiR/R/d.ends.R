#' @title Extreme (both high and low) values are desirable
#'
#' @description Maps a numeric variable to a 0-1 scale such that values at the
#' ends of the distribution are desirable.
#'
#' @details Values less than \code{cut1} and greater than \code{cut4} will have
#' a high desirability. Values between \code{cut2} and \code{cut3} will have a
#' low desirability. Values between \code{cut1} and \code{cut2} and between
#' \code{cut3} and \code{cut4} will have intermediate values. This function is
#' useful when the data represent differences between groups; for example, log2
#' fold-changes in gene expression. In this case, both high an low values are of
#' interest.
#'
#' @param x Vector of numeric or integer values.
#' @param cut1,cut2,cut3,cut4 Values of the original data that define where the
#' desirability function changes.
#' @param des.min,des.max Minimum and maximum desirability values. Defaults to
#' zero and one, respectively.
#' @param scale Controls how steeply the function increases or decreases.
#'
#' @return Numeric vector of desirability values.
#' @seealso \code{\link{d.central}}
#' @examples
#' set.seed(1)
#' x <- rnorm(1000, mean=100, sd =5) # generate data
#' d <- d.ends(x, cut1=90, cut2=95, cut3=105, cut4=110, scale=1)
#'
#' # plot data
#' hist(x, breaks=30)
#' # add line
#' des.line(x, "d.ends", des.args=c(cut1=90, cut2=95, cut3=105,
#' cut4=110, scale=1))
#'
#' hist(x, breaks=30)
#' des.line(x, "d.ends", des.args=c(cut1=90, cut2=95, cut3=105,
#' cut4=110, des.min=0.1, des.max=0.95, scale=1.5))


d.ends <- function(x, cut1, cut2, cut3, cut4, des.min = 0,
                     des.max = 1, scale = 1){

  if(cut1 >= cut2) stop("cut1 must be less than cut2\n")
  if(cut2 >= cut3) stop("cut2 must be less than cut3\n")
  if(cut3 >= cut4) stop("cut3 must be less than cut4\n")
  if(des.min < 0 | des.min > 1) stop("des.min must be between zero and one\n")
  if(des.max < 0 | des.max > 1) stop("des.max must be between zero and one\n")
  if(scale <= 0) stop("scale must be greater than zero\n")
  
  y <- rep(NA, length(x))
  for (i in 1:length(x)){
    if (is.na(x[i])) next
    if (x[i] <= cut1 | x[i] >= cut4)  y[i] <- 1
    if (x[i] >= cut2 & x[i] <= cut3) y[i] <- 0
    if (x[i] > cut1 & x[i] < cut2) y[i] <- ((x[i] - cut2)/(cut1 - cut2))^scale
    if (x[i] > cut3 & x[i] < cut4) y[i] <- ((x[i] - cut3)/(cut4 - cut3))^scale
  }

  # rescale:  des.min to des.max
  y <- (y * (des.max - des.min))  + des.min
  
  return(y)
}
