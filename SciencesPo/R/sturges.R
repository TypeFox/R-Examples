#' @encoding UTF-8
#' @title  Calculate breaks
#'
#' @description Calculate breaks according to the Herbert Sturges' (1926) formula
#'
#' @param x A vector of count values.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @references
#' Sturges, H. (1926) \emph{The choice of a class-interval.} J. Amer. Statist. Assoc., 21, 65-66.
#' @examples
#' set.seed(51)
#' y <- sample(100)
#' sturges(y)
#'
#' @export
sturges <- function(x)
{
  N    = length(x)
  K    = 1+log2(N)
  range = (ceiling(max(x))-floor(min(x)))/K
  cat("N = ",N,", Min = ",min(x),", Max = ",max(x),", Intervals = ",round(K,2),", Range = ",round(range,2),"\n",sep="")
  range = round(range)
  K    = round(K)
  classes =floor(min(x))+(0:K)*range
  if (max(classes ) < max(x)) { classes [length(classes )+1]=max(classes )+range }
  if (length(which(classes  >= max(x))) > 1) { classes  = classes [-length(classes )] }
  classes
}### end -- sturges function
NULL
