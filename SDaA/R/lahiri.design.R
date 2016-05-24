# TODO: add example
#' Draw Samples Using Lahiri's Method
#' 
#' @param relsize vector of relative sizes of population PSUs 
#' @param n desired sample size
#' @param clnames vector of PSU names for population
#' @return clusters vector of n PSUs selected with replacement
#'         and with probability proportional to relsize
#' @author Sharon Lohr, slightly modified by Tobias Verbeke
#' @note Original code from Lohr (1999), p. 452 -- 453.
#' @references Lahiri, D. B. (1951). A method of sample selection 
#'             providing unbiased ratio estimates, \emph{Bulletin of
#'             the International Statistical Institute}, 33: 133 -- 140.
#' @export
lahiri.design <- function(relsize, n, clnames = seq(along = relsize)){
  maxrel <- max(relsize)
  sizeratio <- maxrel / mean(relsize)
  numpsu <- length(relsize)
  size <- 0
  clusters <- NULL
  while (size < n){
    ss <- ceiling((n - size) * sizeratio)
    temp <- sample(1:numpsu, ss, replace = TRUE)
    temp1 <- clnames[temp[relsize[temp] > runif(ss, min = 0, max = maxrel)]]
    clusters <- append(clusters, temp1[!is.na(temp1)])
    size <- length(clusters)
  }
  clusters[1:n]
}
