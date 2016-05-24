
#' @title Relationship between mean metabolic rate (\eqn{\bar{\epsilon}}) and abundance
#'  
#' @description \code{ebar} calculates the relationship between average metabolic rate of a species and that species' abundance. Also known as the Damuth relationship
#' @details
#' See examples.
#' 
#' @param x an object of class meteESF. 
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                abund=arth$count,
#'                power=arth$mass^(.75),
#'                minE=min(arth$mass^(.75)))
#' damuth <- ebar(esf1)
#' 
#' @return An object of class \code{meteRelaT}. The object contains a list with the following elements.
#' \describe{
#'   \item{\code{pred}}{predicted relationship}
#'   \item{\code{obs}}{observed relationship}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteDist, sad.meteESF, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#' @importFrom stats aggregate

ebar <- function(x) {
  if(is.na(x$state.var[3])) stop('must provide metabolic rate data or E0 to calculate power distributions')
  
  dat <- x$data$e
  
  if(is.null(dat)) {
    X <- NULL
  } else {
    X <- aggregate(list(n=x$data$n, e=x$data$e), list(s=x$data$s), sum)
    X$e <- X$e/X$n
    X <- X[, c('n', 'e')]
  }
  
  thr <- data.frame(n=1:min(x$state.var['N0'], max(X$n)), e=1 + 1/(1:min(x$state.var['N0'], max(X$n)) * x$La[2]))
  
  attr(X, 'source') <- 'empirical'
  attr(X, 'type') <- 'damuth'
  class(X) <- 'damuth'
  
  attr(thr, 'source') <- 'theoretical'
  attr(thr, 'type') <- 'damuth'
  class(thr) <- 'damuth'
  
  out <- list(obs=X, pred=thr)
  class(out) <- 'meteRelat'
  
  return(out)
}