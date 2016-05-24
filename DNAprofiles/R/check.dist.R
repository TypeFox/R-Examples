#' Check if dist is properly specified for use with sim.q, exact.q
#' 
#' @param dist list of two numeric vectors \code{x} and \code{fx}: \code{x} lists the outcomes and \code{fx} the probabilities
#' @return list the supplied \code{dist}.
#' @details Several functions deal with computing a probability \deqn{P(X_1 X_2 \cdots X_n > t),}
#' where \eqn{X_1,\ldots,X_n} are non-negative discrete random variables. The function \code{\link{exact.q}} computes this probability by explicitly evaluating the sum, while \code{\link{sim.q}} estimates the probability using importance sampling. An exact approach that can handle relatively large problems is the use of \code{\link{dists.product.pair}} and \code{\link{dist.pair.cdf}}. All functions require the distributions of the \eqn{X_i}'s to be specified as a list with of two numeric vectors \code{x} and \code{fx}, where \code{x} denotes the outcomes and \code{fx} the probabilities. The outcomes should be non-negative. Only the first outcome can be 0 and only the last outcome can be Inf. The current function checks whether the distribution is properly specified.
#' 
#' @examples
#' # proper specification
#' dist <- list(x=c(1,2,3),fx=c(1/3,1/3,1/3))
#' check.dist(dist) # fine
#' 
#' # duplicates are not allowed
#' dist1 <- list(x=c(1,2,2),fx=c(0.1,0.2,0.7))
#' #check.dist(dist1) # would throw an error
#' # remove the duplicate
#' dist1 <- dist.unique.events(dist1)
#' check.dist(dist1) # fine
#' @export
check.dist <- function(dist){
  if (!isTRUE(!is.unsorted(dist$x,strictly = TRUE))){  stop("Not a proper dist: x should be strictly increasing and should not contain NAs")}
  if (isTRUE(is.infinite(dist$x[length(dist$x)-1]))){stop("Not a proper dist: only last element of dist$x may equal Inf")}
  if (!isTRUE(dist$x[1]>=0)){stop("Not a proper dist: elements of dist$x should be non-negative; dist$x[1] can be 0")}
  if (!isTRUE(all(dist$fx>=0))){stop("Not a proper dist: elements of dist$fx should be greather than or equal to 0")}
  if (!isTRUE(all(dist$fx<=1))){stop("Not a proper dist: elements of dist$fx should be smaller than or equal to 1")}
  if (!isTRUE(round(sum(dist$fx),digits = 6)-1==0)){stop("Not a proper dist: sum of dist$fx should be 1")  }  
  if (length(dist$x)!=length(dist$fx)) stop("Not a proper dist: lengths of dist$x and dist$fx must coincide")
  
  if (length(dist$x)==1L){
    if (!isTRUE((dist$x[1]>0)&(dist$x[1]<Inf))) stop("Not a proper dist: dist$x has length one so must be strictly positive and finite")
  }
  
  dist
}