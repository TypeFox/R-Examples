#' @title Combine rfPermute Objects
#' @description Combines two or more ensembles of \code{rfPermute} objects into 
#'   one, combining \code{randomForest} results, null distributions, 
#'   and re-calculating p-values.
#' 
#' @param \dots two or more objects of class \code{rfPermute}, to be combined 
#'   into one.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link[randomForest]{combine}}
#' 
#' @examples
#' 
#' data(iris)
#' rp1 <- rfPermute(Species ~ ., iris, ntree = 50, norm.votes = FALSE, nrep = 100)
#' rp2 <- rfPermute(Species ~ ., iris, ntree = 50, norm.votes = FALSE, nrep = 100)
#' rp3 <- rfPermute(Species ~ ., iris, ntree = 50, norm.votes = FALSE, nrep = 100)
#' rp.all <- rp.combine(rp1, rp2, rp3)
#' plot(rp.all)
#' 
#' @importFrom abind abind
#' @importFrom randomForest combine
#' @export
#' 
rp.combine <- function(...) {
  rp.list <- list(...)
  are.rp <- sapply(rp.list, function(x) inherits(x, "rfPermute"))
  if(any(!are.rp)) stop("some objects in '...' are not rfPermute results.")
  
  rf <- do.call(combine, rp.list)
  rf$null.dist <- sapply(c("unscaled", "scaled"), function(sc) {
   do.call(abind, c(lapply(rp.list, function(x) x$null.dist[[sc]]), along = 3))
  }, simplify = FALSE)
  rf$pval <- .calcImpPval(rf)
  
  class(rf) <- c("rfPermute", "randomForest")
  return(rf)  
}