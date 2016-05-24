##' This function select objects by their position specified by index.
##'
##' This function select object by their position specified by index. When index is vector, a list of selected objects is returned. Main use of this function is to select time-specific covariate with \code{getIndex} function.
##' @title Select objects by their position.
##' @param index a numeric value or vector specifying position.
##' @param ... an arbitrary number of objects to be selected.
##' @return a selected object or a list of selected objects.
##' @seealso \code{\link{getIndex}}, \code{\link{epifit}}
##' @examples
##' id <- 1:10
##' cov1 <- rnorm(10)
##' cov2 <- rbinom(10, 4, 0.5)
##' cov3 <- runif(10)
##' covset1 <- selectByIndex(c(1,2), id, cov1, cov2, cov3)
##' cov <- selectByIndex(3, id, cov1, cov2, cov3)
##' @export
selectByIndex <- function(index, ...){
  args <- list(...)
  if(length(index) == 1){
    if(length(args) < index || index < 1)
      stop("invalid index in select function")
    return(args[[index]])
  } else {
    if(max(index) > length(args) || min(index) < 1)
      stop("invalid index in select function")
    lapply(index, function(x, lst){lst[[x]]}, lst=args)
  }
}
