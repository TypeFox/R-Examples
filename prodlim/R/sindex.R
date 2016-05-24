#' Index for evaluation of step functions.
#' 
#' Returns an index of positions. Intended for evaluating a step function at
#' selected times. The function counts how many elements of a vector, e.g. the
#' jump times of the step function, are smaller or equal to the elements in a
#' second vector, e.g. the times where the step function should be evaluated.
#' 
#' If all \code{jump.times} are greater than a particular \code{eval.time} the
#' sindex returns \code{0}. This must be considered when sindex is used for
#' subsetting, see the Examples below.
#' 
#' @param jump.times Numeric vector: e.g.\ the unique jump times of a step
#' function.
#' @param eval.times Numeric vector: e.g.\ the times where the step function
#' should be evaluated
#' @param strict If TRUE make the comparison of jump times and eval times
#' strict
#' @param comp If "greater" count the number of jump times that are greater
#' (greater or equal when strict==FALSE) than the eval times
#' @return Index of the same length as \code{eval.times} containing the numbers
#' of the \code{jump.times} that are smaller than or equal to
#' \code{eval.times}.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @keywords misc
#' @examples
#' 
#' 
#' test <- list(time = c(1, 1,5,5,2,7,9),
#' 		status = c(1,0,1,0,1,1,0))
#' fit <- prodlim(Hist(time,status)~1,data=test)
#' jtimes <- fit$time
#' etimes <- c(0,.5,2,8,10)
#' fit$surv
#' c(1,fit$surv)[1+sindex(jtimes,etimes)]
#' 
#' @export
"sindex" <- function(jump.times,eval.times,comp="smaller",strict=FALSE) {
  stopifnot(is.numeric(jump.times))
  stopifnot(is.numeric(eval.times))
  N <- length(jump.times)
  if (comp=="greater"){
    N-sindex(jump.times=jump.times,
             eval.times=eval.times,
             comp="smaller",
             strict=!strict)
  }
  else{
    neval <- length(eval.times)
    if (!(neval> 0 && N >0)) stop("missing data")
    new.order <- order(eval.times)
    ind <- .C("sindex",index = integer(neval),as.double(sort(jump.times)),as.double(eval.times[new.order]),as.integer(N),as.integer(neval),as.integer(strict),PACKAGE="prodlim")$index
    ind[order(new.order)]
  }
}

## sindexStrata <- function(jump.times,
                         ## first,
                         ## size,
                         ## eval.times,
                         ## strict=FALSE) {
  ## stopifnot(is.numeric(jump.times))
  ## stopifnot(is.numeric(eval.times))
  ## NK <- length(size)
  ## stopifnot(length(first)==NK)
  ## N <- length(jump.times)
  ## neval <- length(eval.times)
  ## if (!(neval> 0 && N >0)) stop("missing data")
  ## new.order <- order(eval.times)
  ## ind <- .C("sindexStrata",
            ## index = integer(neval),
            ## as.double(sort(jump.times)),
            ## as.double(eval.times[new.order]),
            ## as.integer(N),
            ## as.integer(neval),
            ## as.integer(strict),
            ## DUP=FALSE,
            ## PACKAGE="prodlim")$index
  ## ind[order(new.order)]
## }

