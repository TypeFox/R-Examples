`%nin%` <- function (x, table) is.na(match(x, table))

## TODO 1) Fast version of  !all(x %in% tab) == any(x %nin% tab)
##      2) Fast version of  !any(x %in% tab) == all(x %nin% tab)
##     3a) Fast version of  !(x %in% A  |  x %in  B)
##     3b) Fast version of  !(x %in% A ||  x %in  B) = 3a) where length(x)==1

##' Compute \eqn{\log( (1+r)/(1-r) )} in a numerically stable way.  Note that
##' \eqn{(1+r)/(1-r) = 1 + 2r/(1-r)};  hence the logarithm is
##' \eqn{log(*) = log(1 + 2r/(1-r)) = log1p(2r/(1-r))}.
##'
##' @title  log( (1+r)/(1-r) )  in a Numerically Stable Way
##' @param r numeric (vector) in \eqn{[0,1]}.
##' @return numeric vector \dQuote{as} \code{r}.
##' @author Martin Maechler
log.q1pm <- function(r) log1p(2*r/(1-r))



check.Rgraphviz <- function() {
  if(!requireNamespace("Rgraphviz"))
    stop("Package 'Rgraphviz' (from Bioconductor) must be installed for plotting graphs!")
}

