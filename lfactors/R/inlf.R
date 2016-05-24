#' Implements \%in\% for lfactors
#' @param x same as \code{\link[base]{\%in\%}}
#' @param table same as \code{\link[base]{\%in\%}}
#' @seealso \code{\link[base]{\%in\%}}
#' @export
inlf <- function(x, table) {
  if(inherits(x, "lfactor")) {
  	#m1 <- match(x=as.character(x), table=as.character(table), nomatch=0)
  	#m2 <- match(x=as.character(switchllevels(x)), table=as.character(table), nomatch=0)
    m1 <- mlfactor(x=as.character(x), table=as.character(table), nomatch=0)
    m2 <- mlfactor(x=as.character(switchllevels(x)), table=as.character(table), nomatch=0)
  } else {
  	m1 <- mlfactor(x=as.character(x), table=as.character(table), nomatch=0)
  	m2 <- mlfactor(x=as.character(x), table=as.character(switchllevels(table)), nomatch=0)
  }
  ((m1 > 0) | (m2 > 0))
}

methods::setGeneric("%in%")
methods::setMethod("%in%", methods::signature(x="lfactor"), inlf)
methods::setMethod("%in%", methods::signature(table="lfactor"), inlf)
methods::setMethod("%in%", methods::signature(x="lfactor", table="lfactor"), inlf) 


