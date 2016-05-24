#' Summary methods for ff objects
#'
#' @export
#' @export all.ff
#' @method all ff
#' @param x a \code{ff} object
#' @param ... optional other (\code{ff}) objects
#' @param na.rm should \code{NA} be removed?
#' @param range a \code{ri} or an \code{integer} vector of \code{length==2} giving a range restriction for chunked processing
#' @return \code{TRUE}, \code{FALSE} or \code{NA}
all.ff <- function(x,..., na.rm=FALSE, range=NULL){
  r <- checkRange(range,x)
  all( ...
     , sapply(chunk(x, from=min(r), to=max(r))
             , function(i){
                  all(x[i], na.rm=na.rm)
               }
             )
     )
}

#' Summary methods for ff objects
#' @export
#' @export any.ff
#' @method any ff
#' @param x a \code{ff} object
#' @param ... optional other (\code{ff}) objects
#' @param na.rm should \code{NA} be removed?
#' @param range a \code{ri} or an \code{integer} vector of \code{length==2} giving a range restriction for chunked processing
#' @return \code{TRUE}, \code{FALSE} or \code{NA}
any.ff <- function(x, ..., na.rm=FALSE, range=NULL){
  r <- checkRange(range,x)
  any( ...
     , sapply(chunk(x, from=min(r), to=max(r))
             , function(i){
                  any(x[i], na.rm=na.rm)
               }
             )
     )
}

#' \code{sum} returns the sum of all the values present in its arguments. 
#'
#' @title Sum of \code{ff} vector Elements
#' @method sum ff
#' @export
#' @export sum.ff
#' @param x a \code{ff} object
#' @param ... optional other (\code{ff}) objects
#' @param na.rm should \code{NA} be removed?
#' @param range a \code{ri} or an \code{integer} vector of \code{length==2} giving a range restriction for chunked processing
#' @return sum of elements
sum.ff <- function(x, ..., na.rm=FALSE, range=NULL){
   r <- checkRange(range, x)
   sum( ...
      , sapply( chunk(x, from=min(r), to=max(r))
              , function(i){
	               sum(x[i], na.rm=na.rm)
                }
	           )
	   )
}

#' Minimum, maximum and range of ff vector
#'
#' default behaviour of \code{\link{min}},\code{\link{max}} and \code{\link{range}}
#' @method min ff
#' @method max ff
#' @method range ff
#' @example ../examples/minmaxrange.R
#' @export
#' @export min.ff
#' @aliases min.ff max.ff range.ff
#' @param x a \code{ff} object
#' @param ... optional other (\code{ff}) objects
#' @param na.rm should \code{NA} be removed?
#' @param range a \code{ri} or an \code{integer} vector of \code{length==2} giving a range restriction for chunked processing
#' @return minimun, maximum or range values
min.ff <- function(x, ..., na.rm=FALSE, range=NULL){
   r <- checkRange(range, x)
   
	min( ...    #for all other ff's?
      , sapply( chunk(x, from=min(r), to=max(r))
	           , function(i){
                  #print(x[i])
                  min(x[i], na.rm=na.rm)
	             }
	           )
	   )
}

#' @export
#' @export max.ff
max.ff <- function(x, ..., na.rm=FALSE, range=NULL){
   r <- checkRange(range, x)
   
	max( ...    #for all other ff's?
      , sapply( chunk(x, from=min(r), to=max(r))
	           , function(i){
                  max(x[i], na.rm=na.rm)
	             }
	           )
	   )
}

#' @export
#' @export range.ff
range.ff <- function(x, ..., na.rm=FALSE, range=NULL){
  r <- checkRange(range, x)
	range( ...    #for all other ff's?
      , sapply( chunk(x, from=min(r), to=max(r))
	           , function(i){
                  #print(x[i])
                  range(x[i], na.rm=na.rm)
	             }
	           )
	   )
}
