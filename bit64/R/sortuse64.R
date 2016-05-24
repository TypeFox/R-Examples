# /*
# R-Code for searching and merging
# S3 atomic 64bit integers for R
# (c) 2011 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2011-12-11
# Last changed:  2011-12-11
# */

#! \name{sortnut}
#! \alias{sortnut}
#! \alias{sortnut.integer64}
#! \alias{ordernut}
#! \alias{ordernut.integer64}
#! \alias{sortfin}
#! \alias{sortfin.integer64}
#! \alias{orderpos}
#! \alias{orderpos.integer64}
#! \alias{orderfin}
#! \alias{orderfin.integer64}
#! \alias{sortorderpos}
#! \alias{sortorderpos.integer64}
#! \alias{orderdup}
#! \alias{orderdup.integer64}
#! \alias{sortorderdup}
#! \alias{sortorderdup.integer64}
#! \alias{sortuni}
#! \alias{sortuni.integer64}
#! \alias{orderuni}
#! \alias{orderuni.integer64}
#! \alias{sortorderuni}
#! \alias{sortorderuni.integer64}
#! \alias{orderupo}
#! \alias{orderupo.integer64}
#! \alias{sortorderupo}
#! \alias{sortorderupo.integer64}
#! \alias{ordertie}
#! \alias{ordertie.integer64}
#! \alias{sortordertie}
#! \alias{sortordertie.integer64}
#! \alias{sorttab}
#! \alias{sorttab.integer64}
#! \alias{ordertab}
#! \alias{ordertab.integer64}
#! \alias{sortordertab}
#! \alias{sortordertab.integer64}
#! \alias{orderkey}
#! \alias{orderkey.integer64}
#! \alias{sortorderkey}
#! \alias{sortorderkey.integer64}
#! \alias{orderrnk}
#! \alias{orderrnk.integer64}
#! \alias{sortorderrnk}
#! \alias{sortorderrnk.integer64}
#! \alias{sortqtl}
#! \alias{sortqtl.integer64}
#! \alias{orderqtl}
#! \alias{orderqtl.integer64}
#! \title{
#!    Searching and other uses of sorting for 64bit integers
#! }
#! \description{
#!   This is roughly an implementation of hash functionality but based on sorting instead on a hasmap.
#!   Since sorting is more informative than hashingwe can do some more interesting things.
#! }
#! \usage{
#! sortnut(sorted, \dots)
#! ordernut(table, order, \dots)
#! sortfin(sorted, x, \dots)
#! orderfin(table, order, x, \dots)
#! orderpos(table, order, x, \dots)
#! sortorderpos(sorted, order, x, \dots)
#! orderdup(table, order, \dots)
#! sortorderdup(sorted, order, \dots)
#! sortuni(sorted, nunique, \dots)
#! orderuni(table, order, nunique, \dots)
#! sortorderuni(table, sorted, order, nunique, \dots)
#! orderupo(table, order, nunique, \dots)
#! sortorderupo(sorted, order, nunique, keep.order = FALSE, \dots)
#! ordertie(table, order, nties, \dots)
#! sortordertie(sorted, order, nties, \dots)
#! sorttab(sorted, nunique, \dots)
#! ordertab(table, order, nunique, \dots)
#! sortordertab(sorted, order, \dots)
#! orderkey(table, order, na.skip.num = 0L, \dots)
#! sortorderkey(sorted, order, na.skip.num = 0L, \dots)
#! orderrnk(table, order, na.count, \dots)
#! sortorderrnk(sorted, order, na.count, \dots)
#! \method{sortnut}{integer64}(sorted, \dots)
#! \method{ordernut}{integer64}(table, order, \dots)
#! \method{sortfin}{integer64}(sorted, x, method=NULL, \dots)
#! \method{orderfin}{integer64}(table, order, x, method=NULL, \dots)
#! \method{orderpos}{integer64}(table, order, x, nomatch=NA, method=NULL, \dots)
#! \method{sortorderpos}{integer64}(sorted, order, x, nomatch=NA, method=NULL, \dots)
#! \method{orderdup}{integer64}(table, order, method=NULL, \dots)
#! \method{sortorderdup}{integer64}(sorted, order, method=NULL, \dots)
#! \method{sortuni}{integer64}(sorted, nunique, \dots)
#! \method{orderuni}{integer64}(table, order, nunique, keep.order=FALSE, \dots)
#! \method{sortorderuni}{integer64}(table, sorted, order, nunique, \dots)
#! \method{orderupo}{integer64}(table, order, nunique, keep.order=FALSE, \dots)
#! \method{sortorderupo}{integer64}(sorted, order, nunique, keep.order = FALSE, \dots)
#! \method{ordertie}{integer64}(table, order, nties, \dots)
#! \method{sortordertie}{integer64}(sorted, order, nties, \dots)
#! \method{sorttab}{integer64}(sorted, nunique, \dots)
#! \method{ordertab}{integer64}(table, order, nunique, denormalize=FALSE, keep.order=FALSE, \dots)
#! \method{sortordertab}{integer64}(sorted, order, denormalize=FALSE, \dots)
#! \method{orderkey}{integer64}(table, order, na.skip.num = 0L, \dots)
#! \method{sortorderkey}{integer64}(sorted, order, na.skip.num = 0L, \dots)
#! \method{orderrnk}{integer64}(table, order, na.count, \dots)
#! \method{sortorderrnk}{integer64}(sorted, order, na.count, \dots)
#! \method{sortqtl}{integer64}(sorted, na.count, probs, \dots)
#! \method{orderqtl}{integer64}(table, order, na.count, probs, \dots)
#! }
#! \arguments{
#!   \item{x}{ an \code{\link{integer64}} vector }
#!   \item{sorted}{ a sorted \code{\link{integer64}} vector }
#!   \item{table}{ the original data with original order under the sorted vector }
#!   \item{order}{ an \code{\link{integer}} order vector that turns 'table' into 'sorted' }
#!   \item{nunique}{ number of unique elements, usually we get this from cache or call \code{sortnut} or \code{ordernut} }
#!   \item{nties}{ number of tied values, usually we get this from cache or call \code{sortnut} or \code{ordernut} }
#!   \item{denormalize}{ FALSE returns counts of unique values, TRUE returns each value with its counts }
#!   \item{nomatch}{ the value to be returned if an element is not found in the hashmap }
#!   \item{keep.order}{ determines order of results and speed: \code{FALSE} (the default) is faster and returns in sorted order, \code{TRUE} returns in the order of first appearance in the original data, but this requires extra work } 
#!   \item{probs}{ vector of probabilities in [0..1] for which we seek quantiles }
#!   \item{na.skip.num}{ 0 or the number of \code{NA}s. With 0, \code{NA}s are coded with 1L, with the number of \code{NA}s, these are coded with \code{NA}, the latter needed for \code{\link{as.factor.integer64}} }
#!   \item{na.count}{ the number of \code{NA}s, needed for this low-level function algorithm }
#!   \item{method}{ see details }
#!   \item{\dots}{ further arguments, passed from generics, ignored in methods }
#! }
#! \details{
#! \tabular{rrrrl}{
#!    \bold{sortfun} \tab \bold{orderfun} \tab \bold{sortorderfun} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{sortnut} \tab \code{ordernut} \tab                     \tab  \tab return number of tied and of unique values \cr
#!    \code{sortfin} \tab \code{orderfin} \tab                     \tab \code{\link{\%in\%.integer64}} \tab return logical whether \code{x} is in \code{table} \cr
#!                   \tab \code{orderpos} \tab \code{sortorderpos} \tab \code{\link[=match.integer64]{match}} \tab return positions of \code{x} in \code{table} \cr
#!                   \tab \code{orderdup} \tab \code{sortorderdup} \tab \code{\link[=duplicated.integer64]{duplicated}} \tab return logical whether values are duplicated \cr
#!    \code{sortuni} \tab \code{orderuni} \tab \code{sortorderuni} \tab \code{\link[=unique.integer64]{unique}} \tab return unique values (=dimensiontable) \cr
#!                   \tab \code{orderupo} \tab \code{sortorderupo} \tab \code{\link[=unique.integer64]{unique}} \tab return positions of unique values \cr
#!                   \tab \code{ordertie} \tab \code{sortordertie} \tab  \tab return positions of tied values \cr
#!                   \tab \code{orderkey} \tab \code{sortorderkey} \tab  \tab positions of values in vector of unique values (match in dimensiontable) \cr
#!    \code{sorttab} \tab \code{ordertab} \tab \code{sortordertab} \tab \code{\link[=table.integer64]{table}} \tab tabulate frequency of values  \cr
#!                   \tab \code{orderrnk} \tab \code{sortorderrnk} \tab  \tab rank averaging ties \cr
#!    \code{sortqtl} \tab \code{orderqtl} \tab                     \tab  \tab return quantiles given probabilities \cr
#! }
#! The functions \code{sortfin}, \code{orderfin}, \code{orderpos} and \code{sortorderpos} each offer three algorithms for finding \code{x} in \code{table}.  \cr
#! With \code{method=1L} each value of \code{x} is searched independently using \emph{binary search}, this is fastest for small \code{table}s. \cr
#! With \code{method=2L} the values of \code{x} are first sorted and then searched using \emph{doubly exponential search}, this is the best allround method. \cr
#! With \code{method=3L} the values of \code{x} are first sorted and then searched using simple merging, this is the fastest method if \code{table} is huge and \code{x} has similar size and distribution of values. \cr
#! With \code{method=NULL} the functions use a heuristic to determine the fastest algorithm. \cr
#!
#! The functions \code{orderdup} and \code{sortorderdup} each offer two algorithms for setting the truth values in the return vector.  \cr
#! With \code{method=1L} the return values are set directly which causes random write access on a possibly large return vector. \cr
#! With \code{method=2L} the return values are first set in a smaller bit-vector -- random access limited to a smaller memory region -- and finally written sequentially to the logical output  vector. \cr
#! With \code{method=NULL} the functions use a heuristic to determine the fastest algorithm. \cr
#! }
#! \value{
#!   see details
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ programming }
#! \keyword{ manip }
#! \seealso{ \code{\link[=match.integer64]{match}} }
#! \examples{
#!  message("check the code of 'optimizer64' for examples:")
#!  print(optimizer64)
#! }



sortnut <- function(sorted, ...)UseMethod("sortnut")
sortnut.integer64 <- function(sorted, ...)
{
  ret <- .Call("r_ram_integer64_sortnut", x = sorted, PACKAGE = "bit64")
  names(ret) <- c("nunique","nties")
  ret
}

ordernut <- function(table, order, ...)UseMethod("ordernut")
ordernut.integer64 <- function(table, order, ...)
{
  ret <- .Call("r_ram_integer64_ordernut", table = as.integer64(table), order = as.integer(order), PACKAGE = "bit64")
  names(ret) <- c("nunique","nties")
  ret
}

sortfin <- function(sorted, x, ...)UseMethod("sortfin")
sortfin.integer64 <- function(sorted, x, method=NULL, ...)
{
  n <- length(x)
  if (is.null(method)){
	if (n<2048){
	  method <- 1L
	}else if (n<length(sorted)/128){
	  method <- 2L
	}else{
	  method <- 3L
	}
  }else method <- as.integer(method)
  ret <- logical(n)
  if (method==1L){
	  .Call("r_ram_integer64_sortfin_asc"
	  , x = as.integer64(x)
	  , sorted = as.integer64(sorted)
	  , method= method
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
  }else{
    sx <- clone(as.integer64(x)); o <- seq_along(x); ramsortorder(sx, o, na.last=FALSE, ...)
	ret[o] <- .Call("r_ram_integer64_sortfin_asc"
	  , x = sx
	  , sorted = as.integer64(sorted)
	  , method= method
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
	  ret
  }
}

orderfin <- function(table, order, x, ...)UseMethod("orderfin")
orderfin.integer64 <- function(table, order, x, method=NULL, ...)
{
  n <- length(x)
  if (is.null(method)){
	if (n<4096){
	  method <- 1L
	}else if (n<length(table)/8){
	  method <- 2L
	}else{
	  method <- 3L
	}
  }else method <- as.integer(method)
  ret <- logical(n)
  if (method==1L){
	  .Call("r_ram_integer64_orderfin_asc"
	  , x = as.integer64(x)
	  , table = as.integer64(table)
	  , order = as.integer(order)
	  , method= as.integer(method)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
  }else{
    o <- seq_along(x); ramorder(x, o, na.last=FALSE, ...)
	ret[o] <- .Call("r_ram_integer64_orderfin_asc"
	  , x = x[o]
	  , table = as.integer64(table)
	  , order = as.integer(order)
	  , method= as.integer(method)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
	  ret
  }
}


orderpos <- function(table, order, x, ...)UseMethod("orderpos")
orderpos.integer64 <- function(table, order, x, nomatch=NA, method=NULL, ...)
{
  n <- length(x)
  if (is.null(method)){
	if (n<4096){
	  method <- 1L
	}else if (n<length(table)/8){
	  method <- 2L
	}else{
	  method <- 3L
	}
  }else method <- as.integer(method)
  ret <- integer(n);
  if (method==1L){
	  .Call("r_ram_integer64_orderpos_asc"
	  , x = as.integer64(x)
	  , table = as.integer64(table)
	  , order = as.integer(order)
	  , nomatch = as.integer(nomatch)
	  , method= as.integer(method)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
  }else{
    o <- seq_along(x); ramorder(x, o, na.last=FALSE, ...)
	ret[o] <- .Call("r_ram_integer64_orderpos_asc"
	  , x = x[o]
	  , table = as.integer64(table)
	  , order = as.integer(order)
	  , nomatch = as.integer(nomatch)
	  , method= as.integer(method)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
	  ret
  }
}

sortorderpos <- function(sorted, order, x, ...)UseMethod("sortorderpos")
sortorderpos.integer64 <- function(sorted, order, x, nomatch=NA, method=NULL, ...)
{
  n <- length(x)
  if (is.null(method)){
	if (n<2048){
	  method <- 1L
	}else if (n<length(sorted)/128){
	  method <- 2L
	}else{
	  method <- 3L
	}
  }else method <- as.integer(method)
  ret <- integer(n)
  if (method==1L){
	  .Call("r_ram_integer64_sortorderpos_asc"
	  , x = as.integer64(x)
	  , sorted = as.integer64(sorted)
	  , order = as.integer(order)
	  , nomatch = as.integer(nomatch)
	  , method= as.integer(method)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
  }else{
    sx <- clone(as.integer64(x)); o <- seq_along(x); ramsortorder(sx, o, na.last=FALSE, ...)
	ret[o] <- .Call("r_ram_integer64_sortorderpos_asc"
	  , x = sx
	  , sorted = as.integer64(sorted)
	  , order = as.integer(order)
	  , nomatch = as.integer(nomatch)
	  , method= as.integer(method)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
	  ret
  }
}



orderdup <- function(table, order, ...)UseMethod("orderdup")
orderdup.integer64 <- function(table, order, method=NULL, ...)
{
  if (is.null(method)){
    if (length(table)<4194304)
	  method <- 1L
	else
	  method <- 2L
  }else method <- as.integer(method)
  ret <- logical(length(table))
  .Call("r_ram_integer64_orderdup_asc"
  , table = as.integer64(table)
  , order = as.integer(order)
  , method = method
  , ret = ret
  , PACKAGE = "bit64"
  )
}


sortorderdup <- function(sorted, order, ...)UseMethod("sortorderdup")
sortorderdup.integer64 <- function(sorted, order, method=NULL, ...)
{
  if (is.null(method)){
    if (length(sorted)<4194304)
	  method <- 1L
	else
	  method <- 2L
  }else method <- as.integer(method)
  ret <- logical(length(sorted))
  .Call("r_ram_integer64_sortorderdup_asc"
  , sorted = as.integer64(sorted)
  , order = as.integer(order)
  , method = method
  , ret = ret
  , PACKAGE = "bit64"
  )
}




sortuni <- function(sorted, nunique, ...)UseMethod("sortuni")
sortuni.integer64 <- function(sorted, nunique, ...)
{
  ret <- integer64(nunique)
  .Call("r_ram_integer64_sortuni_asc"
  , sorted = as.integer64(sorted)
  , ret = ret
  , PACKAGE = "bit64"
  )
}

orderuni <- function(table, order, nunique, ...)UseMethod("orderuni")
orderuni.integer64 <- function(table, order, nunique, keep.order=FALSE, ...)
{
  ret <- integer64(nunique)
  .Call("r_ram_integer64_orderuni_asc"
  , table = as.integer64(table)
  , order = as.integer(order)
  , keep.order = as.logical(keep.order)
  , ret = ret
  , PACKAGE = "bit64"
  )
}

sortorderuni <- function(table, sorted, order, nunique, ...)UseMethod("sortorderuni")
sortorderuni.integer64 <- function(table, sorted, order, nunique, ...)
{
  ret <- integer64(nunique)
	  .Call("r_ram_integer64_sortorderuni_asc"
	  , table = as.integer64(table)
	  , sorted = as.integer64(sorted)
	  , order = as.integer(order)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
}

orderupo <- function(table, order, nunique, ...)UseMethod("orderupo")
orderupo.integer64 <- function(table, order, nunique, keep.order=FALSE, ...)
{
	ret <- integer(nunique)
	.Call("r_ram_integer64_orderupo_asc"
	, table = as.integer64(table)
	, order = as.integer(order)
	, keep.order = as.logical(keep.order)
	, ret = ret
	, PACKAGE = "bit64"
	)
}

sortorderupo <- function(sorted, order, nunique, keep.order=FALSE, ...)UseMethod("sortorderupo")
sortorderupo.integer64 <- function(sorted, order, nunique, keep.order=FALSE, ...)
{
	ret <- integer(nunique)
	ret2 <- .Call("r_ram_integer64_sortorderupo_asc"
	, sorted = as.integer64(sorted)
	, order = as.integer(order)
	, keep.order = as.logical(keep.order)
	, ret = ret
	, PACKAGE = "bit64"
	)
	ret2
}


ordertie <- function(table, order, nties, ...)UseMethod("ordertie")
ordertie.integer64 <- function(table, order, nties, ...)
{
  ret <- integer(nties)
  .Call("r_ram_integer64_ordertie_asc"
  , table = as.integer64(table)
  , order = as.integer(order)
  , ret = ret
  , PACKAGE = "bit64"
  )
}

sortordertie <- function(sorted, order, nties, ...)UseMethod("sortordertie")
sortordertie.integer64 <- function(sorted, order, nties, ...)
{
  ret <- integer(nties)
	  .Call("r_ram_integer64_sortordertie_asc"
	  , sorted = as.integer64(sorted)
	  , order = as.integer(order)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
}


sorttab <- function(sorted, nunique, ...)UseMethod("sorttab")
sorttab.integer64 <- function(sorted, nunique, ...)
{
  ret <- integer(nunique)
  .Call("r_ram_integer64_sorttab_asc"
  , sorted = as.integer64(sorted)
  , ret = ret
  , PACKAGE = "bit64"
  )
}

ordertab <- function(table, order, nunique, ...)UseMethod("ordertab")
ordertab.integer64 <- function(table, order, nunique, denormalize=FALSE, keep.order=FALSE, ...)
{
  denormalize <- as.logical(denormalize)
  keep.order <- as.logical(keep.order)
  ret <- integer(if (denormalize || keep.order) length(table) else nunique) 
  .Call("r_ram_integer64_ordertab_asc"
  , table = as.integer64(table)
  , order = as.integer(order)
  , denormalize = denormalize
  , keep.order = keep.order
  , ret = ret
  , PACKAGE = "bit64"
  )
}

sortordertab <- function(sorted, order, ...)UseMethod("sortordertab")
sortordertab.integer64 <- function(sorted, order, denormalize=FALSE, ...)
{
  ret <- integer(length(sorted))
	  .Call("r_ram_integer64_sortordertab_asc"
	  , sorted = as.integer64(sorted)
	  , order = as.integer(order)
	  , denormalize = as.logical(denormalize)
	  , ret = ret
	  , PACKAGE = "bit64"
	  )
}

orderkey <- function(table, order, na.skip.num=0L, ...)UseMethod("orderkey")
orderkey.integer64 <- function(table, order, na.skip.num=0L, ...)
{
	ret <- integer(length(table))
	.Call("r_ram_integer64_orderkey_asc"
	, table = as.integer64(table)
	, order = as.integer(order)
	, na.skip.num=na.skip.num
	, ret = ret
	, PACKAGE = "bit64"
	)
}

sortorderkey <- function(sorted, order, na.skip.num=0L, ...)UseMethod("sortorderkey")
sortorderkey.integer64 <- function(sorted, order, na.skip.num=0L, ...)
{
	ret <- integer(length(sorted))
	.Call("r_ram_integer64_sortorderkey_asc"
	, sorted = as.integer64(sorted)
	, order = as.integer(order)
	, na.skip.num=na.skip.num
	, ret = ret
	, PACKAGE = "bit64"
	)
}


orderrnk <- function(table, order, na.count, ...)UseMethod("orderrnk")
orderrnk.integer64 <- function(table, order, na.count, ...)
{
  ret <- double(length(table))
  .Call("r_ram_integer64_orderrnk_asc"
  , table = as.integer64(table)
  , order = as.integer(order)
  , na.count=as.integer(na.count)
  , ret = ret
  , PACKAGE = "bit64"
  )
}

sortorderrnk <- function(sorted, order, na.count, ...)UseMethod("sortorderrnk")
sortorderrnk.integer64 <- function(sorted, order, na.count, ...)
{
  ret <- double(length(sorted))
  .Call("r_ram_integer64_sortorderrnk_asc"
  , sorted = as.integer64(sorted)
  , order = as.integer(order)
  , na.count=as.integer(na.count)
  , ret = ret
  , PACKAGE = "bit64"
  )
}


sortqtl <- function(sorted, na.count, probs, ...)UseMethod("sortqtl")
sortqtl.integer64 <- function(sorted, na.count, probs, ...){
	n <- length(sorted) - na.count  # nvalid
	ret <- sorted[na.count + round(1L + probs*(n-1L))]
	ret[is.na(probs)] <- NA ## xx this fix only neccessary until we have C-implementation of [.integer64 handling NA
	ret
}

orderqtl <- function(table, order, na.count, probs, ...)UseMethod("orderqtl")
orderqtl.integer64 <- function(table, order, na.count, probs, ...){
	n <- length(table) - na.count  # nvalid
	ret <- table[ order[na.count + round(1L + probs*(n-1L))] ]
	ret[is.na(probs)] <- NA ## xx this fix only neccessary until we have C-implementation of [.integer64 handling NA
	ret
}
