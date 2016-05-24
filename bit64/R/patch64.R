# /*
# R-Code for patching S3 generics
# S3 atomic 64bit integers for R
# (c) 2011 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2011-12-11
# Last changed:  2011-12-11
# */

#! \name{bit64S3}
#! \alias{bit64S3}
#! \alias{:}
#! \alias{:.default}
#! \alias{:.integer64}
#! \alias{is.double}
#! \alias{is.double.default}
#! \alias{is.double.integer64}
#! \alias{match}
#! \alias{match.default}
#! \alias{\%in\%}
#! \alias{\%in\%.default}
#! \alias{rank}
#! \alias{rank.default}
#! %\alias{table}
#! %\alias{table.default}
#! \alias{order}
#! \alias{order.default}
#! \title{
#!   Tunrning base R functions into S3 generics for bit64 
#! }
#! \description{
#! 	Turn those base functions S3 generic which are used in bit64
#! }
#! \usage{
#!	from:to
#!  #--as-cran complains about \method{:}{default}(from, to)
#!  #--as-cran complains about \method{:}{integer64}(from, to)
#!	is.double(x)
#!  \method{is.double}{default}(x)
#!  \method{is.double}{integer64}(x)
#! 	match(x, table, ...)
#!  \method{match}{default}(x, table, ...)
#! 	x \%in\% table
#!  \method{\%in\%}{default}(x, table)
#! 	rank(x, ...)
#!  \method{rank}{default}(x, ...)
#! 	%table(...)
#!  %\method{table}{default}(...)
#! 	order(...)
#!  \method{order}{default}(...)
#! }
#! \arguments{
#!   \item{x}{
#! 	integer64 vector: the values to be matched, optionally carrying a cache created with \code{\link{hashcache}}
#! }
#!   \item{table}{
#! 	integer64 vector: the values to be matched against, optionally carrying a cache created with \code{\link{hashcache}} or \code{\link{sortordercache}}
#! }
#!   \item{from}{ scalar denoting first element of sequence }
#!   \item{to}{ scalar denoting last element of sequence }
#!   \item{\dots}{ ignored }
#! }
#! \details{
#!    The following functions are turned into S3 gernerics in order to dispatch methods for \code{\link{integer64}}: 
#!    \preformatted{
#! 	   \code{\link{:}}
#! 	   \code{\link{is.double}}
#! 	   \code{\link{match}}
#! 	   \code{\link{\%in\%}}
#! 	   %\code{\link{table}}
#! 	   \code{\link{rank}}
#! 	   \code{\link{order}}
#!    }
#! }
#! \value{
#! 	\code{\link{invisible}}
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \note{
#! 	\code{\link{is.double}} returns \code{FALSE} for \code{\link{integer64}} \cr
#! 	\code{\link{:}} currently only dispatches at its first argument, thus \code{as.integer64(1):9} works but \code{1:as.integer64(9)} doesn't
#! 	\code{\link{match}} currently only dispatches at its first argument and expects its second argument also to be integer64, otherwise throws an error. Beware of something like \code{match(2, as.integer64(0:3))}
#! 	\code{\link{\%in\%}} currently only dispatches at its first argument and expects its second argument also to be integer64, otherwise throws an error. Beware of something like \code{2 \%in\% as.integer64(0:3)}
#! 	\code{\link{order}} currently only orders a single argument, trying more than one raises an error
#! }
#! \seealso{
#! 	\code{\link{bit64}}, \code{\link{S3}}
#! }
#! \examples{
#!  is.double(as.integer64(1))
#! 	as.integer64(1):9
#!  match(as.integer64(2), as.integer64(0:3))
#!  as.integer64(2) \%in\% as.integer64(0:3)
#!  
#!  unique(as.integer64(c(1,1,2)))
#!  rank(as.integer64(c(1,1,2)))
#!  
#!  %table(as.integer64(c(1,1,2)))
#!  %table(as.integer64(c(1,1,2)),as.integer64(c(3,4,4)))
#!  %table(as.integer64(c(1,1,2)),c(3,4,4))
#!  %table(c(1,1,2),as.integer64(c(3,4,4)))
#!  
#!  order(as.integer64(c(1,NA,2)))
#!  
#!  \dontshow{
#!  stopifnot(identical(match(as.integer64(2), as.integer64(0:3)), match(2, 0:3)))
#!  stopifnot(identical(as.integer64(2) \%in\% as.integer64(0:3), 2 \%in\% 0:3))
#!  
#!  stopifnot(identical(unique(as.integer64(c(1,1,2))), as.integer64(unique(c(1,1,2)))))
#!  stopifnot(identical(rank(as.integer64(c(1,1,2))), rank(c(1,1,2))))
#!  
#!  %stopifnot(identical(table(as.integer64(c(1,1,2))), table(c(1,1,2))))
#!  %stopifnot(identical(table(as.integer64(c(1,1,2)),as.integer64(c(3,4,4))), table(c(1,1,2),c(3,4,4))))
#!  %stopifnot(identical(table(as.integer64(c(1,1,2)),c(3,4,4)), table(c(1,1,2),c(3,4,4))))
#!  %stopifnot(identical(table(c(1,1,2),as.integer64(c(3,4,4))), table(c(1,1,2),c(3,4,4))))
#!  
#!  stopifnot(identical(order(as.integer64(c(1,NA,2))), order(c(1,NA,2))))
#!  stopifnot(identical(order(as.integer64(c(1,NA,2)), decreasing=TRUE), order(c(1,NA,2), decreasing=TRUE)))
#!  }
#! }
#! \keyword{ methods }

# OCT 2013: bit64S3() at wish of CRAN maintainers replaced by direct conversion to S3 generics
# in order to avoid assigning to globalenv

if (!exists(":.default")){
	":" <- function(from,to) UseMethod(":")
	":.default" <- function(from,to) base::":"(from,to)
}
":.integer64" <- function(from, to)seq.integer64(from=from, to=to)

if (!exists("is.double.default")){
	"is.double" <- function(x) UseMethod("is.double")
	"is.double.default" <- function(x) base::"is.double"(x)
}
"is.double.integer64" <- function(x)FALSE

if (!exists("match.default")){
	"match" <- function(x, table, ...) UseMethod("match")
	"match.default" <- function(x, table, ...) base::"match"(x, table, ...)
}

if (!exists("%in%.default")){
	"%in%" <- function(x, table) UseMethod("%in%")
	"%in%.default" <- function(x, table) base::"%in%"(x, table)
}

if (!exists("rank.default")){
	"rank" <- function(x, ...) UseMethod("rank")
	"rank.default" <- function(x, ...) base::"rank"(x, ...)
}

# not yet able to combinewith other column types - better leave table() as is and hope for as.factor.integer64
#if (!exists("table.default")){
#	"table" <- function(...) UseMethod("table")
#	"table.default" <- function(...) base::"table"(...)
#}

if (!exists("order.default")){
	"order" <- function(...) UseMethod("order")
	"order.default" <- function(...) base::"order"(...)
}

