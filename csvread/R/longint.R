#-------------------------------------------------------------------------------
#
# Package csvread 
#
# int64 class - not all numeric functionality is implemented.  
# 
# Sergei Izrailev, 2011-2014
#-------------------------------------------------------------------------------
# Copyright 2011-2014 Collective, Inc.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-------------------------------------------------------------------------------

#' A very basic 64-bit integer class. 
#' 
#' The \code{int64} class stores 64-bit integers in vectors of doubles and the 
#' base as an attribute \code{base} of the vector for printing and conversion to
#' character. The motivation behind this class is to give R the ability to load
#' 64-bit integers directly, for example, to represent the commonly used 64-bit
#' identifiers in relational and other databases.
#' @name int64
#' @title A very basic 64-bit integer class. 
#' @aliases int64 as.int64 as.int64.default as.int64.factor as.int64.character 
#'          as.int64.numeric as.int64.NULL [.int64 [[.int64 [<-.int64
#' @param x Object to be coerced or tested
#' @param length A non-negative integer specifying the desired length.  
#'          Double values will be coerced to integer: supplying an argument of 
#'          length other than one is an error.
#' @param ... Further arguments passed to or from other methods.
#' @seealso Ops.int64 csvread
int64 <- function(length = 0)
{
   res <- double(length)
   class(res) <- "int64"
   return(res)
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
is.int64 <- function(x) inherits(x, "int64")

#-------------------------------------------------------------------------------

as.int64 <- function(x, ...) UseMethod("as.int64")

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.int64 default
as.int64.default <- function(x, ...)
{
	if(inherits(x, "int64")) return(x)
	stop(gettextf("do not know how to convert '%s' to class \"int64\"",
					deparse(substitute(x))))
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.int64 factor
as.int64.factor <- function(x, ...) as.int64(as.character(x), ...)

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.int64 character
as.int64.character <- function(x, base=10L, ...)
{
   base <- as.integer(base)
   if (base > 36 || base < 2) stop("Can't convert to int64: invalid base.")
	return(.Call("charToInt64", x, base, PACKAGE="csvread"))
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.int64 numeric
as.int64.numeric <- function(x, ...)
{
   return(.Call("doubleToInt64", as.double(x), PACKAGE="csvread"));
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.int64 NULL
as.int64.NULL <- function(x, ...) 
{
   res <- double()
   class(res) <- "int64"
   return(res)
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method format int64
format.int64 <- function(x, ...)
{
   format(as.character(x), ...)
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method print int64
print.int64 <- function(x, ...)
{
	print(format(x), quote = FALSE, ...)
	invisible(x)
}

#-------------------------------------------------------------------------------

#' Operators for the \code{int64} class.
#' 
#' Operators for the \code{int64} class: one of 
#' \code{+}, \code{-}, \code{==}, \code{!=}, \code{<}, \code{<=}, \code{>} or \code{>=}.
#' @rdname Ops.int64
#' @aliases + - <
#' @param e1 int64 object, character vector or numeric vector 
#'        (character and numeric values are converted by \code{as.int64}).
#' @param e2 int64 object, character vector or numeric vector 
#'        (character and numeric values are converted by \code{as.int64}).
#' @usage e1 + e2
#' e1 - e2
#' @seealso int64
#' @export
#' @method Ops int64
Ops.int64 <- function(e1, e2)
{
	if (nargs() == 1)
		stop("unary ", .Generic, " not defined for int64 objects")
	boolean <- switch(.Generic, "<" =, ">" =, "==" =,
			"!=" =, "<=" =, ">=" = TRUE,
			FALSE)
	if (!boolean) stop(.Generic, " not defined for int64 objects")
	## allow character and numeric args to be coerced to int64
   if (is.character(e1)) e1 <- as.int64(e1)
   else if (is.numeric(e1)) e1 <- as.int64(e1)
   if (is.character(e2)) e2 <- as.int64(e2)
   else if (is.numeric(e2)) e2 <- as.int64(e2)
   NextMethod(.Generic)
}

#-------------------------------------------------------------------------------

#' @rdname Ops.int64
#' @export
#' @method + int64
`+.int64` <- function(e1, e2)
{	
	if (nargs() == 1) return(e1)
	if (inherits(e1, "int64") && inherits(e2, "int64"))
   {
      return(.Call("addInt64Int64", e1, e2, PACKAGE="csvread"))
   }
	if (inherits(e1, "int64")) return(.Call("addInt64Int", e1, as.integer(e2), PACKAGE="csvread"))
   if (inherits(e2, "int64")) return(.Call("addInt64Int", e2, as.integer(e1), PACKAGE="csvread"))
}

#-------------------------------------------------------------------------------

#' @rdname Ops.int64
#' @export
#' @method - int64
`-.int64` <- function(e1, e2)
{
   if (nargs() == 1) return(-e1)
   if(inherits(e1, "int64") && inherits(e2, "int64"))
   {
      return(.Call("subInt64Int64", e1, e2, PACKAGE="csvread"))
   }
   if (inherits(e1, "int64")) return(.Call("subInt64Int64", e1, as.int64(e2), PACKAGE="csvread"))
   if (inherits(e2, "int64")) return(.Call("subInt64Int64", as.int64(e1), e2, PACKAGE="csvread"))
}

#-------------------------------------------------------------------------------

`[.int64` <- function(x, ..., drop = TRUE)
{
	cl <- oldClass(x)
	class(x) <- NULL
	val <- NextMethod("[")
	class(val) <- cl
	val
}

#-------------------------------------------------------------------------------

`[[.int64` <- function(x, ..., drop = TRUE)
{
	cl <- oldClass(x)
	class(x) <- NULL
	val <- NextMethod("[[")
	class(val) <- cl
	val
}

#-------------------------------------------------------------------------------

`[<-.int64` <- function(x, ..., value)
{
	if (!length(value)) return(x)
	value <- unclass(as.int64(value))
	cl <- oldClass(x)
	class(x) <- NULL
	x <- NextMethod(.Generic)
	class(x) <- cl
	x
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @param base Specifies the base of the number (default is the base attribute of the object).
#' @method as.character int64
as.character.int64 <- function(x, base = NULL, ...)
{
   if (is.null(base))
   {
      # try attributes first
      i <- match("base", names(attributes(x)))
      if (is.na(i)) base <- 10L
      else base <- as.integer(attr(x, "base"))
   }
   if (base == 10L) return(.Call("int64ToChar", x, PACKAGE="csvread"))
   if (base != 16L) stop(paste("Can't convert int64 to character as base", base))
   return(.Call("int64ToHex", x, PACKAGE="csvread"))
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.double int64
as.double.int64 <- function(x, ...)
{
   return(.Call("int64ToDouble", x, PACKAGE="csvread"))
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.integer int64
as.integer.int64 <- function(x, ...)
{
   return(.Call("int64ToInteger", x, PACKAGE="csvread"))
}

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method is.na int64
is.na.int64 <- function(x, ...)
{
   return(.Call("isInt64NA", x, PACKAGE="csvread"))
}

#-------------------------------------------------------------------------------

# \code{as.data.frame.int64} functions the same as \code{as.data.frame.vector}.
# @param x Object to be coerced or tested
# @param ... Further arguments passed to or from other methods.
# @param row.names \code{NULL} or a character vector giving the row names for the
#           data frame.  Missing values are not allowed.
# @param optional logical. If \code{TRUE}, setting row names and converting column
#           names (to syntactic names: see \code{make.names}) is optional.
#' @rdname int64
#' @export
#' @method as.data.frame int64
as.data.frame.int64 <- function(x, ...) as.data.frame.vector(x, ...)

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method as.list int64
as.list.int64 <- function(x, ...)
	lapply(seq_along(x), function(i) x[i])

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method c int64
c.int64 <- function(...)
	structure(c(unlist(lapply(list(...), unclass))), class="int64")

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method is.numeric int64
is.numeric.int64 <- function(x) FALSE

#-------------------------------------------------------------------------------

#' @rdname int64
#' @export
#' @method rep int64
rep.int64 <- function(x, ...)
{
	y <- NextMethod()
	structure(y, class="int64")
}

#-------------------------------------------------------------------------------

