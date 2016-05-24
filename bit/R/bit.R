# 1-bit boolean vectors for R
# (c) 2008-2009 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk

# currently |.bit and |.bitwhich are bypassed if we ask for bit | bitwhich
# xx explore/write Ops.bit Ops.bitwhich
# xx bit_extract should be comlemented with 

# source("C:/mwp/eanalysis/bit/R/bit.R")

#! \name{bit-package}
#! \alias{bit-package}
#! \alias{bit}
#! \alias{print.bit}
#! \docType{package}
#! \title{
#!    A class for vectors of 1-bit booleans
#! }
#! \description{
#! Package 'bit' provides bitmapped vectors of booleans (no NAs),
#! coercion from and to logicals, integers and integer subscripts;
#! fast boolean operators and fast summary statistics. \cr
#!
#! With bit vectors you can store true binary booleans \{FALSE,TRUE\} at the expense
#! of 1 bit only, on a 32 bit architecture this means factor 32 less RAM and
#! factor 32 more speed on boolean operations. With this speed gain it even
#! pays-off to convert to bit in order to avoid a single boolean operation on
#! logicals or a single set operation on (longer) integer subscripts, the pay-off
#! is dramatic when such components are used more than once. \cr
#!
#! Reading from and writing to bit is approximately as fast as accessing standard
#! logicals - mostly due to R's time for memory allocation. The package allows to
#! work with pre-allocated memory for return values by calling .Call() directly:
#! when evaluating the speed of C-access with pre-allocated vector memory, coping
#! from bit to logical requires only 70\% of the time for copying from logical to
#! logical; and copying from logical to bit comes at a performance penalty of 150\%. \cr
#!
#! Since bit objects cannot be used as subsripts in R, a second class 'bitwhich'
#! allows to store selections as efficiently as possible with standard R types.
#! This is usefull either to represent parts of bit objects or to represent
#! very asymetric selections.  \cr
#!
#! Class 'ri' (range index) allows to select ranges of positions for  chunked processing:
#! all three classes 'bit', 'bitwhich' and 'ri' can be used for subsetting 'ff' objects (ff-2.1.0 and higher).
#! }
#! \usage{
#!  bit(length)
#!  \method{print}{bit}(x, \dots)
#! }
#! \arguments{
#!   \item{length}{ length of vector in bits }
#!   \item{x}{ a bit vector }
#!   \item{\dots}{ further arguments to print }
#! }
#! \details{
#! \tabular{ll}{
#!    Package: \tab bit\cr
#!    Type: \tab Package\cr
#!    Version: \tab 1.1.0\cr
#!    Date: \tab 2012-06-05\cr
#!    License: \tab GPL-2\cr
#!    LazyLoad: \tab yes\cr
#!    Encoding: \tab latin1\cr
#! }
#!
#! Index:
#! \tabular{rrrrl}{
#!    \bold{bit function}           \tab \bold{bitwhich function}          \tab \bold{ri function}                \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{.BITS}                  \tab                                   \tab                                   \tab \code{\link{globalenv}}   \tab variable holding number of bits on this system \cr
#!    \code{\link{bit_init}}        \tab                                   \tab                                   \tab \code{\link{.First.lib}}  \tab initially allocate bit-masks (done in .First.lib) \cr
#!    \code{\link{bit_done}}        \tab                                   \tab                                   \tab \code{\link{.Last.lib}}   \tab finally de-allocate bit-masks (done in .Last.lib) \cr
#!    \code{\link{bit}}             \tab \code{\link{bitwhich}}            \tab \code{\link{ri}}                  \tab \code{\link{logical}}     \tab create bit object \cr
#!    \code{\link{print.bit}}       \tab \code{\link{print.bitwhich}}      \tab \code{\link{print.ri}}            \tab \code{\link{print}}       \tab print bit vector \cr
#!    \code{\link{length.bit}}      \tab \code{\link{length.bitwhich}}     \tab \code{\link{length.ri}}           \tab \code{\link{length}}      \tab get length of bit vector \cr
#!    \code{\link{length<-.bit}}    \tab \code{\link{length<-.bitwhich}}   \tab                                   \tab \code{\link{length<-}}    \tab change length of bit vector \cr
#!    \code{\link{c.bit}}           \tab \code{\link{c.bitwhich}}          \tab                                   \tab \code{\link{c}}           \tab concatenate bit vectors \cr
#!    \code{\link{is.bit}}          \tab \code{\link{is.bitwhich}}         \tab \code{\link{is.ri}}               \tab \code{\link{is.logical}}  \tab test for bit class \cr
#!    \code{\link{as.bit}}          \tab \code{\link{as.bitwhich}}         \tab                                   \tab \code{\link{as.logical}}  \tab generically coerce to bit or bitwhich \cr
#!    \code{\link{as.bit.logical}}  \tab \code{\link{as.bitwhich.logical}} \tab                                   \tab \code{\link{logical}}     \tab coerce logical to bit vector (FALSE => FALSE, c(NA, TRUE) => TRUE) \cr
#!    \code{\link{as.bit.integer}}  \tab \code{\link{as.bitwhich.integer}} \tab                                   \tab \code{\link{integer}}     \tab coerce integer to bit vector (0 => FALSE, ELSE => TRUE) \cr
#!    \code{\link{as.bit.double}}   \tab \code{\link{as.bitwhich.double}}  \tab                                   \tab \code{\link{double}}      \tab coerce double to bit vector (0 => FALSE, ELSE => TRUE) \cr
#!    \code{\link{as.double.bit}}   \tab \code{\link{as.double.bitwhich}}  \tab \code{\link{as.double.ri}}        \tab \code{\link{as.double}}   \tab coerce bit vector to double (0/1) \cr
#!    \code{\link{as.integer.bit}}  \tab \code{\link{as.integer.bitwhich}} \tab \code{\link{as.integer.ri}}       \tab \code{\link{as.integer}}  \tab coerce bit vector to integer (0L/1L) \cr
#!    \code{\link{as.logical.bit}}  \tab \code{\link{as.logical.bitwhich}} \tab \code{\link{as.logical.ri}}       \tab \code{\link{as.logical}}  \tab coerce bit vector to logical (FALSE/TRUE) \cr
#!    \code{\link{as.which.bit}}    \tab \code{\link{as.which.bitwhich}}   \tab \code{\link{as.which.ri}}         \tab \code{\link{as.which}}    \tab coerce bit vector to positive integer subscripts\cr
#!    \code{\link{as.bit.which}}    \tab \code{\link{as.bitwhich.which}}   \tab                                   \tab \code{\link{bitwhich}}    \tab coerce integer subscripts to bit vector \cr
#!    \code{\link{as.bit.bitwhich}} \tab \code{\link{as.bitwhich.bitwhich}}\tab                                   \tab                           \tab coerce from bitwhich  \cr
#!    \code{\link{as.bit.bit}}      \tab \code{\link{as.bitwhich.bit}}     \tab                                   \tab \code{\link{UseMethod}}   \tab coerce from bit \cr
#!    \code{\link{as.bit.ri}}       \tab \code{\link{as.bitwhich.ri}}      \tab                                   \tab                           \tab coerce from range index \cr
#!    \code{\link[ff]{as.bit.ff}}   \tab                                   \tab                                   \tab \code{\link[ff]{ff}}      \tab coerce ff boolean to bit vector \cr
#!    \code{\link[ff]{as.ff.bit}}   \tab                                   \tab                                   \tab \code{\link[ff]{as.ff}}   \tab coerce bit vector to ff boolean \cr
#!    \code{\link[ff]{as.hi.bit}}   \tab \code{\link[ff]{as.hi.bitwhich}}  \tab \code{\link[ff]{as.hi.ri}}        \tab \code{\link[ff]{as.hi}}   \tab coerce to hybrid index (requires package ff) \cr
#!    \code{\link[ff]{as.bit.hi}}   \tab \code{\link[ff]{as.bitwhich.hi}}  \tab                                   \tab                           \tab coerce from hybrid index (requires package ff) \cr
#!    \code{\link{[[.bit}}          \tab                                   \tab                                   \tab \code{\link{[[}}          \tab get single bit (index checked) \cr
#!    \code{\link{[[<-.bit}}        \tab                                   \tab                                   \tab \code{\link{[[<-}}        \tab set single bit (index checked) \cr
#!    \code{\link{[.bit}}           \tab                                   \tab                                   \tab \code{\link{[}}           \tab get vector of bits (unchecked) \cr
#!    \code{\link{[<-.bit}}         \tab                                   \tab                                   \tab \code{\link{[<-}}         \tab set vector of bits (unchecked) \cr
#!    \code{\link{!.bit}}           \tab \code{\link{!.bitwhich}}          \tab (works as second arg in           \tab \code{\link{!}}           \tab boolean NOT on bit \cr
#!    \code{\link{&.bit}}           \tab \code{\link{&.bitwhich}}          \tab  bit and bitwhich ops)            \tab \code{\link{&}}           \tab boolean AND on bit \cr
#!    \code{\link{|.bit}}           \tab \code{\link{|.bitwhich}}          \tab                                   \tab \code{\link{|}}           \tab boolean OR on bit \cr
#!    \code{\link{xor.bit}}         \tab \code{\link{xor.bitwhich}}        \tab                                   \tab \code{\link{xor}}         \tab boolean XOR on bit \cr
#!    \code{\link{!=.bit}}          \tab \code{\link{!=.bitwhich}}         \tab                                   \tab \code{\link{!=}}          \tab boolean unequality (same as XOR) \cr
#!    \code{\link{==.bit}}          \tab \code{\link{==.bitwhich}}         \tab                                   \tab \code{\link{==}}          \tab boolean equality \cr
#!    \code{\link{all.bit}}         \tab \code{\link{all.bitwhich}}        \tab \code{\link{all.ri}}              \tab \code{\link{all}}         \tab aggregate AND \cr
#!    \code{\link{any.bit}}         \tab \code{\link{any.bitwhich}}        \tab \code{\link{any.ri}}              \tab \code{\link{any}}         \tab aggregate OR \cr
#!    \code{\link{min.bit}}         \tab \code{\link{min.bitwhich}}        \tab \code{\link{min.ri}}              \tab \code{\link{min}}         \tab aggregate MIN (first TRUE position) \cr
#!    \code{\link{max.bit}}         \tab \code{\link{max.bitwhich}}        \tab \code{\link{max.ri}}              \tab \code{\link{max}}         \tab aggregate MAX (last TRUE position) \cr
#!    \code{\link{range.bit}}       \tab \code{\link{range.bitwhich}}      \tab \code{\link{range.ri}}            \tab \code{\link{range}}       \tab aggregate [MIN,MAX] \cr
#!    \code{\link{sum.bit}}         \tab \code{\link{sum.bitwhich}}        \tab \code{\link{sum.ri}}              \tab \code{\link{sum}}         \tab aggregate SUM (count of TRUE) \cr
#!    \code{\link{summary.bit}}     \tab \code{\link{summary.bitwhich}}    \tab \code{\link{summary.ri}}          \tab \code{\link{tabulate}}    \tab aggregate c(nFALSE, nTRUE, minRange, maxRange) \cr
#!    \code{\link{regtest.bit}}     \tab                                   \tab                                   \tab                           \tab regressiontests for the package \cr
#!  }
#! }
#! \value{
#!   \code{bit} returns a vector of integer sufficiently long to store 'length' bits
#!   (but not longer) with an attribute 'n' and class 'bit'
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#!
#! Maintainer: Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \note{
#!   Currently operations on bit objects have some overhead from R-calls. Do expect speed gains for vectors
#!   of length ~ 10000 or longer. \cr
#!   Since this package was created for high performance purposes, only positive integer subscripts are allowed:
#!   The '[.bit' and '[<-.bit' methods don't check whether the subscripts are positive integers in the allowed range.
#!   All R-functions behave as expected - i.e. they do not change their arguments and create new return values.
#!   If you want to save the time for return value memory allocation, you must use \code{\link{.Call}} directly
#!   (see the dontrun example in \code{\link{sum.bit}}).
#!   Note that the package has not been tested under 64 bit.
#!   Note also that the mapping of NAs to TRUE differs from the mapping of NAs to FALSE
#!   in \code{\link[ff]{vmode}="boolean"} in package ff (and one of the two may change in the future).
#! }
#! \keyword{ package }
#! \keyword{ classes }
#! \keyword{ logic }
#! \seealso{ \code{\link{logical}} in base R and \code{\link[ff]{vmode}} in package 'ff' }
#! \examples{
#!   x <- bit(12)                                 # create bit vector
#!   x                                            # autoprint bit vector
#!   length(x) <- 16                              # change length
#!   length(x)                                    # get length
#!   x[[2]]                                       # extract single element
#!   x[[2]] <- TRUE                               # replace single element
#!   x[1:2]                                       # extract parts of bit vector
#!   x[1:2] <- TRUE                               # replace parts of bit vector
#!   as.which(x)                                  # coerce bit to subscripts
#!   x <- as.bit.which(3:4, 4)                    # coerce subscripts to bit
#!   as.logical(x)                                # coerce bit to logical
#!   y <- as.bit(c(FALSE, TRUE, FALSE, TRUE))     # coerce logical to bit
#!   is.bit(y)                                    # test for bit
#!   !x                                           # boolean NOT
#!   x & y                                        # boolean AND
#!   x | y                                        # boolean OR
#!   xor(x, y)                                    # boolean Exclusive OR
#!   x != y                                       # boolean unequality (same as xor)
#!   x == y                                       # boolean equality
#!   all(x)                                       # aggregate AND
#!   any(x)                                       # aggregate OR
#!   min(x)                                       # aggregate MIN (integer version of ALL)
#!   max(x)                                       # aggregate MAX (integer version of ANY)
#!   range(x)                                     # aggregate [MIN,MAX]
#!   sum(x)                                       # aggregate SUM (count of TRUE)
#!   summary(x)                                   # aggregate count of FALSE and TRUE
#!
#!   \dontrun{
#!     message("\nEven for a single boolean operation transforming logical to bit pays off")
#!     n <- 10000000
#!     x <- sample(c(FALSE, TRUE), n, TRUE)
#!     y <- sample(c(FALSE, TRUE), n, TRUE)
#!     system.time(x|y)
#!     system.time({
#!        x <- as.bit(x)
#!        y <- as.bit(y)
#!     })
#!     system.time( z <- x | y )
#!     system.time( as.logical(z) )
#!     message("Even more so if multiple operations are needed :-)")
#!
#!     message("\nEven for a single set operation transforming subscripts to bit pays off\n")
#!     n <- 10000000
#!     x <- sample(n, n/2)
#!     y <- sample(n, n/2)
#!     system.time( union(x,y) )
#!     system.time({
#!      x <- as.bit.which(x, n)
#!      y <- as.bit.which(y, n)
#!     })
#!     system.time( as.which.bit( x | y ) )
#!     message("Even more so if multiple operations are needed :-)")
#!
#!     message("\nSome timings WITH memory allocation")
#!     n <- 2000000
#!     l <- sample(c(FALSE, TRUE), n, TRUE)
#!     # copy logical to logical
#!     system.time(for(i in 1:100){  # 0.0112
#!        l2 <- l
#!        l2[1] <- TRUE   # force new memory allocation (copy on modify)
#!        rm(l2)
#!     })/100
#!     # copy logical to bit
#!     system.time(for(i in 1:100){  # 0.0123
#!        b <- as.bit(l)
#!        rm(b)
#!     })/100
#!     # copy bit to logical
#!     b <- as.bit(l)
#!     system.time(for(i in 1:100){  # 0.009
#!        l2 <- as.logical(b)
#!        rm(l2)
#!     })/100
#!     # copy bit to bit
#!     b <- as.bit(l)
#!     system.time(for(i in 1:100){  # 0.009
#!        b2 <- b
#!        b2[1] <- TRUE   # force new memory allocation (copy on modify)
#!        rm(b2)
#!     })/100
#!
#!
#!     l2 <- l
#!     # replace logical by TRUE
#!     system.time(for(i in 1:100){
#!        l[] <- TRUE
#!     })/100
#!     # replace bit by TRUE (NOTE that we recycle the assignment  
#!		 # value on R side == memory allocation and assignment first)
#!     system.time(for(i in 1:100){
#!        b[] <- TRUE
#!     })/100
#!     # THUS the following is faster
#!     system.time(for(i in 1:100){
#!        b <- !bit(n)
#!     })/100
#!
#!     # replace logical by logical
#!     system.time(for(i in 1:100){
#!        l[] <- l2
#!     })/100
#!     # replace bit by logical
#!     system.time(for(i in 1:100){
#!        b[] <- l2
#!     })/100
#!     # extract logical
#!     system.time(for(i in 1:100){
#!        l2[]
#!     })/100
#!     # extract bit
#!     system.time(for(i in 1:100){
#!        b[]
#!     })/100
#!
#!     message("\nSome timings WITHOUT memory allocation (Serge, that's for you)")
#!     n <- 2000000L
#!     l <- sample(c(FALSE, TRUE), n, TRUE)
#!     b <- as.bit(l)
#!     # read from logical, write to logical
#!     l2 <- logical(n)
#!     system.time(for(i in 1:100).Call("R_filter_getset", l, l2, PACKAGE="bit")) / 100
#!     # read from bit, write to logical
#!     l2 <- logical(n)
#!     system.time(for(i in 1:100).Call("R_bit_get", b, l2, c(1L, n), PACKAGE="bit")) / 100
#!     # read from logical, write to bit
#!     system.time(for(i in 1:100).Call("R_bit_set", b, l2, c(1L, n), PACKAGE="bit")) / 100
#!
#!   }
#! }

#was wrong because C-code uses int: .BITS <- 8L * .Machine$sizeof.pointer
.BITS <- 32L



#! \name{bit_init}
#! \alias{bit_init}
#! \alias{bit_done}
#! \alias{.BITS}
#! \title{ Initializing bit masks }
#! \description{
#!   Functions to allocate (and de-allocate) bit masks
#! }
#! \usage{
#!   bit_init()
#!   bit_done()
#! }
#! \details{
#!   The C-code operates with bit masks.
#!   The memory for these is allocated dynamically.
#!   \code{bit_init} is called by \code{\link{.First.lib}}
#!   and \code{bit_done} is called by \code{\link{.Last.lib}}.
#!   You don't need to care about these under normal circumstances.
#! }
#! \value{
#!   NULL
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{bit}}  }
#! \examples{
#!   bit_done()
#!   bit_init()
#! }
#! \keyword{ classes }
#! \keyword{ logic }



# initialize and finalize the bit-mask vectors used in C

bit_init <- function()
  .Call("R_bit_init", .BITS, PACKAGE="bit")

bit_done <- function()
  .Call("R_bit_done", PACKAGE="bit")



# creator for empty bit vector
bit <- function(length){
  length <- as.integer(length)
  if (length %% .BITS)
    n <- length %/% .BITS + 1L
  else
    n <- length %/% .BITS
  x <- integer(n)
  #physical(x) <- list(vmode="boolean")
  #virtual(x)  <- list(Length=length)
  #class(x) <- "bit"
  # tuning
  p <- list()
  v <- list()
  attributes(p) <- list(vmode="boolean", class="physical")
  attributes(v) <- list(Length=length, class="virtual")
  attributes(x) <- list(physical=p, virtual=v, class="bit")
  x
}



print.bit <- function(x, ...){
  n <- length(x)
  cat("bit length=", n, " occupying only ", length(unclass(x)), " integers\n", sep="")
  if (n>16){
    y <- c(x[1:8], "..", x[(n-7L):n])
    names(y) <- c(1:8, "", (n-7L):n)
    print(y, quote=FALSE, ...)
  }else if(n){
    y <- c(x[])
    names(y) <- c(1:n)
    print(y, quote=FALSE, ...)
  }
}


#! \name{bitwhich}
#! \alias{bitwhich}
#! \alias{print.bitwhich}
#! \title{ A class for vectors representing asymetric selections }
#! \description{
#!   A bitwhich object like the result of \code{\link{which}} and \code{\link{as.which}} does represent integer subscript positions,
#!   but bitwhich objects represent some subscripts rather with negative integers, if this needs less space.
#!   The extreme cases of selecting all/none subscripts are represented by TRUE/FALSE.
#!   This needs less RAM compared to \code{\link{logical}} (and often less than \code{\link{as.which}}).
#!   Logical operations are fast if the selection is asymetric (only few or almost all selected).
#! }
#! \usage{
#! bitwhich(maxindex, poslength = NULL, x = NULL)
#! }
#! \arguments{
#!   \item{maxindex}{ the length of the vector (sum of all TRUEs and FALSEs) }
#!   \item{poslength}{ Only use if x is not NULL: the sum of all TRUEs }
#!   \item{x}{ Default NULL or FALSE or unique negative integers or unique positive integers or TRUE}
#! }
#! \value{
#!   An object of class 'bitwhich' carrying two attributes
#!   \item{maxindex}{ see above }
#!   \item{poslength}{ see above }
#! }
#! \details{
#!   class 'bitwhich' represents a boolean selection in one of the following ways
#!   \itemize{
#!    \item FALSE to select nothing
#!    \item TRUE to select everything
#!    \item unique positive integers to select those
#!    \item unique negative integers to exclude those
#!   }
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{as.bitwhich}}, \code{\link{as.which}}, \code{\link{bit}} }
#! \examples{
#!  bitwhich(12, x=c(1,3), poslength=2)
#!  bitwhich(12, x=-c(1,3), poslength=10)
#! }
#! \keyword{ classes }
#! \keyword{ logic }


bitwhich <- function(maxindex, poslength=NULL, x=NULL){
  if (is.null(x)){
    x <- FALSE
    poslength <- 0L
  }else{
    poslength <- as.integer(poslength)
  }
  attr(x, "maxindex") <- as.integer(maxindex)
  attr(x, "poslength") <- poslength
  # NOTE: here we want one (1) copy of x to not modify argument x 
	# therefore we did not replace the oldClass assignment with a call to setttattr
  oldClass(x) <- "bitwhich"
  x
}

print.bitwhich <- function(x, ...){
  cat("bitwhich: ", sum(x), "/", length(x), "\n", sep="")
}



#! \name{is.bit}
#! \alias{is.ri}
#! \alias{is.bit}
#! \alias{is.bitwhich}
#! \title{ Testing for bit, bitwhich and ri selection classes }
#! \description{
#!   Test whether an object inherits from 'ri', 'bit' or 'bitwhich'
#! }
#! \usage{
#! is.ri(x)
#! is.bit(x)
#! is.bitwhich(x)
#! }
#! \arguments{
#!   \item{x}{ an R object of unknown type }
#! }
#! \value{
#!   TRUE or FALSE
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{is.logical}}, \code{\link{bit}}, \code{\link{bitwhich}} }
#! \examples{
#!  is.ri(TRUE)
#!  is.ri(ri(1,4,12))
#!  is.bit(TRUE)
#!  is.bitwhich(TRUE)
#!  is.bit(as.bit(TRUE))
#!  is.bitwhich(as.bitwhich(TRUE))
#! }
#! \keyword{ classes }
#! \keyword{ logic }

is.ri <- function(x)
  inherits(x, "ri")

is.bit <- function(x)
  inherits(x, "bit")

is.bitwhich <- function(x)
  inherits(x, "bitwhich")



#! \name{length.bit}
#! \alias{length.bit}
#! \alias{length.bitwhich}
#! \alias{length.ri}
#! \alias{length<-.bit}
#! \alias{length<-.bitwhich}
#! \title{ Getting and setting length of bit, bitwhich and ri objects }
#! \description{
#!   Query the number of bits in a \code{\link{bit}} vector or change the number of bits in a bit vector. \cr
#!   Query the number of bits in a \code{\link{bitwhich}} vector or change the number of bits in a bit vector. \cr
#! }
#! \usage{
#! \method{length}{ri}(x)
#! \method{length}{bit}(x)
#! \method{length}{bitwhich}(x)
#! \method{length}{bit}(x) <- value
#! \method{length}{bitwhich}(x) <- value
#! }
#! \arguments{
#!   \item{x}{ a \code{\link{bit}}, \code{\link{bitwhich}} or \code{\link{ri}} object }
#!   \item{value}{ the new number of bits }
#! }
#! \details{
#!   NOTE that the length does NOT reflect the number of selected (\code{TRUE}) bits, it reflects the sum of both, \code{TRUE} and \code{FALSE} bits.
#!   Increasing the length of a \code{\link{bit}} object will set new bits to \code{FALSE}.
#!   The behaviour of increasing the length of a \code{\link{bitwhich}} object is different and depends on the content of the object:
#!   \itemize{
#!    \item{TRUE}{all included, new bits are set to \code{TRUE}}
#!    \item{positive integers}{some included, new bits are set to \code{FALSE}}
#!    \item{negative integers}{some excluded, new bits are set to \code{TRUE}}
#!    \item{FALSE}{all excluded:, new bits are set to \code{FALSE}}
#!   }
#!   Decreasing the length of bit or bitwhich removes any previous information about the status bits above the new length.
#! }
#! \value{
#!   the length  A bit vector with the new length
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{length}}, \code{\link[=sum.bit]{sum}}, \code{\link[ff]{poslength}}, \code{\link[ff]{maxindex}} }
#! \examples{
#!   stopifnot(length(ri(1, 1, 32))==32)
#!
#!   x <- as.bit(ri(32, 32, 32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==1)
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==0)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==0)
#!
#!   x <- as.bit(ri(1, 1, 32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==1)
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==1)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==1)
#!
#!   x <- as.bitwhich(bit(32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==0)
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==0)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==0)
#!
#!   x <- as.bitwhich(!bit(32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==32)
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==16)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==32)
#!
#!   x <- as.bitwhich(ri(32, 32, 32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==1)
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==0)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==0)
#!
#!   x <- as.bitwhich(ri(2, 32, 32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==31)
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==15)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==31)
#!
#!   x <- as.bitwhich(ri(1, 1, 32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==1)
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==1)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==1)
#!
#!   x <- as.bitwhich(ri(1, 31, 32))
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==31)
#!   message("NOTE the change from 'some excluded' to 'all excluded' here")
#!   length(x) <- 16
#!   stopifnot(length(x)==16)
#!   stopifnot(sum(x)==16)
#!   length(x) <- 32
#!   stopifnot(length(x)==32)
#!   stopifnot(sum(x)==32)
#! }
#! \keyword{ classes }
#! \keyword{ logic }


length.bit <- function(x)
  virtual(x)$Length

"length<-.bit" <- function(x, value){
  if (value!=length(x)){
    value <- as.integer(value)
    dn <- value %% .BITS
    if (dn){
      n <- value %/% .BITS + 1L
      .Call("R_bit_replace", x, (value+1L):(value+dn), logical(dn), FALSE, PACKAGE="bit")
    }else{
      n <- value %/% .BITS
    }
    #pattr <- physical(x)
    #vattr <- virtual(x)
    pattr <- attr(x, "physical")
    vattr <- attr(x, "virtual")
    cl <- oldClass(x)
    #x <- unclass(x)
    attr(x, "class") <- NULL
    length(x) <- n
    #vattr$Length <- value
    attr(vattr, "Length") <- value
    #physical(x) <- pattr
    #virtual(x) <- vattr
    #class(x) <- cl
    attr(x, "physical") <- pattr
    attr(x, "virtual") <- vattr
    attr(x, "class") <- cl
    x
  }else
    x
}


length.bitwhich <- function(x)
  attr(x, "maxindex")

"length<-.bitwhich" <- function(x, value){
  if (value!=length(x)){
    value <- as.integer(value)
    if (is.integer(x)){
      cl <- oldClass(x)
      oldClass(x) <- NULL
      if (x[1]>0){
        x <- x[x <= value]
        l <- length(x)
        if (l==0)
          x <- FALSE
        else if (l==value)
          x <- TRUE
        else if (l>(value%/%2L))
          x <- -as.integer(seq_len(value))[-x]
        attr(x, "poslength") <- l
      }else{
        x <- x[x >= -value]
        l <- length(x)
        if (l==0)
          x <- TRUE
        else if (l==value)
          x <- FALSE
        else if (!((value-l)>(value%/%2L)))
          x <- -as.integer(seq_len(value))[-x]
        attr(x, "poslength") <- value - l
      }
      oldClass(x) <- cl
    }else if(x){
      attr(x, "poslength") <- value
    }
    attr(x, "maxindex") <- value
  }
  x
}






#! \name{c.bit}
#! \alias{c.bit}
#! \alias{c.bitwhich}
#! \title{ Concatenating bit and bitwhich vectors }
#! \description{
#!   Creating new bit by concatenating bit vectors
#! }
#! \usage{
#! \method{c}{bit}(\dots)
#! \method{c}{bitwhich}(\dots)
#! }
#! \arguments{
#!   \item{\dots}{ bit objects }
#! }
#! \value{
#!   An object of class 'bit'
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{c}}, \code{\link{bit}} , \code{\link{bitwhich}} }
#! \examples{
#!  c(bit(4), bit(4))
#! }
#! \keyword{ classes }
#! \keyword{ logic }

c.bit <- function(...){
  l <- list(...)
  nl <- length(l)
  nold <- sapply(l, length)
  nnew <- sum(nold)
  ncum <- cumsum(nold)
  offsets <- c(0L, ncum[-length(ncum)])
  x <- bit(nnew)
  for (i in as.integer(seq.int(from=1, to=nl, by=1))){
    b <- as.bit(l[[i]])
    .Call("R_bit_shiftcopy", bsource_=b, btarget_=x, otarget_=offsets[i], n_=nold[i], FALSE, PACKAGE="bit")
  }
  x
}

c.bitwhich <- function(...){
  l <- list(...)
  if (length(l)==1)
    l[[1]]
  else
    as.bitwhich(do.call("c", lapply(l, as.bit)))
}



#! \name{as.bit}
#! \alias{as.bit}
#! \alias{as.bit.bit}
#! \alias{as.bit.logical}
#! \alias{as.bit.integer}
#! \alias{as.bit.double}
#! \alias{as.bit.bitwhich}
#! \alias{as.bit.which}
#! \alias{as.bit.ri}
#! \title{ Coercing to bit }
#! \description{
#!   Coercing to bit vector
#! }
#! \usage{
#! as.bit(x, \dots)
#! \method{as.bit}{bit}(x, \dots)
#! \method{as.bit}{logical}(x, \dots)
#! \method{as.bit}{integer}(x, \dots)
#! \method{as.bit}{bitwhich}(x, \dots)
#! \method{as.bit}{which}(x, length, \dots)
#! \method{as.bit}{ri}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an object of class \code{\link{bit}}, \code{\link{logical}}, \code{\link{integer}}, \code{\link{bitwhich}} or an integer from \code{\link{as.which}} or a boolean \code{\link[ff:vmode]{ff}} }
#!   \item{length}{ the length of the new bit vector }
#!   \item{\dots}{ further arguments }
#! }
#! \details{
#!   Coercing to bit is quite fast because we use a double loop that fixes each word in a processor register
#! }
#! \note{
#!   Zero is coerced to FALSE, all other numbers including NA are coerced to TRUE.
#!   This differs from the NA-to-FALSE coercion in package ff and may change in the future.
#! }
#! \value{
#!   \code{is.bit} returns FALSE or TRUE, \code{as.bit} returns a vector of class 'bit'
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{bit}}, \code{\link[bit:as.logical.bit]{as.logical}} }
#! \examples{
#!   x <- as.bit(c(FALSE, NA, TRUE))
#!   as.bit(x)
#!   as.bit.which(c(1,3,4), 12)
#! }
#! \keyword{ classes }
#! \keyword{ logic }

as.bit.bit <- function(x, ...)
  x

as.bit.logical <- function(x, ...){
  n <- length(x)
  b <- bit(n)
  .Call("R_bit_set", b, x, c(1L, n), PACKAGE="bit")
}
as.bit.integer <- function(x, ...){
  n <- length(x)
  b <- bit(n)
  .Call("R_bit_set_integer", b, x, c(1L, n), PACKAGE="bit")
}
as.bit.double <- function(x, ...){
  n <- length(x)
  b <- bit(n)
  .Call("R_bit_set_integer", b, as.integer(x), c(1L, n), PACKAGE="bit")
}

as.bit.bitwhich <- function(x, ...){
  n <- length(x)
  p <- sum(x)
  b <- bit(n)
  if (is.logical(x)){
    if (p==n)
      b[] <- TRUE
  }else{
    oldClass(x) <- NULL
    x <- as.integer(x)
    if (x[1]<0){
      b[-x] <- TRUE  # remember that negative indices are not allowed (and the assignment value is recycled to the length of the index)
      b <- !b
    }else{
      b[x] <- TRUE
    }
  }
  b
}

as.bit.which <- function(x, length, ...){
  b <- bit(length)
  if (length(x)){
    x <- as.integer(x)
    if (x[1]<0){
      b[-x] <- TRUE  # remember that negative indices are not allowed (and the assignment value is recycled to the length of the index)
      b <- !b
    }else{
      b[x] <- TRUE
    }
  }
  b
}


as.bit.ri <- function(x, ...){
  b <- bit(length(x))
  b[x] <- TRUE
  b
}

#! \name{as.logical.bit}
#! \alias{as.logical.bit}
#! \alias{as.integer.bit}
#! \alias{as.double.bit}
#! \alias{as.logical.bitwhich}
#! \alias{as.integer.bitwhich}
#! \alias{as.double.bitwhich}
#! \alias{as.logical.ri}
#! \alias{as.integer.ri}
#! \alias{as.double.ri}
#! \title{ Coercion from bit, bitwhich and ri to logical, integer, double }
#! \description{
#!   Coercing from bit to logical, integer, which.
#! }
#! \usage{
#! \method{as.logical}{bit}(x, \dots)
#! \method{as.logical}{bitwhich}(x, \dots)
#! \method{as.logical}{ri}(x, \dots)
#! \method{as.integer}{bit}(x, \dots)
#! \method{as.integer}{bitwhich}(x, \dots)
#! \method{as.integer}{ri}(x, \dots)
#! \method{as.double}{bit}(x, \dots)
#! \method{as.double}{bitwhich}(x, \dots)
#! \method{as.double}{ri}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an object of class \code{\link{bit}}, \code{\link{bitwhich}} or \code{\link{ri}} }
#!   \item{\dots}{ ignored }
#! }
#! \details{
#!   Coercion from bit is quite fast because we use a double loop that fixes each word in a processor register.
#! }
#! \value{
#!   \code{\link{as.logical}} returns a vector of \code{FALSE, TRUE}, \code{\link{as.integer}} and \code{\link{as.double}} return a vector of \code{0,1}.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{as.bit}}, \code{\link{as.which}}, \code{\link{as.bitwhich}}, \code{\link[ff]{as.ff}}, \code{\link[ff]{as.hi}} }
#! \examples{
#!   x <- ri(2, 5, 10)
#!   y <- as.logical(x)
#!   y
#!   stopifnot(identical(y, as.logical(as.bit(x))))
#!   stopifnot(identical(y, as.logical(as.bitwhich(x))))
#!
#!   y <- as.integer(x)
#!   y
#!   stopifnot(identical(y, as.integer(as.logical(x))))
#!   stopifnot(identical(y, as.integer(as.bit(x))))
#!   stopifnot(identical(y, as.integer(as.bitwhich(x))))
#!
#!   y <- as.double(x)
#!   y
#!   stopifnot(identical(y, as.double(as.logical(x))))
#!   stopifnot(identical(y, as.double(as.bit(x))))
#!   stopifnot(identical(y, as.double(as.bitwhich(x))))
#! }
#! \keyword{ classes }
#! \keyword{ logic }


as.logical.bit <- function(x, ...){
  l <- logical(length(x))
  .Call("R_bit_get", x, l, c(1L, length(x)), PACKAGE="bit")
}
as.integer.bit <- function(x, ...){
  l <- integer(length(x))
  .Call("R_bit_get_integer", x, l, c(1L, length(x)), PACKAGE="bit")
}
as.double.bit <- function(x, ...){
  l <- integer(length(x))
  as.double(.Call("R_bit_get_integer", x, l, c(1L, length(x)), PACKAGE="bit"))
}

as.logical.ri <- function(x, ...){
  if (is.na(x[3]))
    stop("cannot coerce to logical from ri object with unknown maxindex")
  ret <- logical(x[3])
  ret[x[1]:x[2]] <- TRUE
  ret
}

as.integer.ri <- function(x, ...){
  if (is.na(x[3]))
    stop("cannot coerce to integer from ri object with unknown maxindex")
  ret <- integer(x[3])
  ret[x[1]:x[2]] <- 1L
  ret
}

as.double.ri <- function(x, ...){
  if (is.na(x[3]))
    stop("cannot coerce to integer from ri object with unknown maxindex")
  ret <- double(x[3])
  ret[x[1]:x[2]] <- 1
  ret
}




#! \name{as.which}
#! \alias{as.which}
#! \alias{as.which.default}
#! \alias{as.which.bitwhich}
#! \alias{as.which.bit}
#! \alias{as.which.ri}
#! \title{ Coercion to (positive) integer positions }
#! \description{
#!   Coercing to something like the result of which \code{\link{which}}
#! }
#! \usage{
#! as.which(x, \dots)
#! \method{as.which}{default}(x, \dots)
#! \method{as.which}{ri}(x, \dots)
#! \method{as.which}{bit}(x, range = NULL, \dots)
#! \method{as.which}{bitwhich}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an object of classes \code{\link{bit}}, \code{\link{bitwhich}}, \code{\link{ri}} or something on which \code{\link{which}} works }
#!   \item{range}{ a \code{\link{ri}} or an integer vector of length==2 giving a range restriction for chunked processing }
#!   \item{\dots}{ further arguments (passed to \code{\link{which}} for the default method, ignored otherwise) }
#! }
#! \details{
#!   \code{as.which.bit} returns a vector of subscripts with class 'which'
#! }
#! \value{
#!   a vector of class 'logical' or 'integer'
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{as.bit}}, \code{\link{as.logical}}, \code{\link{as.integer}}, \code{\link{as.which}}, \code{\link{as.bitwhich}}, \code{\link[ff]{as.ff}}, \code{\link[ff]{as.hi}} }
#! \examples{
#!   r <- ri(5, 20, 100)
#!   x <- as.which(r)
#!   x
#!
#!   stopifnot(identical(x, as.which(as.logical(r))))
#!   stopifnot(identical(x, as.which(as.bitwhich(r))))
#!   stopifnot(identical(x, as.which(as.bit(r))))
#! }
#! \keyword{ classes }
#! \keyword{ logic }



as.which.default <- function(x, ...){
  ret <- which(x)
  oldClass(ret) <- "which"
  ret
}

as.which.ri <- function(x, ...){
  ret <- x[1]:x[2]
  oldClass(ret) <- "which"
  ret
}

as.which.bit <- function(x, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(x))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(x))
      stop("illegal range")
  }
  s <- sum(x, range=range)
  n <- range[2] - range[1] + 1L
  if (s==0L){
    x <- integer()
  }else if (s==n){
    x <- as.integer(seq.int(from=range[1], to=range[2], by=1))
  }else
    x <- .Call("R_bit_which", x, s, range, negative=FALSE, PACKAGE="bit")
  oldClass(x) <- "which"
  x
}

as.which.bitwhich <- function(x, ...){
  if (is.logical(x)){
    if (unclass(x))
      x <- as.integer(seq_len(length(x)))
    else
      x <- integer()
  }else{
    if (x[[1]]<0)
      x <- as.integer(seq_len(length(x)))[x]
    else{
      attributes(x) <- NULL
    }
  }
  oldClass(x) <- "which"
  x
}



#! \name{as.bitwhich}
#! \alias{as.bitwhich}
#! \alias{as.bitwhich.bit}
#! \alias{as.bitwhich.bitwhich}
#! \alias{as.bitwhich.ri}
#! \alias{as.bitwhich.which}
#! \alias{as.bitwhich.integer}
#! \alias{as.bitwhich.double}
#! \alias{as.bitwhich.logical}
#! \title{ Coercing to bitwhich }
#! \description{
#!   Functions to coerce to bitwhich
#! }
#! \usage{
#! as.bitwhich(x, \dots)
#! \method{as.bitwhich}{bitwhich}(x, \dots)
#! \method{as.bitwhich}{ri}(x, \dots)
#! \method{as.bitwhich}{bit}(x, range=NULL, \dots)
#! \method{as.bitwhich}{which}(x, maxindex, \dots)
#! \method{as.bitwhich}{integer}(x, \dots)
#! \method{as.bitwhich}{double}(x, \dots)
#! \method{as.bitwhich}{logical}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ An object of class 'bitwhich', 'integer', 'logical' or 'bit' or an integer vector as resulting from 'which' }
#!   \item{maxindex}{ the length of the new bitwhich vector }
#!   \item{range}{ a \code{\link{ri}} or an integer vector of length==2 giving a range restriction for chunked processing }
#!   \item{\dots}{ further arguments }
#! }
#! \value{
#!   a value of class \code{\link{bitwhich}}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{bitwhich}}, \code{\link{as.bit}} }
#! \examples{
#!  as.bitwhich(c(FALSE, FALSE, FALSE))
#!  as.bitwhich(c(FALSE, FALSE, TRUE))
#!  as.bitwhich(c(FALSE, TRUE, TRUE))
#!  as.bitwhich(c(TRUE, TRUE, TRUE))
#! }
#! \keyword{ classes }
#! \keyword{ logic }

as.bitwhich.bitwhich <- function(x, ...){
  x
}

as.bitwhich.which <- function(x, maxindex, ...){
  poslength <- length(x)
  if (missing(maxindex))
    stop("you must provide maxindex with as.bitwhich.integer()")

  if (poslength==0)
    bitwhich(maxindex, poslength, FALSE)
  else if (poslength==maxindex)
    bitwhich(maxindex, poslength, TRUE)
  else if (poslength>(maxindex%/%2L)){
    bitwhich(maxindex, poslength, -as.integer(seq_len(maxindex))[-x])
  }else{
    bitwhich(maxindex, poslength, x)
  }
}

as.bitwhich.ri <- function(x, ...){
  poslength <- sum(x)
  maxindex <- length(x)

  if (poslength==0)
    bitwhich(maxindex, poslength, FALSE)
  else if (poslength==maxindex)
    bitwhich(maxindex, poslength, TRUE)
  else if (poslength>(maxindex%/%2L)){
    if (x[1]>1L) a <- 1:(x[1]-1L) else a <- integer()
    if (x[2]<maxindex) b <- (x[2]+1L):maxindex else b <- integer()
    bitwhich(maxindex, poslength, -c(a,b))
  }else{
    bitwhich(maxindex, poslength, x[1]:x[2])
  }
}


as.bitwhich.double <- as.bitwhich.integer <- function(x, ...)
  as.bitwhich(as.logical(x))

as.bitwhich.logical <- function(x, ...){
  poslength <- sum(x, na.rm=TRUE)
  maxindex <- length(x)

  if (poslength==0)
    bitwhich(maxindex, poslength, FALSE)
  else if (poslength==maxindex)
    bitwhich(maxindex, poslength, TRUE)
  else if (poslength>(maxindex%/%2L)){
    bitwhich(maxindex, poslength, -which(!x))
  }else{
    bitwhich(maxindex, poslength, which(x))
  }
}

as.bitwhich.bit <- function(x, range=NULL, ...){
  maxindex <- length(x)
  if (is.null(range))
    range <- c(1L, maxindex)
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>maxindex)
      stop("illegal range")
  }
  poslength <- sum(x, range=range, na.rm=TRUE)
  if (poslength==0)
    bitwhich(maxindex, poslength, FALSE)
  else if (poslength==maxindex)
    bitwhich(maxindex, poslength, TRUE)
  else{
    if (poslength>(maxindex%/%2L)){
      bitwhich(maxindex, poslength, .Call("R_bit_which", x, maxindex - poslength, range=range, negative=TRUE, PACKAGE="bit"))
    }else{
      bitwhich(maxindex, poslength, .Call("R_bit_which", x, poslength, range=range, negative=FALSE, PACKAGE="bit"))
    }
  }
}





as.integer.bitwhich <- function(x, ...){
  n <- length(x)
  if (is.logical(x)){
    if (sum(x)==n)
      rep(1L, n)
    else
      rep(0L, n)
  }else{
    ret <- integer(n)
    ret[x] <- 1L
    ret
  }
}
as.double.bitwhich <- function(x, ...){
  n <- length(x)
  if (is.logical(x)){
    if (sum(x)==n)
      rep(1, n)
    else
      rep(0, n)
  }else{
    ret <- double(n)
    ret[x] <- 1
    ret
  }
}


as.logical.bitwhich <- function(x, ...){
  n <- length(x)
  p <- sum(x)
  if (p==0){
    rep(FALSE, length(x))
  }else if (p==n){
    rep(TRUE, length(x))
  }else{
    ret <- logical(length(x))
    ret[x] <- TRUE
    ret
  }
}

# xx
#"[.bitwhich" <- function(x, i){
#  if (inherits(i, "bitwhich")){
#    nx <- length(x)
#    ni <- length(i)
#    px <- poslength(x)
#    pi <- poslength(i)
#    if (is.logical(x)){
#      if (is.logical(i)){
#        if (unclass(x) && unclass(i))
#          return(bitwhich())
#        else
#          return()
#      }else{
#      }
#    }else{
#      if (is.logical(i)){
#      }else{
#      }
#    }
#  }else
#    stop("not implemented")
#}




#! \name{LogicBit}
#! \alias{LogicBit}
#! \alias{!.bit}
#! \alias{!.bitwhich}
#! \alias{&.bit}
#! \alias{&.bitwhich}
#! \alias{|.bit}
#! \alias{|.bitwhich}
#! \alias{==.bit}
#! \alias{==.bitwhich}
#! \alias{!=.bit}
#! \alias{!=.bitwhich}
#! \alias{xor}
#! \alias{xor.default}
#! \alias{xor.bit}
#! \alias{xor.bitwhich}
#! \title{ Boolean operators and functions for class bit }
#! \description{
#!   Boolean 'negation', 'and', 'or' and 'exclusive or'.
#! }
#! \usage{
#! \method{!}{bit}(x)
#! \method{!}{bitwhich}(x)
#! \method{&}{bit}(e1, e2)
#! \method{&}{bitwhich}(e1, e2)
#! \method{|}{bit}(e1, e2)
#! \method{|}{bitwhich}(e1, e2)
#! \method{==}{bit}(e1, e2)
#! \method{==}{bitwhich}(e1, e2)
#! \method{!=}{bit}(e1, e2)
#! \method{!=}{bitwhich}(e1, e2)
#! xor(x, y)
#! \method{xor}{default}(x, y)
#! \method{xor}{bit}(x, y)
#! \method{xor}{bitwhich}(x, y)
#! }
#! \arguments{
#!   \item{x}{ a bit vector (or one logical vector in binary operators) }
#!   \item{y}{ a bit vector or an logical vector }
#!   \item{e1}{ a bit vector or an logical vector }
#!   \item{e2}{ a bit vector or an logical vector }
#! }
#! \details{
#!   Binary operators and function \code{xor} can combine 'bit' objects and 'logical' vectors.
#!   They do not recycle, thus the lengths of objects must match. Boolean operations on bit vectors are extremely fast
#!   because they are implemented using C's bitwise operators. If one argument is 'logical' it is converted to 'bit'. \cr
#!
#!   Binary operators and function \code{xor} can combine 'bitwhich' objects and other vectors.
#!   They do not recycle, thus the lengths of objects must match. Boolean operations on bitwhich vectors are fast
#!   if the distribution of TRUE and FALSE is very asymetric. If one argument is not 'bitwhich' it is converted to 'bitwhich'. \cr
#!
#!   The \code{xor} function has been made generic and \code{xor.default} has been implemented much faster than R's standard \code{\link[base]{xor}}.
#!   This was possible because actually boolean function \code{xor} and comparison operator \code{!=} do the same (even with NAs), and \code{!=} is much faster than the multiple calls in \code{(x | y) & !(x & y)}
#! }
#! \value{
#!   An object of class 'bit' (or 'bitwhich')
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{bit}}, \code{\link{Logic}} }
#! \examples{
#!   x <- as.bit(c(FALSE, FALSE, FALSE, NA, NA, NA, TRUE, TRUE, TRUE))
#!   yl <- c(FALSE, NA, TRUE, FALSE, NA, TRUE, FALSE, NA, TRUE)
#!   y <- as.bit(yl)
#!   !x
#!   x & y
#!   x | y
#!   xor(x, y)
#!   x != y
#!   x == y
#!   x & yl
#!   x | yl
#!   xor(x, yl)
#!   x != yl
#!   x == yl
#!
#!   x <- as.bitwhich(c(FALSE, FALSE, FALSE, NA, NA, NA, TRUE, TRUE, TRUE))
#!   yl <- c(FALSE, NA, TRUE, FALSE, NA, TRUE, FALSE, NA, TRUE)
#!   y <- as.bitwhich(yl)
#!   !x
#!   x & y
#!   x | y
#!   xor(x, y)
#!   x != y
#!   x == y
#!   x & yl
#!   x | yl
#!   xor(x, yl)
#!   x != yl
#!   x == yl
#! }
#! \keyword{ classes }
#! \keyword{ logic }


"!.bit" <- function(x){
  if (length(x)){
    ret <- x
    ret[1] <- ret[1]  # force duplication
    .Call("R_bit_not", ret, PACKAGE="bit")
  }else{
    x
  }
}


"&.bit" <- function(e1, e2){
  n <- length(e1)
  if(n!=length(e2))
    stop("length(e1) != length(e2)")
  e1 <- as.bit(e1)
  e2 <- as.bit(e2)
  ret <- bit(n)
  .Call("R_bit_and", e1, e2, ret, PACKAGE="bit")
}


"|.bit" <- function(e1, e2){
  n <- length(e1)
  if(n!=length(e2))
    stop("length(e1) != length(e2)")
  e1 <- as.bit(e1)
  e2 <- as.bit(e2)
  ret <- bit(n)
  .Call("R_bit_or", e1, e2, ret, PACKAGE="bit")
}

xor.default <- function(x,y){
  as.logical(x) != as.logical(y)
}

"xor.bit" <- function(x, y){
  n <- length(x)
  if(n!=length(y))
    stop("length(x) != length(y)")
  x <- as.bit(x)
  y <- as.bit(y)
  ret <- bit(n)
  .Call("R_bit_xor", x, y, ret, PACKAGE="bit")
}

"!=.bit" <- function(e1, e2){
  n <- length(e1)
  if(n!=length(e2))
    stop("length(e1) != length(e2)")
  e1 <- as.bit(e1)
  e2 <- as.bit(e2)
  ret <- bit(n)
  .Call("R_bit_xor", e1, e2, ret, PACKAGE="bit")
}

"==.bit" <- function(e1, e2){
  n <- length(e1)
  if(n!=length(e2))
    stop("length(e1) != length(e2)")
  e1 <- as.bit(e1)
  e2 <- as.bit(e2)
  ret <- bit(n)
  .Call("R_bit_equal", e1, e2, ret, PACKAGE="bit")
}



"!.bitwhich" <- function(x){
  n <- length(x)
  p <- sum(x)
  if (is.logical(x)){
    if (p==n){
      bitwhich(maxindex=n, poslength=0L, FALSE)
    }else{
      bitwhich(maxindex=n, poslength=n, TRUE)
    }
  }else{
    bitwhich(maxindex=n, poslength=n-p, -x)
  }
}


"&.bitwhich" <- function(e1, e2){
  e1 <- as.bitwhich(e1)
  e2 <- as.bitwhich(e2)
  n <- c(length(e1), length(e2))
  if(n[1]!=n[2])
    stop("length(e1) != length(e2)")
  p <- c(sum(e1), sum(e2))
  if (p[1]==0 || p[2]==0)
    return(bitwhich(n[1], 0L, FALSE))
  if (p[1]==n[1])
    return(e2)
  if (p[2]==n[2])
    return(e1)
  #negative <- p>(n%/%2L)
  negative <- c(e1[1]<0, e2[1]<0)
  if (negative[1]){
    if (negative[2]){
      ret <- union(e1, e2)
      return( bitwhich(maxindex=n[1], poslength=n[1]-length(ret), ret) )
    }else{
      ret <- setdiff(e2, !e1)
      return( bitwhich(maxindex=n[1], poslength=length(ret), if (length(ret)) ret else FALSE) )
    }
  }else{
    if (negative[2]){
      ret <- setdiff(e1, !e2)
      return( bitwhich(maxindex=n[1], poslength=length(ret), if (length(ret)) ret else FALSE) )
    }else{
      ret <- intersect(e1, e2)
      return( bitwhich(maxindex=n[1], poslength=length(ret), if (length(ret)) ret else FALSE) )
    }
  }
  #as.bitwhich(as.bit(e1) & as.bit(e2))
}


"|.bitwhich" <- function(e1, e2){
  e1 <- as.bitwhich(e1)
  e2 <- as.bitwhich(e2)
  n <- c(length(e1), length(e2))
  if(n[1]!=n[2])
    stop("length(e1) != length(e2)")
  p <- c(sum(e1), sum(e2))
  if (p[1]==n[1] || p[2]==n[2])
    return(bitwhich(n[1], n[1], TRUE))
  if (p[1]==0)
    return(e2)
  if (p[2]==0)
    return(e1)
  #negative <- p>(n%/%2L)
  negative <- c(e1[1]<0, e2[1]<0)
  if (negative[1]){
    if (negative[2]){
      ret <- intersect(e1, e2)
      return( bitwhich(maxindex=n[1], poslength=n[1]-length(ret), if (length(ret)) ret else TRUE) )
    }else{
      ret <- setdiff(e1, !e2)
      return( bitwhich(maxindex=n[1], poslength=n[1]-length(ret), if (length(ret)) ret else TRUE) )
    }
  }else{
    if (negative[2]){
      ret <- setdiff(e2, !e1)
      return( bitwhich(maxindex=n[1], poslength=n[1]-length(ret), if (length(ret)) ret else TRUE) )
    }else{
      ret <- union(e1, e2)
      return( bitwhich(maxindex=n[1], poslength=length(ret), ret) )
    }
  }
  #as.bitwhich(as.bit(e1) | as.bit(e2))
}

"xor.bitwhich" <- function(x, y){
  x <- as.bitwhich(x)
  y <- as.bitwhich(y)
  n <- c(length(x), length(y))
  if(n[1]!=n[2])
    stop("length(x) != length(y)")
  p <- c(sum(x), sum(y))
  if (p[1]==0)
    return(y)
  if (p[1]==n[1])
    return(!y)
  if (p[2]==0)
    return(x)
  if (p[2]==n[2])
    return(!x)
  #negative <- p>(n%/%2L)
  negative <- c(x[1]<0, y[1]<0)
  if (negative[1]){
    if (negative[2]){
      ret <- -union(setdiff(y, x), setdiff(x, y))
      return( bitwhich(maxindex=n[1], poslength=length(ret), if (length(ret)) ret else FALSE) )
    }else{
      ret <- union(-setdiff(y, !x), setdiff(x, !y))
      return( bitwhich(maxindex=n[1], poslength=n[1]-length(ret), if (length(ret)) ret else TRUE) )
    }
  }else{
    if (negative[2]){
      ret <- union(-setdiff(x, !y), setdiff(y, !x))
      return( bitwhich(maxindex=n[1], poslength=n[1]-length(ret), if (length(ret)) ret else TRUE) )
    }else{
      ret <- setdiff(union(x, y), intersect(x, y))
      return( bitwhich(maxindex=n[1], poslength=length(ret), if (length(ret)) ret else FALSE) )
    }
  }
  #as.bitwhich(xor(as.bit(x), as.bit(y)))
}


"!=.bitwhich" <- function(e1, e2)
xor.bitwhich(e1, e2)

"==.bitwhich" <- function(e1, e2)
!xor.bitwhich(e1, e2)







#! \name{Summary}
#! \alias{all.bit}
#! \alias{any.bit}
#! \alias{min.bit}
#! \alias{max.bit}
#! \alias{range.bit}
#! \alias{sum.bit}
#! \alias{summary.bit}
#! \alias{all.bitwhich}
#! \alias{any.bitwhich}
#! \alias{min.bitwhich}
#! \alias{max.bitwhich}
#! \alias{range.bitwhich}
#! \alias{sum.bitwhich}
#! \alias{summary.bitwhich}
#! \alias{all.ri}
#! \alias{any.ri}
#! \alias{min.ri}
#! \alias{max.ri}
#! \alias{range.ri}
#! \alias{sum.ri}
#! \alias{summary.ri}
#! \title{ Summaries of bit vectors }
#! \description{
#!   Fast aggregation functions for bit vectors.
#! }
#! \usage{
#! \method{all}{bit}(x, range = NULL, \dots)
#! \method{any}{bit}(x, range = NULL, \dots)
#! \method{min}{bit}(x, range = NULL, \dots)
#! \method{max}{bit}(x, range = NULL, \dots)
#! \method{range}{bit}(x, range = NULL, \dots)
#! \method{sum}{bit}(x, range = NULL, \dots)
#! \method{summary}{bit}(object, range = NULL, \dots)
#! \method{all}{bitwhich}(x, \dots)
#! \method{any}{bitwhich}(x, \dots)
#! \method{min}{bitwhich}(x, \dots)
#! \method{max}{bitwhich}(x, \dots)
#! \method{range}{bitwhich}(x, \dots)
#! \method{sum}{bitwhich}(x, \dots)
#! \method{summary}{bitwhich}(object, \dots)
#! \method{all}{ri}(x, \dots)
#! \method{any}{ri}(x, \dots)
#! \method{min}{ri}(x, \dots)
#! \method{max}{ri}(x, \dots)
#! \method{range}{ri}(x, \dots)
#! \method{sum}{ri}(x, \dots)
#! \method{summary}{ri}(object, \dots)
#! }
#! \arguments{
#!   \item{x}{ an object of class bit or bitwhich }
#!   \item{object}{ an object of class bit }
#!   \item{range}{ a \code{\link{ri}} or an integer vector of length==2 giving a range restriction for chunked processing }
#!   \item{\dots}{ formally required but not used }
#! }
#! \details{
#!   Bit summaries are quite fast because we use a double loop that fixes each word in a processor register.
#!   Furthermore we break out of looping as soon as possible.
#! }
#! \value{
#!   as expected
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{bit}}, \code{\link{all}}, \code{\link{any}}, \code{\link{min}}, \code{\link{max}}, \code{\link{range}}, \code{\link{sum}}, \code{\link{summary}} }
#! \examples{
#!   x <- as.bit(c(TRUE, TRUE))
#!   all(x)
#!   any(x)
#!   min(x)
#!   max(x)
#!   range(x)
#!   sum(x)
#!   summary(x)
#!
#!   x <- as.bitwhich(c(TRUE, TRUE))
#!   all(x)
#!   any(x)
#!   min(x)
#!   max(x)
#!   range(x)
#!   sum(x)
#!   summary(x)
#!
#!  \dontrun{
#!     n <- .Machine$integer.max
#!     x <- !bit(n)
#!     N <- 1000000L  # batchsize
#!     B <- n \%/\% N   # number of batches
#!     R <- n \%\% N    # rest
#!
#!     message("Batched sum (52.5 sec on Centrino duo)")
#!     system.time({
#!       s <- 0L
#!       for (b in 1:B){
#!         s <- s + sum(x[((b-1L)*N+1L):(b*N)])
#!       }
#!       if (R)
#!         s <- s + sum(x[(n-R+1L):n])
#!     })
#!
#!     message("Batched sum saving repeated memory allocation for the return vector
#!     (44.4 sec on Centrino duo)")
#!     system.time({
#!       s <- 0L
#!       l <- logical(N)
#!       for (b in 1:B){
#!         .Call("R_bit_extract", x, ((b-1L)*N+1L):(b*N), l, PACKAGE = "bit")
#!         s <- s + sum(l)
#!       }
#!       if (R)
#!         s <- s + sum(x[(n-R+1L):n])
#!     })
#!
#!     message("C-coded sum (3.1 sec on Centrino duo)")
#!     system.time(sum(x))
#!  }
#! }
#! \keyword{ classes }
#! \keyword{ logic }


sum.bit <- function(x, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(x))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(x))
      stop("illegal range")
  }
  .Call("R_bit_sum", x, range, PACKAGE="bit")
}

all.bit <- function(x, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(x))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(x))
      stop("illegal range")
  }
  .Call("R_bit_all", x, range, PACKAGE="bit")
}

any.bit <- function(x, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(x))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(x))
      stop("illegal range")
  }
  .Call("R_bit_any", x, range, PACKAGE="bit")
}

min.bit <- function(x, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(x))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(x))
      stop("illegal range")
  }
  .Call("R_bit_min", x, range, PACKAGE="bit")
}

max.bit <- function(x, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(x))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(x))
      stop("illegal range")
  }
  .Call("R_bit_max", x, range, PACKAGE="bit")
}

range.bit <- function(x, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(x))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(x))
      stop("illegal range")
  }
  ret <- integer(2)
  ret[1] <- .Call("R_bit_min", x, range, PACKAGE="bit")
  if (is.na(ret[1]))
    ret[2] <- NA
  else
    ret[2] <- .Call("R_bit_max", x, range, PACKAGE="bit")
  ret
}

summary.bit <- function(object, range=NULL, ...){
  if (is.null(range))
    range <- c(1L, length(object))
  else{
    range <- as.integer(range[1:2])
    if (range[1]<1L || range[2]>length(object))
      stop("illegal range")
  }
  s <- sum(object, range=range)
  r <- range(object, range=range)
  c("FALSE"=range[2]-range[1]+1L-s, "TRUE"=s, "Min."=r[1], "Max."=r[2])
}




sum.bitwhich <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'bitwhich'")
  attr(x, "poslength")
}

all.bitwhich <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'bitwhich'")
  attr(x, "poslength") == attr(x, "maxindex")
}

any.bitwhich <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'bitwhich'")
  attr(x, "poslength") > 0L
}

min.bitwhich <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'bitwhich'")
  n <- attr(x, "maxindex")
  p <- attr(x, "poslength")
  if (p==0)
    return(as.integer(NA))
  if (p==n)
    return(n)
  #negative <- p>(n%/%2L)
  negative <- x[1]<0
  if (negative){
    min(as.bit(x))
  }else{
    min(unclass(x))
  }
}

max.bitwhich <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'bitwhich'")
  n <- attr(x, "maxindex")
  p <- attr(x, "poslength")
  if (p==0)
    return(as.integer(NA))
  if (p==n)
    return(n)
  #negative <- p>(n%/%2L)
  negative <- x[1]<0
  if (negative){
    max(as.bit(x))
  }else{
    max(unclass(x))
  }
}

range.bitwhich <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'bitwhich'")
  n <- attr(x, "maxindex")
  p <- attr(x, "poslength")
  if (p==0)
    return(as.integer(NA))
  if (p==n)
    return(n)
  #negative <- p>(n%/%2L)
  negative <- x[1]<0
  if (negative){
    range(as.bit(x))
  }else{
    range(unclass(x))
  }
}

summary.bitwhich <- function(object, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'bitwhich'")
  n <- attr(object, "maxindex")
  p <- attr(object, "poslength")
  if (p==0)
    return(as.integer(NA))
  if (p==n)
    return(n)
  #negative <- p>(n%/%2L)
  negative <- object[1]<0
  if (negative){
    r <- range(as.bit(object))
  }else{
    r <- range(object)
  }
  c("FALSE"=n-p, "TRUE"=p, "Min."=r[1], "Max."=r[2])
}




if (FALSE){
  library(bit)

  # test correctness of max.bit
  for (n in c(0, 1, 2, 31, 32, 33, 63, 64, 65, 95, 96, 97, 127,128,129)){
    for (to1 in seq_len(n)){
      cat("n", n, "to", to1, "\n")
      for (from1 in seq.int(from=1, to=to1, by=1L)){
      x <- bit(n)
      if (!identical(max(x, from=from1, to=to1), as.integer(NA)))
        stop("wrong")
      for (i in seq_len(n)){
        x[i] <- TRUE
        if (!identical(i, max(x, from=from1, to=to1)))
          stop("wrong")
      }
      }
    }
  }


  # test correctness of min.bit
  for (n in c(0, 1, 2, 31, 32, 33, 63, 64, 65, 95, 96, 97, 127,128,129)){
    for (to1 in seq_len(n)){
      cat("n", n, "to", to1, "\n")
      for (from1 in seq.int(from=1, to=to1, by=1L)){
      x <- bit(n)
      if (!identical(min(x, from=from1, to=to1), as.integer(NA)))
        stop("wrong")
      for (i in rev(seq_len(n))){
        x[i] <- TRUE
        if (!identical(i, min(x, from=from1, to=to1)))
          stop("wrong")
      }
      }
    }
  }

}





#! \name{Extract}
#! \alias{[[.bit}
#! \alias{[[<-.bit}
#! \alias{[.bit}
#! \alias{[<-.bit}
#! \title{ Extract or replace part of an bit vector }
#! \description{
#!   Operators acting on bit objects to extract or replace parts.
#! }
#! \usage{
#! \method{[[}{bit}(x, i)
#! \method{[[}{bit}(x, i) <- value
#! \method{[}{bit}(x, i)
#! \method{[}{bit}(x, i) <- value
#! }
#! \arguments{
#!   \item{x}{ a bit object }
#!   \item{i}{ positive integer subscript }
#!   \item{value}{ new logical or integer values }
#! }
#! \details{
#!   Since this package was created for high performance purposes, only positive integer subscripts make sense.
#!   Negative subscripts are converted to positive ones, beware the RAM consumption.
#!   Further subscript classes allowed for '[' and '[<-' are range indices \code{\link{ri}} and \code{\link{bitwhich}}.
#!   The '[' and '[<-' methods don't check whether the subscripts are positive integers in the allowed range.
#! }
#! \value{
#!   The extractors \code{[[} and \code{[} return a logical scalar or vector.
#!   The replacment functions return a bit object.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{bit}}, \code{\link{Extract}} }
#! \examples{
#!   x <- as.bit(c(FALSE, NA, TRUE))
#!   x[] <- c(FALSE, NA, TRUE)
#!   x[1:2]
#!   x[-3]
#!   x[ri(1,2)]
#!   x[as.bitwhich(c(TRUE,TRUE,FALSE))]
#!   x[[1]]
#!   x[] <- TRUE
#!   x[1:2] <- FALSE
#!   x[[1]] <- TRUE
#! }
#! \keyword{ classes }
#! \keyword{ logic }



"[[.bit" <- function(x, i){
  if (length(i)!=1)
    stop("subscript length not 1")
  if (is.numeric(i)){
    i <- as.integer(i)
    if (is.na(i) || i<1L || i>length(x))
      stop("subscript must be positive integer (or double) within length")
    ret <- logical(1L)
    attr(ret, "vmode") <- "boolean"
    .Call("R_bit_extract", x, i, ret, PACKAGE="bit")
  }else
    stop("subscript must be positive integer (or double) within length")
}


"[[<-.bit" <- function(x, i, value){
  if (length(i)!=1)
    stop("subscript length not 1")
  if (length(value)!=1)
    stop("value length not 1")
  if (is.numeric(i)){
    i <- as.integer(i)
    if (is.na(i) || i<1L || i>length(x))
      stop("subscript must be positive integer (or double) within length")
    value2 <- as.logical(value)
    .Call("R_bit_replace", x, i, value2, FALSE, PACKAGE="bit")
  }else
    stop("subscript must be positive integer (or double) within length")
}

if (FALSE){
  library(ff)
  library(bit)
  a <- bit(100)
  a[1] <- T
  a[100] <- T
  a[]
  a[99:100]
  a[range=c(99,100)]

  a[range=c(1,100)] <- TRUE
  a
  a[range=c(1,100)] <- FALSE
  a


}



"[.bit" <- function(x, i){
  if ( missing(i) ){
    len <- length(x)
    ret <- logical(len)
    .Call("R_bit_get", x, ret, range=c(1L, len), PACKAGE="bit")
  }else if(is.numeric(i)){
    if (inherits(i, "ri")){
      if (i[1]<1L || i[2]>length(x))
        stop("illegal range index 'ri'")
      ret <- logical(i[2]-i[1]+1L)
      .Call("R_bit_get", x, ret, range=i, PACKAGE="bit")
    }else if (inherits(i, "bitwhich")){
			i <- as.which(i)
      n <- length(i)
			ret <- logical(n)
			if (n)
        .Call("R_bit_extract", x, i, ret, PACKAGE="bit")
		}else{
      i <- as.integer(i)
      if (length(i)){
        if (i[1]<0)
          i <- (as.integer(seq_along(x)))[i]
				ret <- logical(length(i))
        .Call("R_bit_extract", x, i, ret, PACKAGE="bit")
      }else{
        ret <- logical()
      }
    }
  }else if(is.logical(i)){
    if (length(i)!=1 || is.na(i)){
      stop("only TRUE or FALSE allowed")
    }else{
      if (i){
        len <- length(x)
        ret <- logical(len)
        .Call("R_bit_get", x, ret, range=c(1L, len), PACKAGE="bit")
      }else{
        ret <- logical()
      }
    }
  }else
      stop("subscript must be integer (or double) or bitwhich")

  attr(ret, "vmode") <- "boolean"
  ret
}


"[<-.bit" <- function(x, i, value){
  if ( missing(i) ){
    len <- length(x)
    if (length(value)==len){
      value2 <- as.logical(value)
    }else{
      value2 <- logical(len)
      value2[] <- value
    }
    .Call("R_bit_set", x, value2, range=c(1L, len), PACKAGE="bit")
  }else if(is.numeric(i)){
    if (inherits(i, "ri")){
      if (i[1]<1L || i[2]>length(x))
        stop("illegal range index 'ri'")
      n <- i[2] - i[1] + 1L
      if (length(value)==n){
        value2 <- as.logical(value)
      }else{
        value2 <- logical(n)
        value2[] <- value
      }
      .Call("R_bit_set", x, value2, range=i, PACKAGE="bit")
    }else{
			if (inherits(i, "bitwhich")){
				i <- as.which(i)
				n <- length(i)
			}else{
				i <- as.integer(i)
				n <- length(i)
				if (n && i[1]<0){
					i <- (as.integer(seq_along(x)))[i]
					n <- length(i)
				}
			} 
      if (length(value)==n){
        value2 <- as.logical(value)
      }else{
        value2 <- logical(n)
        value2[] <- value
      }
      .Call("R_bit_replace", x, i, value2, PACKAGE="bit")
    }
  }else if (is.logical(i)){
    if (length(i)!=1 || is.na(i)){
      stop("only TRUE or FALSE allowed")
    }else{
      if (i){
        len <- length(x)
        if (length(value)==len){
          value2 <- as.logical(value)
        }else{
          value2 <- logical(len)
          value2[] <- value
        }
        .Call("R_bit_set", x, value2, range=c(1L, len), PACKAGE="bit")
      }else{
        x
      }
    }
  }else
      stop("subscript must be integer (or double) within length")
}



#! \name{ri}
#! \alias{ri}
#! \alias{print.ri}
#! \title{ Range index }
#! \description{
#!   A range index can be used to extract or replace a continuous ascending part of the data
#! }
#! \usage{
#! ri(from, to = NULL, maxindex=NA)
#! \method{print}{ri}(x, \dots)
#!
#! }
#! \arguments{
#!   \item{from}{ first position }
#!   \item{to}{ last posistion }
#!   \item{x}{ an object of class 'ri' }
#!   \item{maxindex}{ the maximal length of the object-to-be-subscripted (if known) }
#!   \item{\dots}{ further arguments }
#! }
#! \value{
#!   A two element integer vector with class 'ri'
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[ff]{as.hi.ri}} }
#! \examples{
#!  bit(12)[ri(1,6)]
#! }
#! \keyword{ classes }
#! \keyword{ logic }

ri <- function(from, to=NULL, maxindex=NA){
  if (is.null(to)){
    x <- as.integer(c(from, maxindex))
  }else{
    x <- as.integer(c(from, to, maxindex))
  }
  maxindex = maxindex
  if (length(x)!=3 )
    stop("range must have exactly three elements")
  if (x[[1]]<1L)
    stop("range must at least select one element")
  if (x[[1]]>x[[2]])
    stop("lower bound must be smaller or equal than upper bound")
  if (!is.na(x[[3]]) && x[[2]]>x[[3]])
    stop("lower and upper bound must be smaller or equal to maxindex")
  oldClass(x) <- "ri"
  x
}


print.ri <- function(x, ...)
  cat("range index (ri) from", x[[1]], "to", x[[2]], "maxindex",  x[[3]], "\n")


length.ri <- function(x)
x[[3]]

all.ri <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'ri'")
  x[[1]]<=1L && x[[2]]>=x[[3]]
}

any.ri <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'ri'")
  TRUE
}

min.ri <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'ri'")
  x[[1]]
}

max.ri <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'ri'")
  x[[2]]
}

range.ri <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'ri'")
  x[1:2]
}

sum.ri <- function(x, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'ri'")
  x[[2]] - x[[1]] + 1L
}

summary.ri <- function(object, ...){
  if (any(names(match.call(expand.dots = TRUE))=="range"))
    stop("parameter 'range' allowed only for 'bit' but not for 'ri'")
  s <- object[[2]] - object[[1]] + 1L
   c(`FALSE` = object[[3]] - s, `TRUE` = s, Min. = object[[1]], Max. = object[[2]])
}


#! \name{physical}
#! \alias{physical}
#! \alias{physical<-}
#! \alias{virtual}
#! \alias{virtual<-}
#! \alias{physical.default}
#! \alias{physical<-.default}
#! \alias{virtual.default}
#! \alias{virtual<-.default}
#! \alias{print.physical}
#! \alias{print.virtual}
#! \title{ Physical and virtual attributes }
#! \description{
#!   Compatibility functions (to package ff) for getting and setting physical and virtual attributes.
#! }
#! \usage{
#! physical(x)
#! virtual(x)
#! physical(x) <- value
#! virtual(x) <- value
#! \method{physical}{default}(x)
#! \method{virtual}{default}(x)
#! \method{physical}{default}(x) <- value
#! \method{virtual}{default}(x) <- value
#! \method{print}{physical}(x, \dots)
#! \method{print}{virtual}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ a ff or ram object }
#!   \item{value}{ a list with named elements }
#!   \item{\dots}{ further arguments }
#! }
#! \details{
#!   ff objects have physical and virtual attributes, which have different copying semantics:
#!   physical attributes are shared between copies of ff objects while virtual attributes might differ between copies.
#!   \code{\link[ff]{as.ram}} will retain some physical and virtual atrributes in the ram clone,
#!   such that \code{\link[ff]{as.ff}} can restore an ff object with the same attributes.
#! }
#! \value{
#!   \command{physical} and \command{virtual} returns a list with named elements
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{
#!  \code{\link[ff]{physical.ff}}, \code{\link[ff]{physical.ffdf}}
#! }
#! \examples{
#!   physical(bit(12))
#!   virtual(bit(12))
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }


# this version without vmode() will be overwritte by the version in package ff
physical.default <- function(x){
  p <- attributes(attr(x, "physical"))
  p <- p[is.na(match(names(p), "class"))]
  p
}
"physical<-.default" <- function(x, value){
  attributes(attr(x, "physical")) <- c(value, list(class="physical"))
  x
}


virtual.default <- function(x){
  v <- attributes(attr(x, "virtual"))
  v[is.na(match(names(v), "class"))]
}
"virtual<-.default" <- function(x, value){
  attributes(attr(x, "virtual")) <- c(value, list(class="virtual"))
  x
}


print.physical <- function(x, ...){
  cat("(hidden, use physical(x) to access the physical attributes and vmode(x) for accessing vmode)\n")
  invisible()
}

print.virtual <- function(x, ...){
  cat("(hidden, use virtual(x) to access the virtual attributes)\n")
  invisible()
}


# not exported - just here to avoid cross calling the dll from ff
R_bit_as_hi <- function(x, range, offset)
.Call("R_bit_as_hi", x, range, offset, PACKAGE="bit")



#! \name{regtest.bit}
#! \alias{regtest.bit}
#! \title{ Regressiontests for bit }
#! \description{
#!   Test package bit for correctness
#! }
#! \usage{
#! regtest.bit(N = 100)
#! }
#! \arguments{
#!   \item{N}{ number of random test runs }
#! }
#! \details{
#!   random data of random length are generated and correctness of package functions tested on these
#! }
#! \value{
#!   a vector of class 'logical' or 'integer'
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{bit}}, \code{\link{as.bit}}, \code{\link{as.logical}}, \code{\link{as.integer}}, \code{\link{which}} }
#! \examples{
#!   if (regtest.bit()){
#!     message("regtest.bit is OK")
#!   }else{
#!     message("regtest.bit failed")
#!   }
#!
#!   \dontrun{
#!     regtest.bit(10000)
#!   }
#! }
#! \keyword{ classes }
#! \keyword{ logic }

regtest.bit <- function(
    N = 100  # number of repetitions for random regression tests
)
{
  OK <- TRUE
  pool <- c(FALSE, TRUE)

  for (i in 1:N){
    n <- sample(1:(2*.BITS), 1)
    l <- sample(pool, n, TRUE)
    # check direct coercion
    b <- as.bit(l)
    l2 <- as.logical(b)
    if (!identical(l,l2)){
      message("\nregression test difference between logical")
      print(l)
      message("and as.logical(as.bit(logical))")
      print(l2)
      OK <- FALSE
    }
    # summary functions with logical return
    s <- c(all=all(l), any=any(l))
    s2 <- c(all=all(b), any=any(b))
    if (!identical(s,s2)){
      message("\nregression test difference between logical summaries")
      print(s)
      message("and bit summaries")
      print(s2)
      OK <- FALSE
    }
    # summary functions with integer return
    if (any(l)){
        s <- c(min=min(as.which(l)), max=max(as.which(l)), range=range(as.which(l)), sum=sum(l), summary=c("FALSE"=length(l)-sum(l), "TRUE"=sum(l), "Min."=min(as.which(l)), "Max."=max(as.which(l))))
    }else{
      s <- c( min=as.integer(NA), max=as.integer(NA), range=c(as.integer(NA), as.integer(NA)), sum=sum(l), summary=c("FALSE"=length(l)-sum(l), "TRUE"=sum(l), "Min."=as.integer(NA), "Max."=as.integer(NA)) )
    }
    s2 <- c(min=min(b), max=max(b), range=range(b), sum=sum(b), summary=summary(b))
    if (!identical(s,s2)){
      message("\nregression test difference between logical summaries")
      print(s)
      message("and bit summaries")
      print(s2)
      OK <- FALSE
    }
    # check positive whichs
    w <- as.which(l)
    w2 <- as.which(as.bit.which(w, n))
    if (!identical(w,w2)){
      message("\nregression test difference between which")
      print(w)
      message("and as.which(as.bit.which(which))")
      print(w2)
      OK <- FALSE
    }
    # check automatic whichs (pos or neg whatever shorter)
    s <- sum(l)
    if (s==0){
      w <- FALSE
    }else if (s==n){
      w <- TRUE
    }else if (s>(n%/%2L)){
      w <- -rev(which(!l))
    }else{
      w <- which(l)
    }
    w2 <- as.vector(as.bitwhich(as.bit(l)))
    if (!identical(w,w2)){
      message("\nregression test difference between which")
      print(w)
      message("and as.which(as.bit.which(which))")
      print(w2)
      OK <- FALSE
    }
    # check boolean operators
    l2 <- sample(c(FALSE, TRUE), n, TRUE)
    b2 <- as.bit(l2)
    ops <- c(
      NOT = identical(!l, as.logical(!b))
    , AND = identical(l&l2, as.logical(b&b2))
    , OR = identical(l|l2, as.logical(b|b2))
    , XOR = identical(xor(l,l2), as.logical(xor(b,b2)))
    , NEQ = identical(l!=l2, as.logical(b!=b2))
    , EQ = identical(l==l2, as.logical(b==b2))
    )
    if (!all(ops)){
      message("\nbit differs for boolean operators(s)")
      print(ops)
      print(cbind(l=l, l2=l))
      OK <- FALSE
    }
    w <- as.bitwhich(l)
    w2 <- as.bitwhich(l2)
    ops <- c(
      NOT = identical(!l, as.logical(!w))
    , AND = identical(l&l2, as.logical(w&w2))
    , OR = identical(l|l2, as.logical(w|w2))
    , XOR = identical(xor(l,l2), as.logical(xor(w,w2)))
    , NEQ = identical(l!=l2, as.logical(w!=w2))
    , EQ = identical(l==l2, as.logical(w==w2))
    )
    if (!all(ops)){
      message("\nbitwhich differs for boolean operators(s)")
      print(ops)
      print(cbind(l=l, l2=l))
      OK <- FALSE
    }
    rm(l2,b2,w2)
    # check extractors
    n2 <- sample(1:n, 1)
    j <- sample(1:n, n2)
    if (!identical(l[j], unattr(b[j]))){
      message("\nregression test difference when extracting")
      OK <- FALSE
    }
    # check replacement (index)
    new <- sample(pool, n2, TRUE)
    l[j] <- new
    b[j] <- new
    if (!identical(l, unattr(b[]))){
      message("\nregression test difference when replacing with index")
      OK <- FALSE
    }
    # check replacement (recycle)
    if (n%%2){
      new <- sample(pool, 1)
      l[] <- new
      b[] <- new
    }else{
      l[] <- pool
      b[] <- pool
    }
    if (!identical(l, as.logical(b))){
      message("\nregression test difference when replacing with recylcling")
      OK <- FALSE
    }
  }

  l0 <- c(FALSE, FALSE, FALSE)
  l1 <- c(FALSE, FALSE, TRUE)
  l2 <- c(FALSE, TRUE, TRUE)
  l3 <- c(TRUE, TRUE, TRUE)

  bw0 <- as.bitwhich(l0)
  bw1 <- as.bitwhich(l1)
  bw2 <- as.bitwhich(l2)
  bw3 <- as.bitwhich(l3)

  OK <- OK && identical(l0, as.logical(bw0))
  OK <- OK && identical(l1, as.logical(bw1))
  OK <- OK && identical(l2, as.logical(bw2))
  OK <- OK && identical(l3, as.logical(bw3))

  OK <- OK && identical(l0 & l0, as.logical(bw0 & bw0))
  OK <- OK && identical(l0 & l1, as.logical(bw0 & bw1))
  OK <- OK && identical(l0 & l2, as.logical(bw0 & bw2))
  OK <- OK && identical(l0 & l3, as.logical(bw0 & bw3))

  OK <- OK && identical(l1 & l0, as.logical(bw1 & bw0))
  OK <- OK && identical(l1 & l1, as.logical(bw1 & bw1))
  OK <- OK && identical(l1 & l2, as.logical(bw1 & bw2))
  OK <- OK && identical(l1 & l3, as.logical(bw1 & bw3))

  OK <- OK && identical(l2 & l0, as.logical(bw2 & bw0))
  OK <- OK && identical(l2 & l1, as.logical(bw2 & bw1))
  OK <- OK && identical(l2 & l2, as.logical(bw2 & bw2))
  OK <- OK && identical(l2 & l3, as.logical(bw2 & bw3))

  OK <- OK && identical(l3 & l0, as.logical(bw3 & bw0))
  OK <- OK && identical(l3 & l1, as.logical(bw3 & bw1))
  OK <- OK && identical(l3 & l2, as.logical(bw3 & bw2))
  OK <- OK && identical(l3 & l3, as.logical(bw3 & bw3))


  OK <- OK && identical(l0 | l0, as.logical(bw0 | bw0))
  OK <- OK && identical(l0 | l1, as.logical(bw0 | bw1))
  OK <- OK && identical(l0 | l2, as.logical(bw0 | bw2))
  OK <- OK && identical(l0 | l3, as.logical(bw0 | bw3))

  OK <- OK && identical(l1 | l0, as.logical(bw1 | bw0))
  OK <- OK && identical(l1 | l1, as.logical(bw1 | bw1))
  OK <- OK && identical(l1 | l2, as.logical(bw1 | bw2))
  OK <- OK && identical(l1 | l3, as.logical(bw1 | bw3))

  OK <- OK && identical(l2 | l0, as.logical(bw2 | bw0))
  OK <- OK && identical(l2 | l1, as.logical(bw2 | bw1))
  OK <- OK && identical(l2 | l2, as.logical(bw2 | bw2))
  OK <- OK && identical(l2 | l3, as.logical(bw2 | bw3))

  OK <- OK && identical(l3 | l0, as.logical(bw3 | bw0))
  OK <- OK && identical(l3 | l1, as.logical(bw3 | bw1))
  OK <- OK && identical(l3 | l2, as.logical(bw3 | bw2))
  OK <- OK && identical(l3 | l3, as.logical(bw3 | bw3))


  OK <- OK && identical(xor(l0,l0), as.logical(xor(bw0,bw0)))
  OK <- OK && identical(xor(l0,l1), as.logical(xor(bw0,bw1)))
  OK <- OK && identical(xor(l0,l2), as.logical(xor(bw0,bw2)))
  OK <- OK && identical(xor(l0,l3), as.logical(xor(bw0,bw3)))

  OK <- OK && identical(xor(l1,l0), as.logical(xor(bw1,bw0)))
  OK <- OK && identical(xor(l1,l1), as.logical(xor(bw1,bw1)))
  OK <- OK && identical(xor(l1,l2), as.logical(xor(bw1,bw2)))
  OK <- OK && identical(xor(l1,l3), as.logical(xor(bw1,bw3)))

  OK <- OK && identical(xor(l2,l0), as.logical(xor(bw2,bw0)))
  OK <- OK && identical(xor(l2,l1), as.logical(xor(bw2,bw1)))
  OK <- OK && identical(xor(l2,l2), as.logical(xor(bw2,bw2)))
  OK <- OK && identical(xor(l2,l3), as.logical(xor(bw2,bw3)))

  OK <- OK && identical(xor(l3,l0), as.logical(xor(bw3,bw0)))
  OK <- OK && identical(xor(l3,l1), as.logical(xor(bw3,bw1)))
  OK <- OK && identical(xor(l3,l2), as.logical(xor(bw3,bw2)))
  OK <- OK && identical(xor(l3,l3), as.logical(xor(bw3,bw3)))


  OK <- OK && identical(c(l0,l0), as.logical(c(bw0,bw0)))
  OK <- OK && identical(c(l0,l1), as.logical(c(bw0,bw1)))
  OK <- OK && identical(c(l0,l2), as.logical(c(bw0,bw2)))
  OK <- OK && identical(c(l0,l3), as.logical(c(bw0,bw3)))

  OK <- OK && identical(c(l1,l0), as.logical(c(bw1,bw0)))
  OK <- OK && identical(c(l1,l1), as.logical(c(bw1,bw1)))
  OK <- OK && identical(c(l1,l2), as.logical(c(bw1,bw2)))
  OK <- OK && identical(c(l1,l3), as.logical(c(bw1,bw3)))

  OK <- OK && identical(c(l2,l0), as.logical(c(bw2,bw0)))
  OK <- OK && identical(c(l2,l1), as.logical(c(bw2,bw1)))
  OK <- OK && identical(c(l2,l2), as.logical(c(bw2,bw2)))
  OK <- OK && identical(c(l2,l3), as.logical(c(bw2,bw3)))

  OK <- OK && identical(c(l3,l0), as.logical(c(bw3,bw0)))
  OK <- OK && identical(c(l3,l1), as.logical(c(bw3,bw1)))
  OK <- OK && identical(c(l3,l2), as.logical(c(bw3,bw2)))
  OK <- OK && identical(c(l3,l3), as.logical(c(bw3,bw3)))

  N <- 2L*.BITS
  l <- logical(N)
  b <- bit(N)
  for (i in 1:N){
    l[i] <- TRUE
    b[i] <- TRUE
    if (!identical(l,as.logical(b))){
      message("\nregression test difference when replacing at position", i, "")
      OK <- FALSE
    }
  }

  OK
}
