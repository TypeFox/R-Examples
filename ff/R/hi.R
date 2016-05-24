# hybrid index packing with run length encoding
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-08-24
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/hi.R")

# -- creating a hi object directly ----------------------------------------------------------------------


#! \name{hi}
#! \alias{hi}
#! \alias{print.hi}
#! \alias{str.hi}
#! \title{ Hybrid index class }
#! \description{
#!   Class for hybrid index representation, plain and rle-packed
#! }
#! \usage{
#! hi(from, to, by = 1L, maxindex = NA, vw = NULL, pack = TRUE, NAs = NULL)
#! \method{print}{hi}(x, \dots)
#! \method{str}{hi}(object, nest.lev=0, \dots)
#! }
#! \arguments{
#!   \item{from}{ integer vector of lower sequence bounds }
#!   \item{to}{ integer vector of upper sequence bounds }
#!   \item{by}{ integer of stepsizes }
#!   \item{maxindex}{ maximum indep position (needed for negative indices) }
#!   \item{vw}{ virtual window information, see \code{\link{vw}} }
#!   \item{pack}{ FALSE to suppress rle-packing }
#!   \item{NAs}{ a vector of NA positions (not yet used) }
#!   \item{x}{ an object of class 'hi' to be printed }
#!   \item{object}{ an object of class 'hi' to be str'ed }
#!   \item{nest.lev}{ current nesting level in the recursive calls to str }
#!   \item{\dots}{ further arguments passed to the next method }
#! }
#! \details{
#!   Class \code{hi} will represent index data either as a plain positive or negative index vector or as an rle-packed version thereof.
#!   The current implementation switches from plain index positions \code{i} to rle-packed storage of \code{diff(i)} as soon as the compression ratio is 3 or higher.
#!   Note that sequences shorter than 2 must never be packed (could cause C-side crash).
#!   Furthermore hybrid indices are guaranteed to be sorted ascending, which helps \code{\link{ff}s} access method avoiding to swap repeatedly over the same memory pages (or file positions).
#! }
#! \value{
#!   A list of class 'hi' with components
#!   \item{ x      }{ directly accessed by the C-code: the sorted index as returned by \code{\link[bit]{rlepack}} }
#!   \item{ ix     }{ NULL or positions to restore original order }
#!   \item{ re     }{ logical scalar indicating if sequence was reversed from descending to ascending (in this case \code{is.null(ix)}) }
#!   \item{ minindex  }{ directly accessed by the C-code: represents the lowest positive subscript to be enumerated in case of negative subscripts }
#!   \item{ maxindex  }{ directly accessed by the C-code: represents the highest positive subscript to be enumerated in case of negative subscripts }
#!   \item{ length    }{ number of subscripts, whether negative or positive, not the number of selected elements }
#!   \item{ dim      }{ NULL or dim -- used by \code{\link{as.matrix.hi}} }
#!   \item{ dimorder  }{ NULL or \code{\link{dimorder}} }
#!   \item{ symmetric }{ logical scalar indicating whether we have a symmetric matrix }
#!   \item{ fixdiag   }{ logical scalar indicating whether we have a fixed diagonal (can only be true for symmetric matrices) }
#!   \item{ vw     }{ virtual window information \code{\link{vw}} }
#!   \item{ NAs      }{ NULL or NA positions as returned by \code{\link[bit]{rlepack}} }
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ \command{hi} defines the class structure, however usually \code{\link{as.hi}} is used to acturally Hybrid Index Preprocessing for \code{\link{ff}} }
#! \seealso{ \code{\link{as.hi}} for coercion, \code{\link[bit]{rlepack}}, \code{\link[bit]{intrle}}, \code{\link{maxindex}}, \code{\link{poslength}} }
#! \examples{
#!   hi(c(1, 11, 29), c(9, 19, 21), c(1,1,-2))
#!   as.integer(hi(c(1, 11, 29), c(9, 19, 21), c(1,1,-2)))
#! }
#! \keyword{ IO }
#! \keyword{ data }

hi <- function (from, to, by = 1L, maxindex = NA, vw=NULL, pack = TRUE, NAs = NULL)
{
    minindex <- 1L
    maxindex <- as.integer(maxindex)
    if (is.null(vw)){
      vw.convert <- FALSE
    }else{
      if (is.matrix(vw))
        stop("matrix vw not allowed in hi, use as.hi")
      storage.mode(vw) <- "integer"
      vw.convert <- TRUE
    }
    nspec <- length(from)
    if (nspec > 0) {
        from <- as.integer(from)
        to <- rep(as.integer(to), length = nspec)
        by <- rep(as.integer(by), length = nspec)
        d <- to - from
        N <- d%/%by
        if (any(d != 0 & sign(d) != sign(by)) || any(N * by != d))
            stop("illegal input to hi")
        l <- as.vector(rbind(rep(1L, nspec), N))[-1]
        v <- as.vector(rbind(c(0L, from[-1] - to[-nspec]), by))[-1]
        v <- v[l > 0]
        l <- l[l > 0]
        from <- from[1]
        to <- to[nspec]
        nl <- length(l)
        r <- list(lengths = l, values = v)
        n <- sum(r$lengths) + 1L
        tab <- tabulate(sign(r$values) + 2, 3)
        s <- !tab[1] || !tab[3]
        if (s) {  # sorted
            #if (nl) {
            #    if (pack)
            #      pack <- 2 * nl < n
            #}else
            #  pack <- FALSE
            #if (pack){
                class(r) <- "rle"
            #}else{
            #    r <- as.integer(cumsum(c(from, rep(r$values, r$lengths))))
            #}
            x <- list(first = from, dat = r, last = to)
            ix <- NULL
            re <- tab[1] > 0
            if (re)
                x <- rev.rlepack(x)
        }else{
            re <- FALSE
            x <- as.integer(cumsum(c(from, rep(r$values, r$lengths))))
            x <- sort.int(x, index.return = TRUE, method = "quick")
            ix <- x$ix
            x <- rlepack(x$x, pack = pack)
            #ix <- 1:n
            #radixorder(x, ix)
            #x <- rlepack(x[ix], pack = pack)
        }

        x <- unique.rlepack(x)

        # this ifelse section copied 1L to as.hi.call
        if (x$last < 0) {
          if (is.na(maxindex))
              stop("maxindex is required with negative subscripts")
          if ( -x$first > maxindex )
              stop("negative subscripts out of range")

          re <- FALSE
          ix <- NULL

          if (vw.convert){
            x$first <- x$first - vw[1]
            x$last <- x$last - vw[1]
            if (inherits(x$dat, "rle")){
              n <- sum(x$dat$lengths) + 1L
            }else{
              x$dat <- x$dat - vw[1]
              n <- length(x$dat)
            }
          }else{
            if (inherits(x$dat, "rle")){
              n <- sum(x$dat$lengths) + 1L
            }else{
              n <- length(x$dat)
            }
          }
        }else if (x$first > 0){
          if (!is.na(maxindex) && x$last > maxindex )
              stop("positive subscripts out of range")
          if (vw.convert){
            x$first <- vw[1] + x$first
            x$last <- vw[1] + x$last
            if (inherits(x$dat, "rle")){
              n <- sum(x$dat$lengths) + 1L
            }else{
              x$dat <- vw[1] + x$dat
              n <- length(x$dat)
            }
          }else{
            if (inherits(x$dat, "rle")){
              n <- sum(x$dat$lengths) + 1L
            }else{
              n <- length(x$dat)
            }
          }
        }else{
            stop("0s and mixed positive/negative subscripts not allowed")
        }
    }else{
        x <- list(first = as.integer(NA), dat = integer(), last = as.integer(NA))
        re <- FALSE
        ix <- NULL
        n <- 0L
        minindex <- 1L
        maxindex <- as.integer(maxindex)
    }
    if (!is.null(NAs))
        NAs <- rlepack(as.integer(NAs), pack = pack)

    if (!is.null(vw)){
      minindex <- vw[1] + 1L
      maxindex <- vw[1] + vw[2]
    }
    ret <- list(
      x = x                # directly accessed by the C-code: hybrid index, i.e. either raw or rle
    , ix = ix              # NULL or positions for re-ordering
    , re = re              # logical indicating whether sequence was reversed from descending to ascending
    , minindex = minindex  # directly accessed by the C-code: represents the lowest positive subscript to be enumerated in case of negative subscripts
    , maxindex = maxindex  # directly accessed by the C-code: represents the highest positive subscript to be enumerated in case of negative subscripts
    , length = n           # number of subscripts, whether negative or positive
    , dim = NULL           # NULL or dim
    , dimorder = NULL      # NULL or dimorder
    , symmetric = FALSE    # logical indicating whether we have a symmetric matrix
    , fixdiag = NULL       # logical indicating whether we have a fixed diagonal (can only be true for symmetric matrices)
    , vw = vw            # NULL or OffsetWindowRest definition
    , NAs = NAs            # NULL or positive positions of NAs
    )
    class(ret) <- "hi"
    ret
}


print.hi <- function(x, ...){
  cat("hybrid index (hi) from ", x$x$first, " to ", x$x$last, " over ", if (inherits(x$x$dat, "rle")) "<rle position diffs>" else "<plain positions>", " re=", x$re, " ix=", if(is.null(x$ix)) "NULL" else "<reverse sort info>", "\n", sep="")
  cat("minindex=", x$minindex, " maxindex=", x$maxindex, " length=", x$length, " poslength=", poslength(x), "\n", sep="")
  if (!is.null(x$dim)){
    cat("dim=c(", paste(x$dim, collapse=","), "), dimorder=c(", paste(x$dimorder, collapse=","), ")\n", sep="")
  }
  if (!is.null(x$vw)){
    cat("vw=")
    print(x$vw, ...)
  }
  invisible()
}


str.hi <- function(object, nest.lev=0, ...){
  nest.str <- paste(rep(" ..", nest.lev), collapse="")
  str(unclass(object), nest.lev=nest.lev, ...)
  cat(nest.str, ' - attr(*, "class") = ', sep="")
  str(class(object), nest.lev=nest.lev, ...)
}



# -- coerce to hi object ----------------------------------------------------------------------


#! \name{hiparse}
#! \alias{hiparse}
#! \title{ Hybrid Index, parsing }
#! \description{
#!   \command{hiparse} implements the parsing done in Hybrid Index Preprocessing in order to avoid RAM for expanding index expressions.
#!   \emph{Not to be called directly}
#! }
#! \usage{
#! hiparse(x, envir, first = as.integer(NA), last = as.integer(NA))
#! }
#! \arguments{
#!   \item{x}{ an index expression, precisely: \code{\link{call}} }
#!   \item{envir}{ the environemtn in which to evaluate components of the index expression }
#!   \item{first}{ first index position found so far }
#!   \item{last}{ last index position found so far }
#! }
#! \details{
#!   This primitive parser recognises the following tokens: numbers like 1, symbols like x, the colon sequence operator \code{\link{:}} and the concat operator \code{\link{c}}.
#!   \code{hiparse} will \code{\link{Recall}} until the index expression is parsed or an unknown token is found.
#!   If an unknown token is found, \code{hiparse} evluates it, inspects it and either accepts it or throws an error, catched by \code{\link{as.hi.call}},
#!   which falls back to evaluating the index expression and dispatching (again) an appropriate \code{\link{as.hi}} method.
#!   Reasons for suspending the parsing: if the inspected token is of class 'hi', 'ri', 'bit', 'bitwhich', 'is.logical', 'is.character', 'is.matrix' or has length>16.
#! }
#! \value{
#!   undefined (and redefined as needed by \code{\link{as.hi.call}})
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{hi}}, \code{\link{as.hi.call}} }
#! \keyword{ IO }
#! \keyword{ data }


hiparse <- function(x, envir, first=as.integer(NA), last=as.integer(NA)){
  if (length(x)>1){
    if (x[[1]]=='c'){
      values <- integer()
      lengths <- integer()
      n <- length(x)
      i <- 1
      while(i<n){
        i <- i + 1
        r <- Recall(x[[i]], envir, first=first, last=last)
        first <- r$first
        last <- r$last
        values <- c(values, r$values)
        lengths <- c(lengths, r$lengths)
      }
      return(list(first=first, lengths=lengths, values=values, last=last))
    }else if (x[[1]]==':'){
      from <- eval(x[[2]], envir=envir)
      to <- eval(x[[3]], envir=envir)
      if (is.logical(from) || is.logical(to))
        stop("as.hi.default:hiparse logicals encountered")
      if (length(from)!=1 || length(to)!=1)
        stop("as.hi.default:hiparse: arguments of : have length!=1")
      from <- as.integer(from)
      to <- as.integer(to)
      if ( is.na(from) || is.na(to) || from==0 || to==0 )
        stop("as.hi.default:hiparse NAs or 0s encountered")
      if (is.na(first))
        first <- from
      if (is.na(last)){
        if (from>to)
          return(list(first=first, lengths=from-to, values=as.integer(-1), last=to))
        else
          return(list(first=first, lengths=to-from, values=as.integer(1), last=to))
      }else{
        if (from>to)
          return(list(first=first, lengths=c(as.integer(1), from-to), values=c(from-last, as.integer(-1)), last=to))
        else
          return(list(first=first, lengths=c(as.integer(1), to-from), values=c(from-last, as.integer(1)), last=to))
      }
    }
  }
  x <- eval(x, envir=envir)
  if (inherits(x,"hi"))
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found hi")
  if (inherits(x,"ri"))
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found ri")
  if (inherits(x,"bit"))
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found bit")
  if (inherits(x,"bitwhich"))
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found bitwhich")
  if (is.logical(x))
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found logical")
  if (is.character(x))
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found character")
  if (is.matrix(x))
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found matrix")
  n <- length(x)
  if (n>16)
    stop("DEBUGINFO visible when try(..., silent=FALSE) in as.hi.call: as.hi.default:hiparse found length>16")
  if (n){
    x <- as.integer(x)
    if (is.na(first))
      first <- x[1]
    if (is.na(last)){
      r <- rle(diff(x))           # using standard rle here because x is small and we always need a result
    }else{
      r <- rle(diff(c(last, x)))  # using standard rle here because x is small and we always need a result
    }
    if (is.na(intisasc(r$values)))
      stop("as.hi.default:hiparse found NAs")
    last <- x[n]
    return(list(first=first, lengths=r$lengths, values=r$values, last=last))
  }else{
    return(list(first=first, lengths=integer(), values=integer(), last=last))
  }
}

#! \name{as.hi}
#! \alias{as.hi}
#! \alias{as.hi.hi}
#! \alias{as.hi.ri}
#! \alias{as.hi.bit}
#! \alias{as.hi.bitwhich}
#! \alias{as.hi.call}
#! \alias{as.hi.name}
#! \alias{as.hi.(}
#! \alias{as.hi.integer}
#! \alias{as.hi.which}
#! \alias{as.hi.double}
#! \alias{as.hi.logical}
#! \alias{as.hi.character}
#! \alias{as.hi.matrix}
#! \title{ Hybrid Index, coercion to }
#! \description{
#!   The generic \command{as.hi} and its methods are the main (internal) means for preprocessing index information into the hybrid index class \code{\link{hi}}.
#!   Usually \command{as.hi} is called transparently from \code{\link{[.ff}}. However, you can explicitely do the index-preprocessing,
#!   store the Hybrid Index \code{\link{hi}}, and use the \code{hi} for subscripting.
#! }
#! \usage{
#! as.hi(x, \dots)
#! \method{as.hi}{hi}(x, \dots)
#! \method{as.hi}{ri}(x, maxindex = length(x), \dots)
#! \method{as.hi}{bit}(x, range = NULL, maxindex = length(x), vw = NULL
#! , dim = NULL, dimorder = NULL, pack = TRUE, \dots)
#! \method{as.hi}{bitwhich}(x, maxindex = length(x), pack = FALSE, \dots)
#! \method{as.hi}{call}(x, maxindex = NA, dim = NULL, dimorder = NULL, vw = NULL
#! , vw.convert = TRUE, pack = TRUE, envir = parent.frame(), \dots)
#! \method{as.hi}{name}(x, envir = parent.frame(), \dots)
#! %\method{as.hi}{(}(x, envir = parent.frame(), \dots)
#! \method{as.hi}{integer}(x, maxindex = NA, dim = NULL, dimorder = NULL
#! , symmetric = FALSE, fixdiag = NULL, vw = NULL, vw.convert = TRUE
#! , dimorder.convert  = TRUE, pack = TRUE, NAs = NULL, \dots)
#! \method{as.hi}{which}(x, \dots)
#! \method{as.hi}{double}(x, \dots)
#! \method{as.hi}{logical}(x, maxindex = NA, dim = NULL, vw = NULL, pack = TRUE, \dots)
#! \method{as.hi}{character}(x, names, vw = NULL, vw.convert = TRUE, \dots)
#! \method{as.hi}{matrix}(x, dim, dimorder = NULL, symmetric = FALSE, fixdiag = NULL
#! , vw = NULL, pack = TRUE, \dots)
#! }
#! \arguments{
#!   \item{x}{ an appropriate object of the class for which we dispatched }
#!   \item{envir}{ the environment in which to evaluate components of the index expression }
#!   \item{maxindex}{ maximum positive indexposition \code{maxindex}, is needed with negative indices, if vw or dim is given, maxindex is calculated automatically }
#!   \item{names}{ the \code{\link[ff:names.ff]{names}} of the indexed vector for character indexing }
#!   \item{dim}{ the \code{\link[ff:dim.ff]{dim}} of the indexed matrix to be stored within the \code{\link{hi}} object }
#!   \item{dimorder}{ the \code{\link{dimorder}} of the indexed matrix to be stored within the \code{\link{hi}} object, may convert interpretation of \code{x} }
#!   \item{symmetric}{ the \code{\link{symmetric}} of the indexed matrix to be stored within the \code{\link{hi}} object }
#!   \item{fixdiag}{ the \code{\link{fixdiag}} of the indexed matrix to be stored within the \code{\link{hi}} object }
#!   \item{vw}{ the virtual window \code{\link{vw}} of the indexed vector or matrix to be stored within the \code{\link{hi}} object, see details }
#!   \item{vw.convert}{ FALSE to prevent doubly virtual window conversion, this is needed for some internal calls that have done the virtual window conversion already, see details }
#!   \item{dimorder.convert}{ FALSE to prevent doubly dimorder conversion, this is needed for some internal calls that have done the dimorder conversion already, see details }
#!   \item{NAs}{ a vector of NA positions to be stored \code{\link[bit]{rlepack}ed}, not fully supported yet }
#!   \item{pack}{ FALSE to prevent \code{\link[bit]{rlepack}ing}, note that this is a hint rather than a guarantee, \code{as.hi.bit} might ignore this }
#!   \item{range}{ NULL or a vector with two elements indicating first and last position to be converted from 'bit' to 'hi' }
#!   \item{\dots}{ further argument passed from generic to method or from wrapper method to \code{as.hi.integer} }
#! }
#! \details{
#!   The generic dispatches appropriately, \code{as.hi.hi} returns an \code{\link{hi}} object unchanged,
#!   \code{as.hi.call} tries to \code{\link{hiparse}} instead of evaluate its input in order to save RAM.
#!   If parsing is successfull \code{as.hi.call} will ignore its argument \code{pack} and always pack unless the subscript is too small to do so.
#!   If parsing fails it evaluates the index expression and dispatches again to one of the other methods.
#!   \code{as.hi.name} and \code{as.hi.(} are wrappers to \code{as.hi.call}.
#!   \code{as.hi.integer} is the workhorse for coercing evaluated expressions
#!   and \code{as.hi.which} is a wrapper removing the \code{which} class attribute.
#!   \code{as.hi.double}, \code{as.hi.logical} and \code{as.hi.character} are also wrappers to \code{as.hi.integer},
#!   but note that \code{as.hi.logical} is not memory efficient because it expands \emph{all} positions and then applies logical subscripting.
#!   \cr
#!   \code{as.hi.matrix} calls \code{\link{arrayIndex2vectorIndex}} and then \code{as.hi.integer} to interpret and preprocess matrix indices.
#!   \cr
#!   If the \code{dim} and \code{dimorder} parameter indicate a non-standard dimorder (\code{\link{dimorderStandard}}), the index information in \code{x} is converted from a standard dimorder interpretation to the requested \code{\link{dimorder}}.
#!   \cr
#!   If the \code{vw} parameter is used, the index information in \code{x} is interpreted relative to the virtual window but stored relative to the abolute origin.
#!   Back-coercion via \code{\link{as.integer.hi}} and friends will again return the index information relative to the virtual window, thus retaining symmetry and transparency of the viurtual window to the user.
#!   \cr
#!   You can use \code{\link[ff:length.hi]{length}} to query the index length (possibly length of negative subscripts),
#!   \code{\link[ff:poslength.hi]{poslength}} to query the number of selected elements (even with negative subscripts),
#!   and \code{\link[ff:maxindex.hi]{maxindex}} to query the largest possible index position (within virtual window, if present)
#!   \cr
#!   Duplicated negative indices are removed and will not be recovered by \code{\link{as.integer.hi}}.
#! }
#! \value{
#!   an object of class \code{\link{hi}}
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ Avoid changing the Hybrid Index representation, this might crash the \code{\link{[.ff}} subscripting. }
#! \seealso{ \code{\link{hi}} for the Hybrid Index class, \code{\link{hiparse}} for parsing details, \code{\link{as.integer.hi}} for back-coercion, \code{\link{[.ff}} for ff subscripting }
#! \examples{
#!   message("integer indexing with and without rel-packing")
#!   as.hi(1:12)
#!   as.hi(1:12, pack=FALSE)
#!   message("if index is double, the wrapper method just converts to integer")
#!   as.hi(as.double(1:12))
#!   message("if index is character, the wrapper method just converts to integer")
#!   as.hi(c("a","b","c"), names=letters)
#!   message("negative index must use maxindex (or vw)")
#!   as.hi(-(1:3), maxindex=12)
#!   message("logical index can use maxindex")
#!   as.hi(c(FALSE, FALSE, TRUE, TRUE))
#!   as.hi(c(FALSE, FALSE, TRUE, TRUE), maxindex=12)
#!
#!   message("matrix index")
#!   x <- matrix(1:12, 6)
#!   as.hi(rbind(c(1,1), c(1,2), c(2,1)), dim=dim(x))
#!
#!   message("first ten positions within virtual window")
#!   i <- as.hi(1:10, vw=c(10, 80, 10))
#!   i
#!   message("back-coerce relativ to virtual window")
#!   as.integer(i)
#!   message("back-coerce relativ to absolute origin")
#!   as.integer(i, vw.convert=FALSE)
#!
#!   message("parsed index expressions save index RAM")
#!     as.hi(quote(1:1000000000))
#! \dontrun{
#!   message("compare to RAM requirement when the index experssion is evaluated")
#!     as.hi(1:1000000000)
#! }
#!
#! message("example of parsable index expression")
#!   a <- seq(100, 200, 20)
#!   as.hi(substitute(c(1:5, 4:9, a)))
#!   hi(c(1,4, 100),c(5,9, 200), by=c(1,1,20))
#!
#! message("two examples of index expression temporarily expanded to full length due to 
#! non-supported use of brackets '(' and mathematical operators '+' accepting token")
#! message("example1: accepted token but aborted parsing because length>16")
#!   as.hi(quote(1+(1:16)))
#! message("example1: rejected token and aborted parsing because length>16")
#!   as.hi(quote(1+(1:17)))
#! }
#! \keyword{ IO }
#! \keyword{ data }



as.hi.hi <- function(x, ...){
  x
}


as.hi.name <- function(x, envir=parent.frame(), ...){
  as.hi(eval(x, envir=envir), ...)
}

"as.hi.(" <- function(x, envir=parent.frame(), ...){
  as.hi.call(x[[2]], envir=envir, ...)
}

as.hi.call <- function(
  x
, maxindex    = NA
, dim         = NULL
, dimorder    = NULL
, vw          = NULL
, vw.convert  = TRUE
, pack        = TRUE
, envir       = parent.frame()
, ...
){

  if ((!is.null(dim) && !dimorderStandard(dimorder)) || !is.null(dim(vw)))
    return(as.hi(eval(x, envir=envir), maxindex=maxindex, dim=dim, dimorder=dimorder, vw=vw, vw.convert=vw.convert, pack=pack, ...))


  #message("DEBUGINFO: trying hiparse")
  r <- try(hiparse(x, envir=envir), silent=TRUE)
  if (inherits(r,"try-error")){
    #message("DEBUGINFO: hiparse failed, evaluating the index expression and dispatching again")
    return(as.hi(eval(x, envir=envir), maxindex=maxindex, dim=dim, dimorder=dimorder, vw=vw, vw.convert=vw.convert, pack=pack, ...))
  }

  if (is.null(vw))
    vw.convert <- FALSE
  else{
    storage.mode(vw) <- "integer"
  }

  minindex <- 1L
  if (is.na(maxindex)){
    if(is.null(dim))
      maxindex <- maxindex(x)
    else
      maxindex <- as.integer(prod(dim))
  }else{
    maxindex <- as.integer(maxindex)
  }

  if (is.na(r$first)){
    x <- list(first=as.integer(NA), dat=integer(), last=as.integer(NA))
    ix <- NULL
    re <- FALSE
    n <- 0L
  }else{
    nl <- length(r$lengths)
    n <- sum(r$lengths) + 1L
    # test for sorted
    tab <- tabulate(sign(r$values)+2, 3)
    if (tab[1] && tab[3]){ # unsorted in both directions
      re <- FALSE
      x <- as.integer(cumsum(c(r$first, rep(r$values, r$lengths))))
      x <- sort.int(x, index.return=TRUE, method="quick")
      ix <- x$ix
      x <- rlepack(x$x, pack=pack)
      #ix <- 1:n
      #radixorder(x, ix)
      #x <- rlepack(x[ix], pack = pack)

    }else{
      if (nl){
        pack <- 2*length(r$lengths)<n
      }else
        pack <- FALSE
      if (pack){
        dat <- list(lengths=r$lengths, values=r$values)
        class(dat) <- "rle"
      }else{
        dat <- as.integer(cumsum(c(r$first, rep(r$values, r$lengths))))
      }
      x <- list(first=r$first, dat=dat, last=r$last)
      ix <- NULL
      if (tab[1]){  # sorted descending
        re <- TRUE
        x <- rev.rlepack(x)
      }else{        # sorted ascending
        re <- FALSE
      }
    }
    # xx this ifelse section copied 1L from hi
    if (x$last < 0) {
      if (is.na(maxindex))
          stop("maxindex is required with negative subscripts")
      if ( -x$first > maxindex )
          stop("negative subscripts out of range")

      re <- FALSE
      ix <- NULL

      if (vw.convert){
        x$first <- x$first - vw[1]
        x$last <- x$last - vw[1]
        if (inherits(x$dat, "rle")){
          n <- sum(x$dat$lengths) + 1L
        }else{
          x$dat <- x$dat - vw[1]
          n <- length(x$dat)
        }
      }else{
        if (inherits(x$dat, "rle")){
          n <- sum(x$dat$lengths) + 1L
        }else{
          n <- length(x$dat)
        }
      }
    }else if (x$first > 0){
      if (!is.na(maxindex) && x$last > maxindex )
          stop("positive subscripts out of range")
      if (vw.convert){
        x$first <- vw[1] + x$first
        x$last <- vw[1] + x$last
        if (inherits(x$dat, "rle")){
          n <- sum(x$dat$lengths) + 1L
        }else{
          x$dat <- vw[1] + x$dat
          n <- length(x$dat)
        }
      }else{
        if (inherits(x$dat, "rle")){
          n <- sum(x$dat$lengths) + 1L
        }else{
          n <- length(x$dat)
        }
      }
    }else{
        stop("0s and mixed positive/negative subscripts not allowed")
    }
  }

  if (!is.null(vw)){
    if (is.null(dim)){
      # minindex..maxindex represents the window of allowed values (used by the C-code in case of negative subscripts for enumerating all positive subscripts)
      minindex <- vw[1] + 1L
      maxindex <- vw[1] + vw[2]
    }else{
      # NOTE that negative subscripts cannot be handled in a (vw && dim)-context (enumerating all positive subscripts is not simply minindex..maxindex)
      minindex <- 1L
      maxindex <- as.integer(prod(colSums(vw)))
    }
  }
  ret <- list(
    x         = x
  , ix        = ix
  , re        = re
  , minindex  = minindex
  , maxindex  = maxindex
  , length    = n
  , dim       = NULL
  , dimorder  = NULL
  , symmetric = FALSE
  , fixdiag   = NULL
  , vw        = vw
  , NAs       =  NULL
  )
  class(ret) <- "hi"
  return(ret)
}


as.hi.integer <- function(
  x
, maxindex    = NA
, dim         = NULL
, dimorder    = NULL
, symmetric   = FALSE
, fixdiag     = NULL
, vw          = NULL
, vw.convert  = TRUE        # as.hi.matrix sets this to false in order to avoid applying vw twice
, dimorder.convert  = TRUE  # as.hi.matrix sets this to false in order to avoid dimorder conversion twice
, pack        = TRUE
, NAs         = NULL
, ... # dummy to keep R CMD check quiet
){

  n <- length(x)

  if (is.null(vw))
    vw.convert <- FALSE
  else{
    storage.mode(vw) <- "integer"
    if (is.null(dim) && !is.null(dim(vw)))
      dim <- vw[2,]
  }

  # these are still the limits to be checked from user perspective
  minindex <- 1L
  if (is.na(maxindex)){
    if(is.null(dim))
      maxindex <- maxindex(x)
    else
      maxindex <- as.integer(prod(dim))
  }else{
    maxindex <- as.integer(maxindex)
  }

  if (n){
    if (is.null(dim) || dimorderStandard(dimorder))
      dimorder.convert <- FALSE
    prechecked <- dimorder.convert || (vw.convert && !( is.null(dim) || dimorderStandard(dimorder) ))
    if (prechecked){
      # need dimorder conversion, i.e., sorting pre-conversion != sorting post-conversion
      # since we have no sorting pre-conversion, we cannot simply check the most extreme value, we have to check all values
      # therefor we check all values pre-conversion
      if (all(x<0, na.rm=TRUE)){
        if (any(x < -maxindex, na.rm=TRUE))
          stop("negative subscripts out of range")
        x <- (1:maxindex)[x]  # convert to positive indexes because we cannot enumerate
      }else if (all(x>0, na.rm=TRUE)){
        if (any(x > maxindex, na.rm=TRUE))
          stop("positive subscripts out of range")
      }else
        stop("0s and mixed positive/negative subscripts not allowed")
      x <- arrayIndex2vectorIndex(vectorIndex2arrayIndex(x, dim=dim), dim=dim, dimorder=dimorder, vw=vw)
      vw.convert <- FALSE
      # these are already the final limits from file perspective (make sure the converted values pass the (redundant) test further down)
      if (is.null(vw))
        maxindex <- prod(dim)
      else
        maxindex <- prod(colSums(vw))
    }

    isasc <- intisasc(x)
    if (is.na(isasc))
      stop("NAs in as.hi.integer")
    if (isasc){
      ix <- NULL
      re <- FALSE
    }else{
      if (intisdesc(x)){
        x <- rev(x)
        ix <- NULL
        re <- TRUE
      }else{
        x <- sort.int(x, index.return=TRUE, method="quick")
        ix <- x$ix
        x <- x$x
        #ix <- 1:n
        #radixorder(x, ix)
        #x <- x[ix]
        re <- FALSE
      }
    }

    # after sorting the range-checks can be done on the extremes only
    if (x[n]<0){  # not possible after prechecked
      if (is.na(maxindex)){
        if (vw.convert && is.null(dim))
          maxindex <- vw[[2]]
        else
          stop("maxindex is required with negative subscripts")
      }
      if ( -x[1] > maxindex )
        stop("negative subscripts out of range")

      ix <- NULL
      re <- FALSE

      x <- unique(x)
      n <- length(x)
      if (vw.convert){
        # convert window positions to absolute positions
        if (is.null(dim)){
          x <- x - vw[1]
        }else{
          # (vw && dim)-context: convert negative indices to positive ones
          x <- (1:maxindex)[x]
          n <- length(x)
          if (n)
            x <- arrayIndex2vectorIndex(vectorIndex2arrayIndex(x, dim=dim, dimorder=dimorder), dimorder=dimorder, vw=vw)
        }
      }
    }else if (x[1]>0){
      if ( !is.na(maxindex) && x[n] > maxindex )
        stop("positive subscripts out of range")
      if (vw.convert){
        # convert window positions to absolute positions
        if (is.null(dim)){
          x <- vw[1] + x
        }else{
          x <- arrayIndex2vectorIndex(vectorIndex2arrayIndex(x, dim=dim, dimorder=dimorder), dimorder=dimorder, vw=vw)
        }
      }
    }else{
      stop("0s and mixed positive/negative subscripts not allowed")
    }

    x <- rlepack(x, pack=pack)
    #We could restrict compression to the case where we spent time and RAM on sorting anyhow
    #x <- rlepack(x, pack=if (is.null(ix)) FALSE else pack)

  }else{ # no data
    x  <- list(first=as.integer(NA), dat=integer(), last=as.integer(NA))
    ix <- NULL
    re <- FALSE
  }

  # fix the final limits
  if (!is.null(vw)){
    if (is.null(dim)){
      # minindex..maxindex represents the window of allowed values (used by the C-code in case of negative subscripts for enumerating all positive subscripts)
      minindex <- vw[1] + 1L
      maxindex <- vw[1] + vw[2]
    }else{
      # NOTE that negative subscripts cannot be handled in a (vw && dim)-context or in a non-standard dimorder-context (enumerating all positive subscripts is not simply minindex..maxindex)
      maxindex <- as.integer(prod(colSums(vw)))
    }
  }

  r <- list(
    x         = x
  , ix        = ix
  , re        = re
  , minindex  = minindex
  , maxindex  = maxindex
  , length    = n
  , dim       = dim
  , dimorder  = dimorder
  , symmetric = symmetric
  , fixdiag   = fixdiag
  , vw        = vw
  , NAs       = NAs
  )
  class(r) <- "hi"
  r
}


as.hi.which <-  function(x, ...)
  as.hi.integer(unclass(x), ...)


if (FALSE){
  dim <- 3:4
  dimorder <- 1:2
  vw <- rbind(c(1,1), dim, c(1,1))
  i <- 1:prod(dim)
  m <- vectorIndex2arrayIndex(i, dim=dim)
  p <- arrayIndex2vectorIndex(m, dim=dim, dimorder=dimorder, vw=vw)

  m
  vectorIndex2arrayIndex(p, dim=dim, dimorder=dimorder, vw=vw)

  h <- as.hi(m, dim=dim, dimorder=dimorder, vw=vw)
  str(h)
  p

  as.integer(h)
  i

}



as.hi.matrix <- function(x, dim, dimorder=NULL, symmetric=FALSE, fixdiag=NULL, vw=NULL, pack=TRUE
, ... # dummy to keep R CMD check quiet
){
  if (is.null(vw)){
    maxindex <- as.integer(prod(dim))
  }else{
    maxindex <- as.integer(prod(colSums(vw)))
  }
  if (nrow(x)){
    if (x[1]<0)
      stop("matrix subscripts must be positive")
    if (symmetric){
      i <- symmIndex2vectorIndex(x, dim=dim, fixdiag=fixdiag)
      if (is.null(fixdiag)){
        ret <- as.hi.integer(i, maxindex=maxindex, dim=dim, symmetric=symmetric, fixdiag=fixdiag, vw=vw, pack=pack)
      }else{
        isna <- is.na(i)
        NAs <- (1:length(i))[isna]
        if (length(NAs))
          ret <- as.hi.integer(i[!isna], maxindex=maxindex, dim=dim, symmetric=symmetric, fixdiag=fixdiag, vw=vw, pack=pack, NAs=rlepack(NAs))
        else
          ret <- as.hi.integer(i, maxindex=maxindex, dim=dim, symmetric=symmetric, fixdiag=fixdiag, vw=vw, pack=pack)
      }
    }else{
      ret <- as.hi.integer(
        arrayIndex2vectorIndex(x, dim=dim, dimorder=dimorder, vw=vw)
      , maxindex=maxindex
      , dim=dim
      , dimorder=dimorder
      , symmetric=symmetric
      , fixdiag=fixdiag
      , vw=vw
      , vw.convert=FALSE
      , dimorder.convert=FALSE
      , pack=pack
      )
    }
  }else{
    ret <- as.hi.integer(integer(), maxindex=maxindex, dim=dim, dimorder=dimorder, symmetric=symmetric, fixdiag=fixdiag, vw=vw, pack=pack)
  }
  ret
}

as.hi.logical <- function(
  x
, maxindex  = NA
, dim       = NULL
, vw        = NULL
, pack      = TRUE
, ... # dummy to keep R CMD check quiet
){
  if(is.null(dim)){
    if (is.na(maxindex))
      maxindex <- length(x)
    else
      maxindex <- as.integer(maxindex)
  }else{
    maxindex <- as.integer(prod(dim))
  }
  if (length(x)>maxindex)
    stop("as.hi.logical longer than maxindex")
  if (maxindex>0){
    x <- (1:maxindex)[rep(x, length=maxindex)]
  }else{
    x <- integer()
  }
  return(as.hi.integer(
    x
  , maxindex  = maxindex  # if !is.null(vw) maxindex is wrong, thus we rely on as.hi.integer ignoring maxindex in this case !!
  , dim       = dim
  , vw        = vw
  , pack      = pack
  ))
}

as.hi.double <- function(x, ...){
  #silent as usually, thus not: warning("converting doubles to integer in as.hi")
  as.hi.integer(as.integer(x), ...)
}

# used for character subsetting
# needs names:
as.hi.character <- function(x
, names           # either character vector or some object of class that has "[" defined and returns integer positions (such as a named integer vector or class index)
, vw = NULL
, vw.convert=TRUE   # if names refers to the vw-window, the default vw.convert=TRUE is fine, if names refers to the total object, set vw.convert = FALSE
, ...
){
  #if (inherits(names, "fffc"))
  #  as.hi.integer(match.fffc(x, names), vw=vw, vw.convert=vw.convert, ...)
  #else
  if (is.atomic(names) && is.character(names))
    as.hi.integer(match(x, names), vw=vw, vw.convert=vw.convert, ...)
  else
    as.hi.integer(names[x], vw=vw, vw.convert=vw.convert, ...)
}


# -- reverting hi to original (as far as possible) -----------------------------------------------------------

#! \name{as.integer.hi}
#! \alias{as.which.hi}
#! \alias{as.bitwhich.hi}
#! \alias{as.bit.hi}
#! \alias{as.integer.hi}
#! \alias{as.logical.hi}
#! \alias{as.character.hi}
#! \alias{as.matrix.hi}
#! \title{ Hybrid Index, coercing from }
#! \description{
#!   Functions that (back-)convert an \code{\link{hi}} object to the respective subscripting information.
#! }
#! \usage{
#! \method{as.which}{hi}(x, \dots)
#! \method{as.bitwhich}{hi}(x, \dots)
#! \method{as.bit}{hi}(x, \dots)
#! \method{as.integer}{hi}(x, vw.convert = TRUE, \dots)
#! \method{as.logical}{hi}(x, maxindex = NULL, \dots)
#! \method{as.character}{hi}(x, names, vw.convert = TRUE, \dots)
#! \method{as.matrix}{hi}(x, dim = x$dim, dimorder = x$dimorder
#! , vw = x$vw, symmetric = x$symmetric, fixdiag = x$fixdiag, \dots)
#! }
#! \arguments{
#!   \item{x}{ an object of class \code{\link{hi}} }
#!   \item{maxindex}{ the \code{\link{length}} of the subscripted object (needed for logical output) }
#!   \item{names}{ the \code{\link{names}} vector of the subscripted object }
#!   \item{dim}{ the \code{\link{dim}} of the subscripted object }
#!   \item{dimorder}{ the \code{\link{dimorder}} of the subscripted object }
#!   \item{vw}{ the virtual window \code{\link{vw}} of the subscripted object }
#!   \item{vw.convert}{ \code{vw.convert} }
#!   \item{symmetric}{ TRUE if the subscripted matrix is \code{\link{symmetric}} }
#!   \item{fixdiag}{ TRUE if the subscripted matrix has \code{\link{fixdiag}} }
#!   \item{\dots}{ further arguments passed }
#! }
#! \value{
#!   \command{as.integer.hi} returns an integer vector, see \code{\link{as.hi.integer}}.
#!   \command{as.logical.hi} returns an logical vector, see \code{\link{as.hi.logical}}.
#!   \command{as.character.hi} returns a character vector, see \code{\link{as.hi.character}}.
#!   \command{as.matrix.hi} returns a matrix index, see \code{\link{as.hi.matrix}}.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{hi}}, \code{\link{as.hi}} }
#! \examples{
#!   x <- 1:6
#!   names(x) <- letters[1:6]
#!   as.integer(as.hi(c(1:3)))
#!   as.logical(as.hi(c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE)))
#!   as.character(as.hi(letters[1:3], names=names(x)), names=names(x))
#!   x <- matrix(1:12, 6)
#!   as.matrix(as.hi(rbind(c(1,1), c(1,2), c(2,1)), dim=dim(x)), dim=dim(x))
#! }
#! \keyword{ IO }
#! \keyword{ data }

# note that the result may be changed for negative subscripts (sorted ascending, duplicates removed)
# original order of subscripts NEED NOT be restored and MUST NOT be restored in order to not create long hi$ix component
as.integer.hi <- function(
  x
, vw.convert=TRUE   # set to FALSE when called from as.matrix.hi in order to avoid double conversion
, ... # dummy to keep R CMD check quiet
){
  if (x$length){
    ret <- unsort.hi(rleunpack(x$x), x)
    if (is.null(x$dim)){
      if (!is.null(x$vw) && vw.convert){
        if (ret[1]<0){
          ret <- ret + x$vw[1]
        }else{
          ret <- ret - x$vw[1]
        }
      }
    }else{
      if (!is.null(x$vw) && vw.convert){
        # we know that subscripts must be positive in this case
          ret <- arrayIndex2vectorIndex(vectorIndex2arrayIndex(ret, dimorder=x$dimorder, vw=x$vw), dim=x$vw[2,])
      }else{
        if (!dimorderStandard(x$dimorder))
          ret <- arrayIndex2vectorIndex(vectorIndex2arrayIndex(ret, dim=x$dim, dimorder=x$dimorder), dim=x$dim)
      }
    }
  }else{
    ret <- integer()
  }
  ret
}

as.which.hi <- function(x, ...){
  i <- as.integer(x, ...)
  if (length(i) && i[[1]]<0)
    i <- (1:maxindex(x))[i]
  class(i) <- "which"
  i
}


as.matrix.hi <- function(
  x
, dim       = x$dim
, dimorder  = x$dimorder
, vw       = x$vw
, symmetric = x$symmetric
, fixdiag   = x$fixdiag
, ... # dummy to keep R CMD check quiet
){
  if (x$length){
    if (is.null(dim))
      stop("need dim to return matrix subscripts")
    if (x$x$first<0)
      stop("matrix subscripts must be positive")
    if (symmetric){
      if (is.null(fixdiag)){
        stop("not yet implemented for symmetric matices with fixdiag")
      }else{
        stop("not yet implemented for symmetric matices without fixdiag (redundant diagonal)")
      }
    }else{
      ret <- unsort.hi(rleunpack(x$x), x)
      ret <- vectorIndex2arrayIndex(ret, dim=dim, dimorder=dimorder, vw=vw)
    }
    ret
  }else{
    matrix(integer(), 0, length(x$dim))
  }
}


# note that result may be longer due to recycling
as.logical.hi <- function(
  x
, maxindex=NULL
, ... # dummy to keep R CMD check quiet
){
  if (is.null(maxindex))
    maxindex <- maxindex(x)
  if (is.na(maxindex))
    stop("can't make logical without knowing vector length")
  ret <- rep(FALSE, maxindex)
  ret[(1:maxindex)[as.integer.hi(x)]] <- TRUE
  ret
}

as.character.hi <- function(
  x
, names
, vw.convert=TRUE  # set to FALSE if the names do not correspond to window but to total object
, ... # dummy to keep R CMD check quiet
){
  names[as.integer.hi(x, vw.convert=vw.convert)]
}


# -- querying the 'length' of a hi object ----------------------------------------------


#! \name{length.hi}
#! \alias{length.hi}
#! \alias{poslength}
#! \alias{poslength.hi}
#! \alias{poslength.ri}
#! \alias{poslength.bit}
#! \alias{poslength.bitwhich}
#! \alias{poslength.logical}
#! \alias{poslength.default}
#! \alias{maxindex}
#! \alias{maxindex.hi}
#! \alias{maxindex.ri}
#! \alias{maxindex.bit}
#! \alias{maxindex.bitwhich}
#! \alias{maxindex.logical}
#! \alias{maxindex.default}
#! \title{ Hybrid Index, querying }
#! \description{
#!   Functions to query some index attributes
#! }
#! \usage{
#! \method{length}{hi}(x)
#! maxindex(x, \dots)
#! poslength(x, \dots)
#! \method{maxindex}{hi}(x, \dots)
#! \method{maxindex}{ri}(x, \dots)
#! \method{maxindex}{bit}(x, \dots)
#! \method{maxindex}{bitwhich}(x, \dots)
#! \method{maxindex}{logical}(x, \dots)
#! \method{maxindex}{default}(x, \dots)
#! \method{poslength}{hi}(x, \dots)
#! \method{poslength}{ri}(x, \dots)
#! \method{poslength}{bit}(x, \dots)
#! \method{poslength}{bitwhich}(x, \dots)
#! \method{poslength}{logical}(x, \dots)
#! \method{poslength}{default}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an object of class \code{\link{hi}} }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!   \command{length.hi} returns the number of the subsript elements in the index (even if they are negative).
#!   By contrast the generic \command{poslength} returns the number of selected elements (which for negative indices is \code{maxindex(x) - length(unique(x))}).
#!   The generic \command{maxindex} returns the highest possible index position.
#! }
#! \value{
#!   an integer scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ duplicated negative indices are removed }
#! \seealso{ \code{\link{hi}}, \code{\link{as.hi}}, \code{\link{length.ff}}, \code{\link{length}} }
#! \examples{
#!   length(as.hi(-1, maxindex=12))
#!   poslength(as.hi(-1, maxindex=12))
#!   maxindex(as.hi(-1, maxindex=12))
#!   message("note that")
#!   length(as.hi(c(-1, -1), maxindex=12))
#!   length(as.hi(c(1,1), maxindex=12))
#! }
#! \keyword{ IO }
#! \keyword{ data }


# be aware that this is is the length of the pos/neg integer index, not necessarily the length of an original logical (no FALSEs, recycled) or negative (no duplicates)
length.hi <- function(x){
  x$length
}

# this give the length of the object to be subscripted (if known)
# must always be known after as.hi.logical
maxindex.hi <- function(
  x
, ... # dummy to keep R CMD check quiet
)
{
  if (is.null(x$vw))
    x$maxindex
  else{
    if (is.null(x$dim))
      x$vw[2]
    else
      as.integer(prod(x$vw[2,]))
  }
}
poslength.hi <- function(
  x
, ... # dummy to keep R CMD check quiet
){
  if (is.na(x$x$first))
    0L
  else if (x$x$first<0){
    if (is.na(x$maxindex))
      stop("poslength.hi requires maxindex")
    maxindex.hi(x) - x$length
  }else
    x$length
}

if (!exists("maxindex.default"))
  maxindex.default <- function(x
  , ... # dummy to keep R CMD check quiet
  ){
    mi <- attr(x, "maxindex")
    if (is.null(mi))
      as.integer(NA)
    else
      mi
  }
if (!exists("poslength.default"))
  poslength.default <- function(x
  , ... # dummy to keep R CMD check quiet
  ){
    mi <- maxindex(x)
    if (is.null(mi))
      as.integer(NA)
    else
      mi - length(x)
  }

if (!exists("maxindex.logical"))
  maxindex.logical <- function(x
  , ... # dummy to keep R CMD check quiet
  )
  length(x)
if (!exists("poslength.logical"))
  poslength.logical <- function(x
  , ... # dummy to keep R CMD check quiet
  ){
    sum(x, na.rm=TRUE)
  }




#! \name{unsort}
#! \alias{unsort}
#! \alias{unsort.hi}
#! \alias{unsort.ahi}
#! \alias{subscript2integer}
#! \title{ Hybrid Index, internal utilities }
#! \description{
#!   Non-documented internal utilities that might change
#! }
#! \usage{
#! unsort(x, ix)
#! unsort.hi(x, index)
#! unsort.ahi(x, index, ixre = any(sapply(index, function(i) {
#!     if (is.null(i$ix)) {
#!         if (i$re) TRUE else FALSE
#!     } else {
#!         TRUE
#!     }
#! })), ix = lapply(index, function(i) {
#!     if (is.null(i$ix)) {
#!         if (i$re)
#!             orig <- rev(1:poslength(i))
#!         else orig <- 1:poslength(i)
#!     }
#!     else {
#!         orig <- i$ix
#!     }
#!     orig
#! }))
#! subscript2integer(x, maxindex = NULL, names = NULL)
#! }
#! \arguments{
#!   \item{x}{ \code{x} }
#!   \item{ix}{ \code{ix} }
#!   \item{ixre}{ \code{ixre} }
#!   \item{index}{ \code{index} }
#!   \item{maxindex}{ \code{maxindex} }
#!   \item{names}{ \code{names} }
#! }
#! \details{
#!   These are utility functions for restoring original order after sorting.
#!   For now we 'mimic' the intuitive but wrong argument order of match()
#!   which should rather have the 'table' argument as its first argument,
#!   then one could properly method-dispatch on the type of table.
#!   xx We might change to proper 'unsort' generic, but then we have to change argument order.
#! }
#! \value{
#!  undefined
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{hi}}, \code{\link{as.hi}} }
#! \keyword{ IO }
#! \keyword{ data }



# not actually used
unsort <- function(
  x   # sorted values
, ix  # needed to restore orig order
 ){
  orig <- vector(mode=storage.mode(x), length=length(x))
  orig[ix] <- x
  orig
}

unsort.hi <- function(
  x       # a vector of values in hi sorting
, index   # a hi object
){
  if (is.null(index$ix)){
    if (index$re)
      orig <- rev(x)
    else
      orig <- x
  }else{
    orig <- vector(mode=storage.mode(x), length=length(x))
    orig[index$ix] <- x
  }
  orig
}

unsort.ahi <- function(
  x       # an m-array of values in multi-hi sorting
, index   # a list of m hi indices
, ixre = any(sapply(index, function(i){
    if (is.null(i$ix)){
      if (i$re)
        TRUE
      else
        FALSE
    }else{
      TRUE
    }
  }))
, ix = lapply(index, function(i){
      if (is.null(i$ix)){
        if (i$re)
          orig <- rev(1:poslength(i))
        else
          orig <- 1:poslength(i)
      }else{
        orig <- i$ix
      }
      orig
    })
){
  if (ixre){
    x <- do.call("[<-", c(list(x=x), ix, list(value=x)))
  }
  x
}

# used in subset.ff_dist
# does not handle matrix subscripts
subscript2integer <- function(
  x
, maxindex=NULL
, names=NULL
){
  if(any(is.na(x)))
    stop("NAs not allowed in ff subscripting")
  if (is.character(x)){
    if (is.null(names))
      stop("need names")
    match(x, names)
  }else if(is.logical(x)){
    if (is.null(maxindex))
      stop("need maxindex with logical subscripts")
    (1:maxindex)[x]
  }else{
    if (is.double(x))
      x <- as.integer(x)
    tab <- tabulate(sign(x)+2, 3)
    if (tab[[2]])
      stop("no zeros allowed in ff subscripts")
    if (tab[[1]] && tab[[2]])
      stop("mixing negative and positive subscripts is not alllowed")
    if (tab[[1]]){
    if (is.null(maxindex))
      stop("need maxindex with negative subscripts")
      (1:maxindex)[x]
    }else{
      x
    }
  }
}




# Example
if (FALSE){
  a <- seq(100, 200, 20)
  as.hi(substitute(c(1:5, 4:9, a)))
  hi(c(1,4, 100),c(5,9, 200), by=c(1,1,20))

  as.hi(c(1:5, 4:9, a))
  x <- c(1:5, 4:9, a)
  as.hi(x)
  as.hi(substitute(x))

  as.integer(as.hi(x))
  as.logical(as.hi(x))
  as.logical(as.hi(x, maxindex=200))
  length(as.hi.integer(x))
  maxindex(as.hi(x))
  poslength(as.hi(x, maxindex=200))

  library(regtest)
  # parsing has some overhead ...
  timefactor(as.hi(substitute(c(1:4, 5:9, a))),  hi(c(1,5,100),c(4,9, 200), by=c(1,1,20)), 1000, 1000)
  # ... but with long sequences
  timefactor(as.hi(substitute(c(1:4, a, 500:999999))),  as.hi(c(1:4, a, 500:999999)), 100, 1)

  s1 <- hi(c(1,4, 200),c(5,9, 100), by=c(1,1,-20))
  s2 <- as.hi(substitute(c(1:5, 4:9, a)))
  s3 <- as.hi(c(1:5, 4:9, a))
  identical(s1, s2)
  identical(s3, s2)
  identical(as.integer(c(1:5, 4:9, a)), as.integer(s1))
  identical(as.integer(c(1:5, 4:9, a)), as.integer(s2))
  identical(as.integer(c(1:5, 4:9, a)), as.integer(s3))


  library(ff)
  n <- 10000000
  a <- ff(0L, length=n)
  #i <- c(1,n:2)
  #i <- sample(n)
  #save(i, file="c:/tmp/i.RData")
  load(file="c:/tmp/i.RData")
  memory.size(max=T)

  j <- rlepack(i) # x10
  memory.size(max=T)

  debug(as.hi.integer)
  j <- as.hi.integer(i)
  memory.size(max=T)

  system.time(j <- as.hi(quote(i)))

  x <- 20:29
  as.hi(quote((c(1, 3:10, x))))


  # C-coded with 33% trick
  load(file="c:/tmp/i.RData")
  memory.size(max=T)
  gc()
  system.time(j <- intrle(i))
  memory.size(max=T)


  # minus structure() = 7x RAM incl. input/output
  load(file="c:/tmp/i.RData")
  rle <-
  function (x)
  {
      if (!is.vector(x) && !is.list(x))
          stop("'x' must be an atomic vector")
      n <- length(x)
      if (n == 0)
          return(list(lengths = integer(0), values = x))
      y <- x[-1] != x[-n]
      i <- c(which(y | is.na(y)), n)
      ret <- list(lengths = diff(c(0L, i)), values = x[i])
      class(ret) = "rle"
      ret
  }
  gc()
  j <- rle(i)
  memory.size(max=T)

  # original = 9x RAM incl. input/output
  load(file="c:/tmp/i.RData")
  gc()
  j <- rle(i)
  memory.size(max=T)



}



