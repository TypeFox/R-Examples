# R layer of ff
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-09-03
# Last changed: 2007-10-25

# source("c:/mwp/eanalysis/ff/R/ff.R")
# package.skeleton("fftest", path="c:/tmp", list=ls("package:ff"))


if(getRversion() < "2.11.0")
    .POSIXct <- function(xx, tz = NULL)
    structure(xx, class = c("POSIXt", "POSIXct"), tzone = tz)


# --- ff info -----------------------------------------------------------

#! \name{ffxtensions}
#! \alias{ffxtensions}
#! \alias{ffsymmxtensions}
#! \title{ Test for availability of ff extensions }
#! \description{
#!   checks if this version of package ff supports ff extensions.
#! }
#! \usage{
#!  ffxtensions()
#!  ffsymmxtensions()
#! }
#! \value{
#!   logical scalar
#! }
#! \details{
#!   ff extensions are needed for certain bitcompressed vmodes and ff symm extensions for symmetric matrices.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{vmode}}%, \code{\link{symm}}
#! }
#! \examples{
#!   ffxtensions()
#!   ffsymmxtensions()
#! }
#! \keyword{ IO }
#! \keyword{ data }

ffxtensions <- function()
  .Call("ffxtensions", PACKAGE="ff")

ffsymmxtensions <- function()
  .Call("ffsymmxtensions", PACKAGE="ff")


#vector of those physical ff attributes that are kept in the ram version
ramphysical_includes <- c("pattern","filename","pagesize","caching","finalizer","finonexit")
#vector of those virtual ff attributes that are kept in the ram version
ramvirtual_includes <- c("Dimorder")


# vector of attributes that are not handled by the ramattribs mechanism
ramattribs_excludes <- c(
  "length"
, "dim"
, "dimnames"
, "names"
, "levels"
, "class"
, "vmode"
, "physical"
, "virtual"
)

# vector of classes that are not handled by the ramclass mechanism
ramclass_excludes <- c(
  "matrix"
, "array"
, "symm"
, "logical"
, "integer"
, "numeric"
)

caching_schemes <- c("mmnoflush","mmeachflush")

#! \name{is.ff}
#! \alias{is.ff}
#! \title{ Test for class ff }
#! \description{
#!   checks if x inherits from class "ff"
#! }
#! \usage{
#! is.ff(x)
#! }
#! \arguments{
#!   \item{x}{ any object }
#! }
#! \value{
#!   logical scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{inherits}}, \code{\link{as.ff}}, \code{\link{is.ffdf}} }
#! \examples{
#!   is.ff(integer())
#! }
#! \keyword{ IO }
#! \keyword{ data }

is.ff <- function(x)
{
  inherits(x,"ff")
}


#! \name{geterror.ff}
#! \alias{geterror.ff}
#! \alias{geterrstr.ff}
#! \title{ Get error and error string }
#! \description{
#!   Get last error code and error string that occured on an ff object.
#! }
#! \usage{
#! geterror.ff(x)
#! geterrstr.ff(x)
#! }
#! \arguments{
#!   \item{x}{ an ff object }
#! }
#! \value{
#!   \command{geterror.ff} returns an error integer code (no error = 0) and \command{geterrstr.ff} returns the error message (no error = "no error").
#! }
#! \author{ Jens Oehlschlägel, Daniel Adler (C++ back-end) }
#! \seealso{  \code{\link{ff}} }
#! \examples{
#!   x <- ff(1:12)
#!   geterror.ff(x)
#!   geterrstr.ff(x)
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

geterror.ff <- function(x)
{
  .Call("geterror", attr(x, "physical"), PACKAGE="ff")
}
geterrstr.ff <- function(x)
{
  .Call("geterrstr", attr(x, "physical"), PACKAGE="ff")
}


#! \name{filename}
#! \alias{filename}
#! \alias{filename.default}
#! \alias{filename.ff_pointer}
#! \alias{filename.ffdf}
#! \alias{filename<-}
#! \alias{filename<-.ff}
#! \alias{pattern}
#! \alias{pattern<-}
#! \alias{pattern.ff}
#! \alias{pattern<-.ff}
#! \alias{pattern<-.ffdf}
#! \title{ Get or set filename }
#! \description{
#!   Get or set filename from ram or \code{\link{ff}} object via the \code{filename} and \code{filename<-} generics
#!   or rename all files behind a \code{\link{ffdf}} using the \code{pattern<-} generic.
#! }
#! \usage{
#! filename(x, \dots)
#! filename(x, \dots) <- value
#! \method{filename}{default}(x, \dots)
#! \method{filename}{ff_pointer}(x, \dots)
#! \method{filename}{ffdf}(x, \dots)
#! \method{filename}{ff}(x, \dots) <- value
#! pattern(x, \dots)
#! pattern(x, \dots) <- value
#! \method{pattern}{ff}(x, \dots)
#! \method{pattern}{ff}(x, \dots) <- value
#! \method{pattern}{ffdf}(x, \dots) <- value
#! }
#! \arguments{
#!   \item{x}{ a ram or ff  object, or for pattern assignment only - a ffdf object }
#!   \item{value}{ a new filename }
#!   \item{\dots}{ dummy to keep R CMD CHECK quiet }
#! }
#! \value{
#!   \code{filename} and \code{pattern} return a character filename or pattern.
#!   For \code{\link{ffdf}} returns a list with one filename element for each \code{\link[=physical.ffdf]{physical}} component.
#!   The assignment functions return the changed object, which will keep the change even without re-assigning the return-value
#! }
#! \details{
#!   Assigning a \code{filename<-} means renaming the corresponding file on disk - even for ram objects. If that fails, the assignment fails.
#!   If a file is moved in or out of \code{getOption("fftempdir")} the \code{\link{finalizer}}  is changed accordingly to 'delete' in \code{getOption("fftempdir")} and 'close' otherwise.
#!   \cr
#!   A \code{pattern} is an incomplete filename (optional path and optional filename-prefix) that is turned to filenames by
#!   adding a random string using and optionally an extension from optionally an extension from \code{getOption("ffextension")} (see \code{\link{fftempfile}}).
#!   \code{filename<-} exhibits R's standard behaviour of considering "filename" and "./filename" both to be located in \code{\link{getwd}}.
#!   By constrast \code{pattern<-} will create "filename" without path in \code{getOption("fftempdir")} and only "./filename" in \code{\link{getwd}}.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{fftempfile}}, \code{\link{finalizer}}, \code{\link{ff}}, \code{\link{as.ff}}, \code{\link{as.ram}}, \code{\link{update.ff}},  \code{\link{file.move}}}
#! \examples{
#!   \dontrun{
#!   message("Neither giving pattern nor filename gives a random filename 
#! with extension ffextension in fftempdir")
#!   x <- ff(1:12)
#!   finalizer(x)
#!   filename(x)
#!   message("Giving a pattern with just a prefix moves to a random filename 
#! beginning with the prefix in fftempdir")
#!   pattern(x) <- "myprefix_"
#!   filename(x)
#!   message("Giving a pattern with a path and prefix moves to a random filename 
#! beginning with prefix in path (use . for getwd) ")
#!   pattern(x) <- "./myprefix"
#!   filename(x)
#!   message("Giving a filename moves to exactly this filename and extension 
#! in the R-expected place) ")
#!   if (!file.exists("./myfilename.myextension")){
#!     filename(x) <- "./myfilename.myextension"
#!     filename(x)
#!   }
#!
#!   message("NOTE that the finalizer has changed from 'delete' to 'close': 
#! now WE are responsible for deleting the file - NOT the finalizer")
#!   finalizer(x)
#!   delete(x)
#!   rm(x)
#!   }
#! }
#! \keyword{ IO }
#! \keyword{ data }

filename.ff_pointer <- function(x
, ... # dummy to keep R CMD check quiet
)
  attr(x, "filename")

filename.default <- function(x
, ... # dummy to keep R CMD check quiet
)
  attr(attr(x, "physical"), "filename")

"filename<-.ff" <- function(x
, ... # dummy to keep R CMD check quiet
, value
){
  isopen <- is.open(x)
  tmpdirnam <- getOption("fftempdir")

  oldnam <- filename(x)
  olddirnam <- dirname(oldnam)

  splitted <- splitPathFile(value)
  if (splitted$path=="."){
    splitted$path <- getwd()
    value <- unsplitPathFile(splitted)
  }else if(splitted$path=="" && splitted$fsep==""){
    splitted$path <- tmpdirnam
    value <- unsplitPathFile(splitted)
  }else{
    # convert to absolute path
    cwd <- getwd()
    on.exit(setwd(cwd))
    dfile <- dirname(value)
    bfile <- basename(value)
    setwd(dfile)
    dfile <- getwd()
    value <- file.path(dfile, bfile)
    # fix problem in file.path
    value <- gsub("/+","/",value)
  }

  if (olddirnam==tmpdirnam){
    if (splitted$path!=tmpdirnam){
      if (isopen)
        attr(attr(x,"physical"),"finalizer") <- "close"
      else
        attr(attr(x,"physical"),"finalizer") <- NULL
    }
  }else{
    if (splitted$path==tmpdirnam)
      finalizer(x) <- "delete"
  }

  if (isopen){
    close(x)
    on.exit(open(x), add=TRUE)
  }

  if(file.move(oldnam, value))
    physical(x)$filename <- value
  else
    stop("changing ff filename from '", oldnam, "' to '", value, "' failed")
  x
}


pattern.ff <- function(x, ...){
  attr(attr(x,"physical"),"pattern")
}

"pattern<-.ff" <- function(x, ...,value){
  filename <- fftempfile(value)
  filename(x) <- filename
  x
}


filename.ffdf <- function(x, ...)
  lapply(physical(x), filename)


"pattern<-.ffdf" <- function(x, ..., value){
  for (i in seq_len(ncol(x)))
    pattern(x[[i]]) <- value
  x
}


#! \name{is.readonly}
#! \alias{is.readonly}
#! \alias{is.readonly.ff}
#! \title{ Get readonly status }
#! \description{
#!   Get readonly status of an ff object
#! }
#! \usage{
#! is.readonly(x, \dots)
#! \method{is.readonly}{ff}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ \code{x} }
#!   \item{\dots}{ \code{\dots} }
#! }
#! \details{
#!   ff objects can be created/opened with \code{readonly=TRUE}.
#!   After each opening of the ff file readonly status is stored in the \code{\link[=physical.ff]{physical}} attributes and serves as the default for the next opening.
#!   Thus querying a closed ff object gives the last readonly status.
#! }
#! \value{
#!   logical scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{open.ff}}, \code{\link[=physical.ff]{physical}} }
#! \examples{
#!   x <- ff(1:12)
#!   is.readonly(x)
#!   close(x)
#!   open(x, readonly=TRUE)
#!   is.readonly(x)
#!   close(x)
#!   is.readonly(x)
#!   rm(x)
#! }
#! \keyword{ IO }
#! \keyword{ data }


is.readonly.ff <- function(x
, ... # dummy to keep R CMD check quiet
){
  attr(attr(x, "physical"), "readonly")
}


#! \name{is.open}
#! \alias{is.open}
#! \alias{is.open.ff}
#! \alias{is.open.ffdf}
#! \alias{is.open.ff_pointer}
#! \title{ Test if object is opened }
#! \description{
#!   Test whether an ff or ffdf object or a \code{ff_pointer} is opened.
#! }
#! \usage{
#! is.open(x, \dots)
#! \method{is.open}{ff}(x, \dots)
#! \method{is.open}{ffdf}(x, \dots)
#! \method{is.open}{ff_pointer}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an \code{\link{ff}} or \code{\link{ffdf}} object }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!   ff objects open automatically if accessed while closed.
#!   For ffdf objects we test all of their \code{\link[=physical.ffdf]{physical}} components including their \code{\link[=row.names.ffdf]{row.names}} if they are \code{\link{is.ff}}
#! }
#! \value{
#!   TRUE or FALSE (or NA if not all components of an ffdf object are opened or closed)
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{is.readonly}}, \code{\link{open.ff}}, \code{\link{close.ff}} }
#! \examples{
#!   x <- ff(1:12)
#!   is.open(x)
#!   close(x)
#!   is.open(x)
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

is.open.ff_pointer <- function(x
, ... # dummy to keep R CMD check quiet
){
  .Call("is_open", x, PACKAGE="ff")
}

is.open.ff <- function(x
, ... # dummy to keep R CMD check quiet
){
  .Call("is_open", attr(x, "physical"), PACKAGE="ff")
}


#! \name{pagesize}
#! \alias{pagesize}
#! \alias{pagesize.ff}
#! \title{ Pagesize of ff object }
#! \description{
#!   Returns current pagesize of ff object
#! }
#! \usage{
#! pagesize(x, \dots)
#! \method{pagesize}{ff}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an \code{\link{ff}} object }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \value{
#!   integer number of bytes
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{getpagesize}} }
#! \examples{
#!   x <- ff(1:12)
#!   pagesize(x)
#! }
#! \keyword{ IO }
#! \keyword{ data }



pagesize.ff <- function(x, ...){
  attr(attr(x, "physical"), "pagesize")
}


#! \name{maxlength}
#! \alias{maxlength}
#! \alias{maxlength.ff}
#! \alias{maxlength.default}
#! \title{ Get physical length of an ff or ram object }
#! \description{
#!   \command{maxlength} returns the physical length of an ff or ram object
#! }
#! \usage{
#! maxlength(x, \dots)
#! \method{maxlength}{ff}(x, \dots)
#! \method{maxlength}{default}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ ff or ram object }
#!   \item{\dots}{ additional arguments (not used) }
#! }
#! \value{
#!   integer scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{length.ff}}, \code{\link{maxindex}} }
#! \examples{
#!   x <- ff(1:12)
#!   length(x) <- 10
#!   length(x)
#!   maxlength(x)
#!   x
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }




# since we have separated length from maxlength, we might allow length to be reduced (while maxlength remains the same)
maxlength.ff  <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  #we no longer call .Call("maxlength", x, PACKAGE="ff")
  attr(attr(x, "physical"), "maxlength")
}
maxlength.default <- function(x
, ... # dummy to keep R CMD check quiet
){
  length(x)
}


#! \name{symmetric}
#! \alias{symmetric}
#! \alias{symmetric.ff}
#! \alias{symmetric.default}
#! \alias{symmetric.dist}
#! \title{ Test for symmetric structure }
#! \description{
#!   Check if an object is inherently symmetric (its structure, not its data)
#! }
#! \usage{
#! symmetric(x, \dots)
#! \method{symmetric}{ff}(x, \dots)
#! \method{symmetric}{default}(x, \dots)
#! \method{symmetric}{dist}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an ff or ram object }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!   ff matrices can be declared symmetric at creation time. Compatibility function \command{symmetric.default} returns FALSE, \command{symmetric.dist} returns TRUE.
#! }
#! \value{
#!   TRUE or FALSE
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{symmetric}}, \code{\link{ff}} %, \code{\link{symm}}
#!         , \code{\link{dist}}, \code{\link{isSymmetric}} }
#! \examples{
#!   symmetric(matrix(1:16, 4, 4))
#!   symmetric(dist(rnorm(1:4)))
#! }
#! \keyword{ IO }
#! \keyword{ data }

symmetric.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  attr(attr(x, "virtual"),"Symmetric")
}
symmetric.default <- function(x
, ... # dummy to keep R CMD check quiet
)FALSE
symmetric.dist <- function(x
, ... # dummy to keep R CMD check quiet
)TRUE



#! \name{fixdiag}
#! \alias{fixdiag}
#! \alias{fixdiag<-}
#! \alias{fixdiag.ff}
#! \alias{fixdiag.default}
#! \alias{fixdiag.dist}
#! \title{ Test for fixed diagonal }
#! \description{
#!   Check if an object has fixed diagonal
#! }
#! \usage{
#! fixdiag(x, \dots)
#! fixdiag(x, \dots) <- value
#! \method{fixdiag}{ff}(x, \dots)
#! \method{fixdiag}{default}(x, \dots)
#! \method{fixdiag}{dist}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an ff or ram object }
#!   \item{value}{ assignement value }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!   ff symmetric matrices can be declared to have fixed diagonal at creation time. Compatibility function \command{fixdiag.default} returns NULL, \command{fixdiag.dist} returns 0.
#! }
#! \value{
#!   NULL or the scalar representing the fixed diagonal
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{fixdiag}}, \code{\link{ff}} %, \code{\link{symm}}
#!         , \code{\link{dist}} }
#! \examples{
#!   fixdiag(matrix(1:16, 4, 4))
#!   fixdiag(dist(rnorm(1:4)))
#! }
#! \keyword{ IO }
#! \keyword{ data }


fixdiag.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  attr(attr(x, "virtual"),"Fixdiag")
}
fixdiag.default <- function(x
, ... # dummy to keep R CMD check quiet
)
  NULL
fixdiag.dist <- function(x
, ... # dummy to keep R CMD check quiet
)
  0

#! \name{is.sorted}
#! \alias{is.sorted.default}
#! \alias{is.sorted<-.default}
#! \title{ Getting and setting 'is.sorted' physical attribute }
#! \description{
#!   Functions to mark an ff or ram object as 'is.sorted' and query this. Responsibility to maintain this attribute is with the user.
#! }
#! \usage{
#! \method{is.sorted}{default}(x, \dots)
#! \method{is.sorted}{default}(x, \dots) <- value
#! }
#! \arguments{
#!   \item{x}{ an ff or ram object }
#!   \item{\dots}{ ignored }
#!   \item{value}{ NULL (to remove the 'is.sorted' attribute) or TRUE or FALSE }
#! }
#! \details{
#!   Sorting is slow, see \code{\link[base]{sort}}.
#!   Checking whether an object is sorted can avoid unnessary sorting -- see \code{\link[base]{is.unsorted}}, \code{\link[bit]{intisasc}} -- but still takes too much time with large objects stored on disk.
#!   Thus it makes sense to maintain an attribute, that tells us whether sorting can be skipped.
#!   Note that -- though you change it yourself -- \code{is.sorted} is a \code{\link[=physical.ff]{physical}} attribute of an object,
#!   because it represents an attribute of the \emph{data}, which is shared between different \code{\link[=physical.ff]{virtual}} views of the object.
#! }
#! \value{
#!   TRUE (if set to TRUE) or FALSE (if set to NULL or FALSE)
#! }
#! \note{
#!   \command{ff} will set \code{is.sorted(x) <- FALSE} if \code{\link{clone}} or \code{\link{length<-.ff}} have increased length.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{is.ordered.ff}} for testing factor levels, \code{\link[base]{is.unsorted}} for testing the data, \code{\link[bit]{intisasc}} for a quick version thereof, \code{\link{na.count}} for yet another \code{\link[=physical.ff]{physical}} attribute }
#! \examples{
#!   x <- 1:12
#!   is.sorted(x) <- !( is.na(is.unsorted(x)) || is.unsorted(x))
#!   is.sorted(x)
#!   x[1] <- 100L
#!   message("don't forget to maintain once it's no longer TRUE")
#!   is.sorted(x) <- FALSE
#!   message("check whether as 'is.sorted' attribute is maintained")
#!   !is.null(physical(x)$is.sorted)
#!   message("remove the 'is.sorted' attribute")
#!   is.sorted(x) <- NULL
#!   message("NOTE that querying 'is.sorted' still returns FALSE")
#!   is.sorted(x)
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ arith }


# !is.sorted does not imply is.unsorted
is.sorted.default <- function(x, ...)
{
  s <- physical(x)$is.sorted
  if (is.null(s) || !s)
    FALSE
  else
    TRUE
}
"is.sorted<-.default" <- function(x
, ...
, value
)
{
  if (is.null(value))
    physical(x)$is.sorted <- NULL
  else if (is.na(value) || !value)
    physical(x)$is.sorted <- FALSE
  else
    physical(x)$is.sorted <- TRUE
  x
}



#! \name{na.count}
#! \alias{na.count.ff}
#! \alias{na.count.default}
#! \alias{na.count<-.ff}
#! \alias{na.count<-.default}
#! \title{ Getting and setting 'na.count' physical attribute }
#! \description{
#!   The 'na.count' physical attribute gives the current number of NAs if properly initialized and properly maintained, see details.
#! }
#! \usage{
#! \method{na.count}{ff}(x, \dots)
#! \method{na.count}{default}(x, \dots)
#! \method{na.count}{ff}(x, \dots) <- value
#! \method{na.count}{default}(x, \dots) <- value
#! }
#! \arguments{
#!   \item{x}{ an ff or ram object }
#!   \item{\dots}{ further arguments (not used) }
#!   \item{value}{ NULL (to remove the 'na.count' attribute) or TRUE to activate or an integer value }
#! }
#! \details{
#!   The 'na.count' feature is activated by assigning the current number of NAs to \code{na.count(x) <- currentNA} and deactivated by assigning NULL.
#!   The 'na.count' feature is maintained by the, \code{\link{getset.ff}}, \code{\link{readwrite.ff}} and \code{\link{swap}},
#!   other ff methods for writing -- \code{\link{set.ff}}, \code{\link{[[<-.ff}}, \code{\link{write.ff}}, \code{\link{[<-.ff}} -- will stop if 'na.count' is activated.
#!   The functions \command{na.count} and \command{na.count<-} are generic.
#!   For ram objects, the default method for \command{na.count} calculates the number of NAs on the fly, thus no maintenance restrictions apply.
#! }
#! \value{
#!   NA (if set to NULL or NA) or an integer value otherwise
#! }
#! \author{ Jens Oehlschlägel, Daniel Adler (C++ back-end) }
#! \seealso{ \code{\link{getset.ff}}, \code{\link{readwrite.ff}} and \code{\link{swap}} for methods that support maintenance of 'na.count', \code{\link[base]{NA}}, \code{\link{is.sorted}} for yet another \code{\link[=physical.ff]{physical}} attribute }
#! \examples{
#!   message("--- ff examples ---")
#!   x <- ff(1:12)
#!   na.count(x)
#!   message("activate the 'na.count' physical attribute and set the current na.count manually")
#!   na.count(x) <- 0L
#!   message("add one NA with a method that maintains na.count")
#!   swap(x, NA, 1)
#!   na.count(x)
#!   message("remove the 'na.count' physical attribute (and stop automatic maintenance)")
#!   na.count(x) <- NULL
#!   message("activate the 'na.count' physical attribute and have ff automatically 
#! calculate the current na.count")
#!   na.count(x) <- TRUE
#!   na.count(x)
#!   message("--- ram examples ---")
#!   x <- 1:12
#!   na.count(x)
#!   x[1] <- NA
#!   message("activate the 'na.count' physical attribute and have R automatically 
#! calculate the current na.count")
#!   na.count(x) <- TRUE
#!   na.count(x)
#!   message("remove the 'na.count' physical attribute (and stop automatic maintenance)")
#!   na.count(x) <- NULL
#!   na.count(x)
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

na.count.ff <- function(x
, ... # dummy to keep R CMD check quiet
){
  physical <- physical(x)
  if (!is.null(physical$VW))
    return(as.integer(NA))
  nc <- physical$na.count
  if (is.null(nc)){
    as.integer(NA)
  }else{
    nc
  }
}
"na.count<-.ff" <- function(x
, ...     # dummy to keep R CMD check quiet
, value   # assigning NA deactivates na.count, assigning integer value sets na.count to this value, assigning TRUE calculates na.count and sets it.
){
  if (!is.null(attr(attr(x, "physical"),"VW")))
    stop("you can't set na.count during existence of a virtual window, see ?vw")
  if (is.null(value))
    physical(x)$na.count <- NULL
  else{
    if (is.logical(value)){
      if (value){
        i1 <- i2 <- 0L  # dummy assignment to shut up R CMD CHECK NOTE about no visible binding
        value <- ffvecapply(sum(is.na(x[i1:i2])), X=x, RETURN=TRUE, CFUN="sum")
      }else{
        stop("assign value=NULL (not value=FALSE) to deactivate na.count")
      }
    }
    physical(x)$na.count <- as.integer(value)
  }
  x
}
# for ram objects we always calculate it on-the-fly (when 'activated')
na.count.default <- function(x
, ... # dummy to keep R CMD check quiet
){
  if (is.null(physical(x)$na.count))
    as.integer(NA)
  else
    as.integer(sum(is.na(x[])))
}
"na.count<-.default" <- function(x
, ... # dummy to keep R CMD check quiet
, value
){
  if (is.null(value))
    physical(x)$na.count <- NULL
  else
    physical(x)$na.count <- TRUE  # == 'activation' of na.count
  x
}


#! \name{physical.ff}
#! \alias{physical.ff}
#! \alias{physical<-.ff}
#! \alias{virtual.ff}
#! \alias{virtual<-.ff}
#! \title{ Getting and setting physical and virtual attributes of ff objects }
#! \description{
#!   Functions for getting and setting physical and virtual attributes of ff objects.
#! }
#! \usage{
#! \method{physical}{ff}(x)
#! \method{virtual}{ff}(x)
#! \method{physical}{ff}(x) <- value
#! \method{virtual}{ff}(x) <- value
#! }
#! \arguments{
#!   \item{x}{ an ff object }
#!   \item{value}{ a list with named elements }
#! }
#! \details{
#!   ff objects have physical and virtual attributes, which have different copying semantics:
#!   physical attributes are shared between copies of ff objects while virtual attributes might differ between copies.
#!   \code{\link{as.ram}} will retain some physical and virtual atrributes in the ram clone,
#!   such that \code{\link{as.ff}} can restore an ff object with the same attributes.
#! }
#! \value{
#!   \command{physical} and \command{virtual} returns a list with named elements
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{
#!  \code{\link[bit]{physical}}, \code{\link{physical.ffdf}}, \code{\link{ff}}, \code{\link{as.ram}}; \cr
#!  \code{\link{is.sorted}} and \code{\link{na.count}} for applications of physical attributes; \cr
#!  \code{\link{levels.ff}} and \code{\link{ramattribs}} for applications of virtual attributes
#! }
#! \examples{
#!   x <- ff(1:12)
#!   x
#!   physical(x)
#!   virtual(x)
#!   y <- as.ram(x)
#!   physical(y)
#!   virtual(y)
#!   rm(x,y); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }


physical.ff <- function(x){
  p <- attributes(attr(x, "physical"))
  p <- p[is.na(match(names(p), "class"))]
  if (is.na(match("vmode", names(p))))
    c(list(vmode=vmode(x)), p)
  else
    p
}


"physical<-.ff" <- function(x, value){
  attributes(attr(x, "physical")) <- c(value, list(class="ff_pointer"))
  x
}



# -- now virtual attributes follow -- NOTE that the class of the ff object also has copy by value semantics, thus you can have two objects of different classes pointing to the same ff file

virtual.ff <- function(x){
  v <- attributes(attr(x, "virtual"))
  v[is.na(match(names(v), "class"))]
}

"virtual<-.ff" <- function(x, value){
  attributes(attr(x, "virtual")) <- c(value, list(class="virtual"))
  x
}



#! \name{ramattribs}
#! \alias{ramclass}
#! \alias{ramclass.ff}
#! \alias{ramclass.default}
#! \alias{ramclass_excludes}
#! \alias{ramattribs}
#! \alias{ramattribs.ff}
#! \alias{ramattribs.default}
#! \alias{ramattribs_excludes}
#! \title{ Get ramclass and ramattribs }
#! \description{
#!   Functions \command{ramclass} and \command{ramattribs} return the respective virtual attributes, that determine which class (and attributes) an ff object receives when subscripted (or coerced) to ram.
#! }
#! \usage{
#! ramclass(x, \dots)
#! \method{ramclass}{ff}(x, \dots)
#! \method{ramclass}{default}(x, \dots)
#! ramattribs(x, \dots)
#! \method{ramattribs}{ff}(x, \dots)
#! \method{ramattribs}{default}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ \code{x} }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!  \command{ramclass} and \command{ramattribs} provide a general mechanism to store atomic classes in ff objects,
#!   for example \code{\link{factor}} -- see \code{\link{levels.ff}} -- and \code{\link[base:DateTimeClasses]{POSIXct}}, see the example.
#! }
#! \value{
#!   \command{ramclass} returns a character vector with classnames and \command{ramattribs} returns a list with names elemens just like \code{\link[base]{attributes}}.
#!   The vectors \code{ramclass_excludes} and \code{ramattribs_excludes} name those attributes, which are not exported from ff to ram objects when using \code{\link{as.ram}}.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ff}}, \code{\link[=physical.ff]{virtual}}, \code{\link{as.ram}}, \code{\link{levels.ff}}, \code{\link{attributes}}, \code{\link[base]{DateTimeClasses}} }
#! \examples{
#!   x <- ff(as.POSIXct(as.POSIXlt(Sys.time(), "GMT")), length=12)
#!   x
#!   ramclass(x)
#!   ramattribs(x)
#!   class(x[])
#!   attributes(x[])
#!   virtual(x)$ramattribs$tzone = NULL
#!   attributes(x[])
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }
#! \keyword{ classes }

ramclass.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
  attr(attr(x, "virtual"), "ramclass")
ramattribs.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
  attr(attr(x, "virtual"), "ramattribs")

ramclass.default <- function(x
, ... # dummy to keep R CMD check quiet
){
  cl <- class(x)
  cl <- cl[is.na(match(cl, ramclass_excludes))]
  if (length(cl))
    cl
  else
    NULL
}
ramattribs.default <- function(x
, ... # dummy to keep R CMD check quiet
){
  a <- attributes(x)
  a <- a[is.na(match(names(a), ramattribs_excludes))]
  if (length(a))
    a
  else
    NULL
}


#! \name{length.ff}
#! \alias{length.ff}
#! \alias{length<-.ff}
#! \title{ Getting and setting length }
#! \description{
#!   Gets and sets length of ff objects.
#! }
#! \usage{
#! \method{length}{ff}(x)
#! \method{length}{ff}(x) <- value
#! }
#! \arguments{
#!   \item{x}{ object to query }
#!   \item{value}{ new object length }
#! }
#! \details{
#!   Changing the length of ff objects is only allowed if no \code{\link{vw}} is used.
#!   Changing the length of ff objects will remove any \code{\link{dim.ff}} and \code{\link{dimnames.ff}} attribute.
#!   Changing the length of ff objects will remove any \code{\link{na.count}} or \code{\link{is.sorted}} attribute and warn about this.
#!   New elements are usually zero, but it may depend on OS and filesystem what they really are.
#!   If you want standard R behaviour: filling with NA ,you need to do this yourself.
#!   As an exception to this rule, ff objects with \code{\link{names.ff}} will be filled with NA's automatically,
#!   and the length of the names will be adjusted (filled with position numbers where needed, which can easily consume a lot of RAM,
#!   therefore removing 'names' will help to faster increase length without RAM problems).
#! }
#! \note{
#!   Special care needs to be taken with regard ff objects that represent factors.
#!   For ff factors based on UNSIGNED \code{\link{vmode}s}, new values of zero are silently interpreted as the first factor level.
#!   For ff factors based on SIGNED \code{\link{vmode}s}, new values of zero result in illegal factor levels.
#!   See \code{\link{nrow<-}}.
#! }
#! \value{
#!   Integer scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{length}}, \code{\link{maxlength}}, \code{\link{file.resize}}, \code{\link[ff:dim.ff]{dim}}, \code{\link[=physical.ff]{virtual}} }
#! \examples{
#!   x <- ff(1:12)
#!   maxlength(x)
#!   length(x)
#!   length(x) <- 10
#!   maxlength(x)
#!   length(x)
#!   length(x) <- 16
#!   maxlength(x)
#!   length(x)
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }



length.ff <- function(x)
{
  attr(attr(x, "virtual"), "Length")
}


"length<-.ff" <- function(x, value)
{
  virtual <- virtual(x)
  if (!is.null(virtual$VW))
    stop("you can't set length during existence of a virtual window, see ?vw")
  if (!is.null(virtual$Dim))
    stop("you can't set length of arrays, either remove the dim attribute or use 'dim<-' instead")

  value <- as.integer(value)
  oldlen <- virtual$Length

  if (value==oldlen)
    return(x)


  vm <- vmode(x)
  newbytes <- ceiling(value * .ffbytes[vm])

  io <- is.open(x)
  if (io)
    close(x)
  file.resize( filename(x), newbytes )
  open(x)

  physical <- physical(x)
  physical$maxlength <- value
  if (!is.null(physical$na.count)){
    physical$na.count <- NULL
    warning("changing length(ff) removed attribute 'na.count'")
  }
  if (!is.null(physical$is.sorted) && value>oldlen){
    physical$is.sorted <- NULL
    warning("increasing length(ff) removed attribute 'is.sorted'")
  }
  physical(x) <- physical

  virtual$Length <- value

  if (!is.null(virtual$Names)){
    if (value>oldlen){
      if (!.vunsigned[vm] && vm!="raw")
        x[hi(oldlen+1L, value)] <- NA
      virtual$Names <- c(virtual$Names, (oldlen+1L):value)
    }else
      virtual$Names <- virtual$Names[1:value]
  }

  virtual(x) <- virtual

  if (!io)
    close(x)

  x
}


#! \name{levels.ff}
#! \alias{levels.ff}
#! \alias{levels<-.ff}
#! \alias{is.factor}
#! \alias{is.factor.default}
#! \alias{is.factor.ff}
#! \alias{is.ordered}
#! \alias{is.ordered.default}
#! \alias{is.ordered.ff}
#! \title{ Getting and setting factor levels }
#! \description{
#!   \code{levels.ff<-} sets factor levels, \code{levels.ff} gets factor levels
#! }
#! \usage{
#! \method{levels}{ff}(x)
#! \method{levels}{ff}(x) <- value
#!  is.factor(x)
#!  is.ordered(x)
#! \method{is.factor}{ff}(x)
#! \method{is.ordered}{ff}(x)
#! \method{is.factor}{default}(x)
#! \method{is.ordered}{default}(x)
#! }
#! \arguments{
#!   \item{x}{ an ff object }
#!   \item{value}{ the new factor levels, if NA is an allowed level it needs to be given explicitely, nothing is excluded }
#! }
#! \details{
#!   The ff object must have an integer vmode, see \code{\link{.rammode}}.
#!   If the mode is unsigned -- see \code{\link{.vunsigned}} -- the first factor level is coded with 0L instead of 1L in order to maximize the number of codable levels.
#!   Usually the internal ff coding -- see \code{\link{ram2ffcode}} -- is invisible to the user: when subscripting from an ff factor, unsigend codings are automatically converted to R's standard factor codes starting at 1L.
#!   However, you need to be aware of the internal ff coding in two situtations. \cr
#!   1. If you convert an ff integer object to an ff factor object and vice versa by assigning levels and \code{is.null(oldlevels)!=is.null(newlevels)}.  \cr
#!   2. Assigning data that does not match any level usually results in NA, however, in unsigned types there is no NA and all unknown data are mapped to the first level.
#! }
#! \value{
#!   \command{levels} returns a character vector of levels (possibly including \code{as.cha racter(NA)}).
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ When levels as assigned to an ff object that formerly had not levels, we assign automatically \code{\link{ramclass}} == "factor". If you want to change to an ordered factor, use \code{\link[=virtual.ff]{virtual}$ramclass <- c("ordered", "factor")} }
#! \seealso{ \code{\link{ramclass}}, \code{\link{factor}}, \code{\link[=physical.ff]{virtual}} }
#! \examples{
#!   message("--- create an ff factor including NA as last level")
#!   x <- ff("a", levels=c(letters, NA), length=99)
#!   message('    we expect a warning because "A" is an unknown level')
#!   x[] <- c("a", NA,"A")
#!   x
#!   levels(x)
#!
#!   message("--- create an ff ordered factor")
#!   x <- ff(letters, levels=letters, ramclass=c("ordered","factor"), length=260)
#!   x
#!   levels(x)
#!
#!   message("    make it a non-ordered factor")
#!   virtual(x)$ramclass <- "factor"
#!   x
#!   rm(x); gc()
#!
#!  \dontrun{
#!   message("--- create an unsigned quad factor")
#!   x <- ff(c("A","T","G","C"), levels=c("A","T","G","C"), vmode="quad", length=100)
#!   x
#!   message("  0:3 coding usually invisible to the user")
#!   unclass(x[1:4])
#!   message("    after removing levels, the 0:3 coding becomes visible to the user")
#!   message("    we expect a warning here")
#!   levels(x) <- NULL
#!   x[1:4]
#!   rm(x); gc()
#!  }
#!
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }

levels.ff <- function(x){
  attr(attr(x, "virtual"), "Levels")
}

"levels<-.ff" <- function(x, value){
  v <- attr(attr(x, "physical"), "vmode")
  if (.rammode[v]!="integer")
    stop("factors must be integer")
  if (length(value)>.vmax[v]+.vunsigned[v])
    stop("vmode '", v, "' can carry max ", .vmax[v]+.vunsigned[v], " levels")
  oldlevels <- attr(attr(x, "virtual"), "Levels")
  if (is.null(oldlevels)){
    if (!is.null(value)){
      if (.vunsigned[v])
        warning("assigning levels to unsigned vmode interprets 0 as first level !!")
      attr(attr(x, "virtual"), "Levels") <- value
      attr(attr(x, "virtual"), "ramclass") <- "factor" # make it factor, we don't know if it is ordered
    }
  }else{
    if (is.null(value)){
      attr(attr(x, "virtual"), "ramclass") <- NULL
      attr(attr(x, "virtual"), "Levels") <- NULL
      if (.vunsigned[v])
        warning("removing levels from unsigned vmode leaves first level as 0 !!")
    }else{
      if (length(value)<length(oldlevels))
        warning("lengths of levels was reduced")
      attr(attr(x, "virtual"), "Levels") <- value
    }
  }
  x
}


is.factor.ff <- function(x){
  ramclass <- attr(attr(x, "virtual"), "ramclass")
  !is.null(ramclass) && !is.na(match("factor",ramclass))
}

is.ordered.ff <- function(x){
  ramclass <- attr(attr(x, "virtual"), "ramclass")
  !is.null(ramclass) && !is.na(match("ordered",ramclass))
}


#! \name{names.ff}
#! \alias{names.ff}
#! \alias{names<-.ff}
#! \alias{names.ff_array}
#! \alias{names<-.ff_array}
#! \title{ Getting and setting names }
#! \description{
#!   For \code{ff_vector}s you can set names, though this is not recommended for large objects.
#! }
#! \usage{
#!   \method{names}{ff}(x)
#!   \method{names}{ff}(x) <- value
#!   \method{names}{ff_array}(x)
#!   \method{names}{ff_array}(x) <- value
#! }
#! \arguments{
#!   \item{x}{ a ff vector }
#!   \item{value}{ a character vector }
#! }
#! \details{
#!   If \code{\link{vw}} is set, \command{names.ff} returns the appropriate part of the names, but you can't set names while \command{vw} is set.
#!   \command{names.ff\_array}
#!   returns NULL and setting names for
#!   \code{ff_array}s is not allowed,
#!   but setting \code{\link[ff:dimnames.ff_array]{dimnames}} is.
#! }
#! \value{
#!   \command{names} returns a character vector (or NULL)
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{names}}, \code{\link{dimnames.ff_array}}, \code{\link{vw}}, \code{\link[=physical.ff]{virtual}} }
#! \examples{
#!   x <- ff(1:26, names=letters)
#!   names(x)
#!   names(x) <- LETTERS
#!   names(x)
#!   names(x) <- NULL
#!   names(x)
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

names.ff <- function(x){
  virtual <- attr(x, "virtual")
  vw <- attr(virtual, "VW")
  nam <- attr(virtual, "Names")
  if (is.null(vw)){
    nam
  }else{
    if (vw[2]>0){
      nam <- nam[(vw[1]+1L):(vw[1]+vw[2])]
      nam
    }else{
      character()
    }
  }
}

"names<-.ff" <- function(x, value){
  if (is.null(value)){
    attr(attr(x, "virtual"), "Names") <- NULL
  }else{
    virtual <- attr(x, "virtual")
    vw <- attr(virtual, "VW")
    if (is.null(vw)){
      if (length(value)!=length(x))
        stop("names attributes has wrong length")
      attr(attr(x, "virtual"), "Names") <- as.character(as.vector(value))
    }else{
      nam <- attr(virtual, "Names")
      if (is.null(nam))
        nam <- character(sum(vw))
      if (vw[2]>0)
        nam[(vw[1]+1L):(vw[1]+vw[2])] <- as.character(as.vector(value))
      attr(attr(x, "virtual"), "Names") <- nam
    }
  }
  x
}

"names.ff_array" <- function(x)
  NULL

"names<-.ff_array" <- function(x, value){
  if (!is.null(value))
    stop("assigning names (other than NULL) to ff_array not supported (not useful and very complicated due to dimorder)")
  else
    x
}


#! \name{dimnames.ff_array}
#! \alias{dimnames.ff}
#! \alias{dimnames.ff_array}
#! \alias{dimnames<-.ff_array}
#! \title{ Getting and setting dimnames }
#! \description{
#!   For \code{ff_array}s you can set dimnames.
#! }
#! \usage{
#!   \method{dimnames}{ff_array}(x)
#!   \method{dimnames}{ff_array}(x) <- value
#! }
#! \arguments{
#!   \item{x}{ a ff array (or matrix) }
#!   \item{value}{ a list with length(dim(x)) elements (either NULL of character vector of length of dimension }
#! }
#! \details{
#!   if \code{\link{vw}} is set, \command{dimnames.ff\_array} returns the appropriate part of the names, but you can't set \command{dimnames} while \command{vw} is set.
#!   \command{dimnames} returns NULL for \code{ff_vectors} and setting \code{dimnames} for \code{ff_vector} is not allowed, but setting \code{\link[ff:names.ff]{names}} is.
#! }
#! \value{
#!   \command{dimnames} returns a list, see \code{\link[base]{dimnames}}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{dimnames}}, \code{\link{dim.ff}} , \code{\link{names.ff}} , \code{\link{vw}}, \code{\link[=physical.ff]{virtual}} }
#! \examples{
#!   x <- ff(1:12, dim=c(3,4), dimnames=list(letters[1:3], LETTERS[1:4]))
#!   dimnames(x)
#!   dimnames(x) <- list(LETTERS[1:3], letters[1:4])
#!   dimnames(x)
#!   dimnames(x) <- NULL
#!   dimnames(x)
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

dimnames.ff_array <- function(x){
  vw <- attr(attr(x, "virtual"), "VW")
  if (is.null(vw))
    attr(attr(x, "virtual"), "Dimnames")
  else{
    vw <- vw[1,]
    d <- attr(attr(x, "virtual"),"Dim")
    ii <- 1:length(d)
    names(ii) <- names(d)
    dn <- attr(attr(x, "virtual"), "Dimnames")
    lapply(ii, function(i){
      dn[[i]][(vw[i]+1L):(vw[i]+d[i])]
    })
  }
}

"dimnames<-.ff_array" <- function(x, value){
  if (!is.null(value)){
    if (!is.null(attr(attr(x, "virtual"), "VW")))
      stop("vw must be NULL before you can set dimnames")
    if (!is.list(value))
      stop("dimnames must be NULL or list")
    ffdim <- dim(x)
    ndim <- length(ffdim)
    if (length(value)!=ndim)
      stop("length of dimnames does not match length of dim")
    dimcheck <- sapply(1:ndim, function(i){
      is.null(value[[i]]) || length(value[[i]]) == ffdim[i]
    })
    if (!all(dimcheck))
      stop("dimnames[[i]] is neither NULL nor matches dim[i]")
  }
  attr(attr(x, "virtual"), "Dimnames") <- value
  x
}



#! \name{dim.ff}
#! \alias{dim.ff}
#! \alias{dim.ffdf}
#! \alias{dim<-.ff}
#! \alias{dim<-.ffdf}
#! \alias{dimorder}
#! \alias{dimorder.default}
#! \alias{dimorder.ff_array}
#! \alias{dimorder.ffdf}
#! \alias{dimorder<-}
#! \alias{dimorder<-.ff_array}
#! \alias{dimorder<-.ffdf}
#! \title{ Getting and setting dim and dimorder }
#! \description{
#!   Assigning \code{dim} to an \code{ff_vector} changes it to an \code{ff_array}.
#!   Beyond that \code{dimorder} can be assigned to change from column-major order to row-major order or generalizations for higher order \code{ff_array}.
#! }
#! \usage{
#!   \method{dim}{ff}(x)
#!   \method{dim}{ffdf}(x)
#!   \method{dim}{ff}(x) <- value
#!   \method{dim}{ffdf}(x) <- value
#!    dimorder(x, \dots)
#!    dimorder(x, \dots) <- value
#!   \method{dimorder}{default}(x, \dots)
#!   \method{dimorder}{ff_array}(x, \dots)
#!   \method{dimorder}{ffdf}(x, \dots)
#!   \method{dimorder}{ff_array}(x, \dots) <- value
#!   \method{dimorder}{ffdf}(x, \dots) <- value  # just here to catch forbidden assignments
#! }
#! \arguments{
#!   \item{x}{ a ff object }
#!   \item{value}{ an appropriate integer vector }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!    \command{dim} and \command{dimorder} are \code{\link[=physical.ff]{virtual}} attributes. Thus two copies of an R ff object can point to the same file but interpret it differently.
#!    \command{dim} has the usual meaning, \command{dimorder} defines the dimension order of storage, i.e. \code{c(1,2)} corresponds to R's standard column-major order,
#!    \code{c(1,2)} corresponds to row-major order, and for higher dimensional arrays dimorder can also be used. Standard dimorder is \code{1:length(dim(x))}. \cr
#!    For \code{\link{ffdf}} \code{dim} returns the number of rows and virtual columns. With \code{dim<-.ffdf} only the number of rows can be changed. For convenience you can assign \code{NA} to the number of columns. \cr
#!    For \code{\link{ffdf}} the dimorder returns non-standard dimorder if any of its columns contains a ff object with non-standard dimorder (see \code{\link{dimorderStandard}})
#!    An even higher level of virtualization is available using virtual windows, see \code{\link{vw}}.
#! }
#! \note{
#!   \code{x[]} returns a matrix like \code{x[,]} and thus respects dimorder, while \code{x[i:j]} returns a vector and simply returns elements in the stored order.
#!   Check the corresponding example twice, in order to make sure you understand that for non-standard dimorder \code{x[1:length(x)]} is \emph{not the same} as \code{as.vector(x[])}.
#! }
#! \value{
#!   \command{names} returns a character vector (or NULL)
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{dim}}, \code{\link{dimnames.ff_array}}, \code{\link{dimorderStandard}}, \code{\link{vw}}, \code{\link[=physical.ff]{virtual}} }
#! \examples{
#!   x <- ff(1:12, dim=c(3,4), dimorder=c(2:1))
#!   y <- x
#!   dim(y) <- c(4,3)
#!   dimorder(y) <- c(1:2)
#!   x
#!   y
#!   x[]
#!   y[]
#!   x[,bydim=c(2,1)]
#!   y[,bydim=c(2,1)]
#!
#!   message("NOTE that x[] like x[,] returns a matrix (respects dimorder),")
#!   message("while x[1:12] returns a vector IN STORAGE ORDER")
#!   message("check the following examples twice to make sure you understand this")
#!   x[,]
#!   x[]
#!   as.vector(x[])
#!   x[1:12]
#!   rm(x,y); gc()
#!
#!   \dontshow{
#!     message("some regression test with regard to different dimorders")
#!     k <- 24
#!     d <- 3:5
#!     n <- prod(d)
#!     for (i in 1:k){
#!       a <- array(sample(n), dim=sample(d))
#!       x <- as.ff(a, dimorder=sample(seq_along(d)))
#!       if (!identical(a[1:n], x[1:n]))
#!         stop("error in caclulating access positions")
#!       if (!identical(a[1:dim(a)[1],,], x[1:dim(a)[1],,]))
#!         stop("error in caclulating access positions")
#!     }
#!     rm(x); gc()
#!   }
#!   \dontrun{
#!     message("some performance comparison between different dimorders")
#!     n <- 100
#!     m <- 100000
#!     a <- ff(1L,dim=c(n,m))
#!     b <- ff(1L,dim=c(n,m), dimorder=2:1)
#!     system.time(lapply(1:n, function(i)sum(a[i,])))
#!     system.time(lapply(1:n, function(i)sum(b[i,])))
#!     system.time(lapply(1:n, function(i){i<-(i-1)*(m/n)+1; sum(a[,i:(i+m/n-1)])}))
#!     system.time(lapply(1:n, function(i){i<-(i-1)*(m/n)+1; sum(b[,i:(i+m/n-1)])}))
#!     rm(a,b); gc()
#!   }
#! }
#! \keyword{ IO }
#! \keyword{ data }


dim.ff <- function(x)
{
  attr(attr(x, "virtual"),"Dim")
}


# Attention, assigning dim sets dimorder to 1:ndim
"dim<-.ff" <- function(
  x
, value
)
{
  virtual <- attr(x, "virtual")
  if (!is.null(attr(virtual, "VW")))
    stop("you can't set dim during existence of a virtual window, see ?vw")
  if (is.null(value)){
    attr(virtual, "Dim") <- NULL
    attr(virtual, "Dimnames") <- NULL
    attr(virtual, "Dimorder") <- NULL
    attr(x, "virtual") <- virtual
    class(x) <- c("ff_vector","ff")
  }else{
    value <- as.integer(value)
    d <- attr(virtual, "Dim")
    if (identical(d, value))
      return(x)
    n <- attr(virtual, "Length")
    nvalue <- prod(value)
    nvaluedim <- length(value)
    if (nvalue==n){
      dimorder <- 1:nvaluedim
    }else{
      # we allow to grow or shrink the fastest rotating dim (given dimorder)
      dimorder <- attr(virtual, "Dimorder")
      nd <- length(d)
      if ( nvaluedim==nd && ( nd==1 || all(d[dimorder][-nd]==value[dimorder][-nd]) ) ){
        dim(x) <- NULL
        length(x) <- nvalue
        dim(x) <- value
        dimorder(x) <- dimorder
        return(x)
      }else
        stop("you can only change the fastest rotating dim")
    }
    attr(virtual, "Dim") <- value
    attr(virtual, "Dimnames") <- NULL
    attr(virtual, "Dimorder") <- dimorder
    attr(virtual, "Names") <- NULL
    attr(x, "virtual") <- virtual
    if (nvaluedim==2)
      class(x) <- c("ff_matrix", "ff_array","ff")
    else
      class(x) <- c("ff_array","ff")
  }
  x
}


dimorder.ff_array <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  attr(attr(x, "virtual"),"Dimorder")
}
# Attention: you can not arbitrarily choose dimorder
# if you assign dimorder<-, dimorder MUST correspond to the physical layout
"dimorder<-.ff_array" <- function(x
, ... # dummy to keep R CMD check quiet
, value
)  # we assume you know what you do if you use this !!
{
  if (!is.null(attr(attr(x, "virtual"),"VW")))
    stop("you can't set dimorder during existence of a virtual window, see ?vw")
  do <- attr(attr(x, "virtual"),"Dimorder")
  value <- as.vector(as.integer(value))
  if (!identical(sort(do), sort(value)))
    stop("illegal dimorder, do you know what you are doing?")
  attr(attr(x, "virtual"),"Dimorder") <- value
  x
}

dimorder.default <- function(x
, ... # dummy to keep R CMD check quiet
){
  d <- dim(x)
  if (is.null(d))
    NULL
  else
    1:length(d)
}



#! \name{vw}
#! \alias{vw}
#! \alias{vw<-}
#! \alias{vw.ff}
#! \alias{vw.default}
#! \alias{vw<-.ff_vector}
#! \alias{vw<-.ff_array}
#! \title{ Getting and setting virtual windows }
#! \description{
#!   The virtual window \command{vw} function allows one to define a virtual window into an \code{ff_vector} or \code{ff_array}.
#!   The ff object will behave like a smaller array and it is mapped into the specified region of the complete array.
#!   This allows for example to execute recursive divide and conquer algorithms that work on parts of the full object,
#!   without the need to repeatedly create subfiles.
#! }
#! \usage{
#! vw(x, \dots)
#! vw(x, \dots) <- value
#! \method{vw}{ff}(x, \dots)
#! \method{vw}{default}(x, \dots)
#! \method{vw}{ff_vector}(x, \dots) <- value
#! \method{vw}{ff_array}(x, \dots) <- value
#! }
#! \arguments{
#!   \item{x}{ an \code{ff_vector} or \code{ff_array} }
#!   \item{\dots}{ further arguments (not used) }
#!   \item{value}{ a vector or matrix with an Offset, Window and Rest component, see details and examples }
#! }
#! \details{
#!   Each dimension of an ff array (or vector) is decomposed into three components, an invisible Offset, a visibe Window and an invisible Rest.
#!   For each dimension the sum of the vw components must match the dimension (or length).
#!   For an \code{ff_vector}, \code{vw} is simply a vector[1:3], for an array is is a \code{matrix[1:3,1:length(dim(x))]}.
#!   \code{vw} is a \code{\link[=physical.ff]{virtual}} attribute. \cr
#! }
#! \value{
#!   NULL or a vw specification, see details
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{length.ff}}, \code{\link{dim.ff}}, \code{\link[=physical.ff]{virtual}} }
#! \examples{
#!   x <- ff(1:26, names=letters)
#!   y <- x
#!   vw(x) <- c(0, 13, 13)
#!   vw(y) <- c(13, 13, 0)
#!   x
#!   y
#!   x[1] <- -1
#!   y[1] <- -2
#!   vw(x) <- NULL
#!   x[]
#!
#!   z <- ff(1:24, dim=c(4,6), dimnames=list(letters[1:4], LETTERS[1:6]))
#!   z
#!   vw(z) <- rbind(c(1,1), c(2,4), c(1,1))
#!   z
#!
#!   rm(x,y,z); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ array }


vw.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  attr(attr(x, "virtual"),"VW")
}

"vw<-.ff_vector" <- function(x
, ... # dummy to keep R CMD check quiet
, value
)
{
  virtual <- attr(x, "virtual")
  vw <- attr(virtual,"VW")
  if (is.null(value)){
    if (!is.null(vw)){
      # restore original length when disabling vw
      attr(virtual,"Length") <- sum(vw)
      attr(virtual,"VW") <- NULL
    }
  }else{
    value <- as.integer(value)
    if (length(value)!=3)
      stop("vw(ff_vector) must be vector[3]")
    if (is.null(vw)){
      if (sum(value)!=attr(virtual, "Length"))
        stop("sum(vw(ff_vector)) must equal length(ff_vector)")
    }else{
      if (sum(value)!=sum(vw))
        stop("sum(vw(ff_vector)) must equal length(ff_vector)")
    }
    attr(virtual,"Length") <- value[2]
    attr(virtual,"VW") <- value
    attr(attr(x, "physical"),"na.count") <- NULL
  }
  attr(x, "virtual") <- virtual
  x
}

"vw<-.ff_array" <- function(x
, ... # dummy to keep R CMD check quiet
, value
)
{
  virtual <- attr(x, "virtual")
  vw <- attr(virtual,"VW")
  if (is.null(value)){
    if (!is.null(vw)){
      # restore original length when disabling vw
      d <- as.integer(colSums(attr(attr(x, "virtual"),"VW")))
      attr(virtual,"Length") <- as.integer(prod(d))
      attr(virtual,"Dim") <- d
      attr(virtual,"VW") <- NULL
    }
  }else{
    d <- attr(virtual,"Dim")
    storage.mode(value) <- "integer"
    if (!identical(dim(value), c(3L, length(d))))
      stop("vw(ff_array) must be matrix[3, length(dim(x))]")
    dimnames(value) <- NULL
    if (is.null(vw)){
      if (!identical(as.integer(colSums(value)), d))
        stop("colSums(vw(ff_array)) must equal dim(x)")
    }else{
      if (!identical(colSums(value), colSums(vw)))
        stop("colSums(vw(ff_array)) must equal dim(x)")
    }
    attr(virtual,"Length") <- as.integer(prod(value[2,]))
    attr(virtual,"Dim") <- value[2,]
    attr(virtual,"VW") <- value
    attr(attr(x, "physical"),"na.count") <- NULL
  }
  attr(x, "virtual") <- virtual
  x
}

vw.default <- function(x
, ... # dummy to keep R CMD check quiet
){
  NULL
}


#! \name{print.ff}
#! \alias{print.ff}
#! \alias{print.ffdf}
#! \alias{print.ff_vector}
#! \alias{print.ff_matrix}
#! \alias{str.ff}
#! \alias{str.ffdf}
#! \title{ Print and str methods }
#! \description{
#!   printing ff objects and compactly showing their structure
#! }
#! \usage{
#! \method{print}{ff}(x, \dots)
#! \method{print}{ff_vector}(x, maxlength = 16, \dots)
#! \method{print}{ff_matrix}(x, maxdim = c(16, 16), \dots)
#! \method{str}{ff}(object, nest.lev=0, \dots)
#! \method{str}{ffdf}(object, nest.lev=0, \dots)
#! }
#! \arguments{
#!   \item{x}{ a ff object }
#!   \item{object}{ a ff object }
#!   \item{nest.lev}{ current nesting level in the recursive calls to str }
#!   \item{maxlength}{ max number of elements to print from an \code{ff_vector} }
#!   \item{maxdim}{ max number of elements to print from each dimension from an \code{ff_array} }
#!   \item{\dots}{ further arguments to print }
#! }
#! \details{
#!   The print methods just print a few exmplary elements from the beginning and end of the dimensions.
#! }
#! \value{
#!   \code{invisible()}
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ff}}, \code{\link{print}}, \code{\link{str}} }
#! \examples{
#!   x <- ff(1:10000)
#!   x
#!   print(x, maxlength=30)
#!   dim(x) <- c(100,100)
#!   x
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ print }


print.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  vw <- vw(x)
  d <- dim(x)
  l <- levels(x)
  cat("ff (", if (file.exists(filename(x))) {if (is.open(x)) "open" else "closed"} else "deleted", ") "
  , vmode(x)
  , if (is.null(d) && !is.null(vw))  paste(" offset=", vw[1], " length=", vw[2], " totallength=", sum(vw), sep="")
    else paste(" length=", length(x), sep="")
  , " (", maxlength(x), ")"
  , if (is.sorted(x)) " (sorted)"
  , if (!is.null(physical(x)$na.count)) paste(" (na.count=",na.count(x),")",sep="")
  , if (!is.null(d))
      paste(
        if (is.null(vw))
          paste(" dim=c(", paste(dim(x), collapse=","), ")", sep="")
        else
          paste(" offset=c(", paste(vw[1,], collapse=","), ") dim=c(", paste(vw[2,], collapse=","), ") totaldim=c(", paste(colSums(vw), collapse=","), ")", sep="")
        , if (!is.null(dimorder(x))) paste(" dimorder=c(", paste(dimorder(x), collapse=","), ")", sep="")
        , if (symmetric(x)) paste(" symmetric", if (!is.null(fixdiag(x))) paste(" (fixdiag=", fixdiag(x), ")", sep=""), sep="")
    ,  sep="")
  , if (!is.null(l)) paste(" levels: ", paste(ifelse(is.na(l), "<NA>", l), sep="", collapse=if (is.ordered(x)) " < " else " "), sep="")
  , "\n", sep="")
}

print.ff_vector <- function(x, maxlength=16, ...){
  NextMethod("print")
  if(is.open(x)){
    print(vecprint(x, maxlength=maxlength), ...)
  }
}

print.ff_matrix <- function(x, maxdim=c(16,16), ...){
  NextMethod("print")
  if(is.open(x))
    print(matprint(x, maxdim=maxdim), ...)
}


str.ff <- function(object, nest.lev=0, ...){
  nest.str <- paste(rep(" ..", nest.lev), collapse="")
  str(unclass(object), nest.lev=nest.lev, ...)
  cat(nest.str, ' - attr(*, "class") = ', sep="")
  str(class(object), nest.lev=nest.lev, ...)
}



# --- ff new / update / clone / open / close / delete / deleteIfOpen -----------------------------------------------------------

#! \name{ff}
#! \alias{ff}
#! \alias{ff_pointer}
#! \title{ ff classes for representing (large) atomic data }
#! \description{
#!   The ff package provides atomic data structures that are stored on disk but behave (almost) as if they were in RAM by
#!   mapping only a section (pagesize) into main memory (the effective main memory consumption per ff object).
#!   Several access optimization techniques such as Hyrid Index Preprocessing (\code{\link{as.hi}}, \code{\link{update.ff}}) and Virtualization (\code{\link[=physical.ff]{virtual}}, \code{\link{vt}}, \code{\link{vw}}) are implemented to achieve good performance even with large datasets.
#!   In addition to the basic access functions, the ff package also provides compatibility functions that facilitate writing code for ff and ram objects (\code{\link{clone}}, \code{\link{as.ff}}, \code{\link{as.ram}}) and very basic support for operating on ff objects (\code{\link{ffapply}}).
#!   While the (possibly packed) raw data is stored on a flat file, meta
#!   informations about the atomic data structure such as its dimension,
#!   virtual storage mode (\code{\link{vmode}}), factor level encoding,
#!   internal length etc.. are stored as an ordinary R object (external
#!   pointer plus attributes) and can be saved in the workspace.
#!   The raw flat file data encoding is always in native machine format for
#!   optimal performance and provides several packing schemes for different
#!   data types such as logical, raw, integer and double (in an extended version
#!   support for more tighly packed virtual data types is supported).
#!   flatfile data files can be shared among ff objects in the same R process or
#!   even from different R processes due to Memory-Mapping, although the
#!   caching effects have not been tested extensively.
#!   \cr
#!   Please do read and understand the limitations and warnings in \code{\link{LimWarn}} before you do anything serious with package ff.
#! }
#! \usage{
#! ff( initdata  = NULL
#! , length      = NULL
#! , levels      = NULL
#! , ordered     = NULL
#! , dim         = NULL
#! , dimorder    = NULL
#! , bydim       = NULL
#! , symmetric   = FALSE
#! , fixdiag     = NULL
#! , names       = NULL
#! , dimnames    = NULL
#! , ramclass    = NULL
#! , ramattribs  = NULL
#! , vmode       = NULL
#! , update      = NULL
#! , pattern     = NULL
#! , filename    = NULL
#! , overwrite   = FALSE
#! , readonly    = FALSE
#! , pagesize    = NULL  # getOption("ffpagesize")
#! , caching     = NULL  # getOption("ffcaching")
#! , finalizer   = NULL
#! , finonexit   = NULL  # getOption("fffinonexit")
#! , FF_RETURN   = TRUE
#! , BATCHSIZE   = .Machine$integer.max
#! , BATCHBYTES  = getOption("ffbatchbytes")
#! , VERBOSE     = FALSE
#! )
#! }
#! \arguments{
#!   \item{initdata}{ scalar or vector of the \code{\link{.vimplemented}} \code{\link{vmode}}s, recycled if needed, default 0, see also \code{\link{as.vmode}} and \code{\link{vector.vmode}} }
#!   \item{length}{ optional vector \code{\link[base]{length}} of the object (default: derive from 'initdata' or 'dim'), see \code{\link{length.ff}} }
#!   \item{levels}{ optional character vector of levels if (in this case initdata must be composed of these) (default: derive from initdata) }
#!   \item{ordered}{ indicate whether the levels are ordered (TRUE) or non-ordered factor (FALSE, default) }
#!   \item{dim}{ optional array \code{\link{dim}}, see \code{\link{dim.ff}} and \code{\link{array}} }
#!   \item{dimorder}{ physical layout (default 1:length(dim)), see \code{\link{dimorder}} and \code{\link{aperm}} }
#!   \item{bydim}{ dimorder by which to interpret the 'initdata', generalization of the 'byrow' paramter in \code{\link{matrix}} }
#!   \item{symmetric}{ extended feature: TRUE creates symmetric matrix (default FALSE) %, see \code{\link{symm}}, \code{\link{ff_symm}}, \code{\link{ff_dist}}
#!        }
#!   \item{fixdiag}{ extended feature: non-NULL scalar requires fixed diagonal for symmetric matrix (default NULL is free diagonal) }
#!   \item{names}{ NOT taken from initdata, see \code{\link{names}} }
#!   \item{dimnames}{ NOT taken from initdata, see \code{\link{dimnames}} }
#!   \item{ramclass}{ class attribute attached when moving all or parts of this ff into ram, see \code{\link{ramclass}} }
#!   \item{ramattribs}{ additional attributes attached when moving all or parts of this ff into ram, see \code{\link{ramattribs}} }
#!   \item{vmode}{ virtual storage mode (default: derive from 'initdata'), see \code{\link{vmode}} and \code{\link{as.vmode}} }
#!   \item{update}{ set to FALSE to avoid updating with 'initdata' (default TRUE) (used by \code{\link{ffdf}}) }
#!   \item{pattern}{ root pattern with or without path for automatic ff filename creation (default NULL translates to "ff"), see also argument 'filename' }
#!   \item{filename}{ ff \code{\link{filename}} with or without path (default tmpfile with 'pattern' prefix); without path the file is created in \code{getOption("fftempdir")}, with path '.' the file is created in \code{\link{getwd}}. Note that files created in \code{getOption("fftempdir")} have default finalizer "delete" while other files have default finalizer "close". See also arguments 'pattern' and 'finalizer' and \code{\link[=physical.ff]{physical}} }
#!   \item{overwrite}{ set to TRUE to allow overwriting existing files (default FALSE) }
#!   \item{readonly}{ set to TRUE to forbid writing to existing files }
#!   \item{pagesize}{ pagesize in bytes for the memory mapping (default from \code{getOptions("ffpagesize")} initialized by \code{\link{getdefaultpagesize}}), see also \code{\link[=physical.ff]{physical}} }
#!   \item{caching}{ caching scheme for the backend, currently 'mmnoflush' or 'mmeachflush' (flush mmpages at each swap, default from \code{getOptions("ffcaching")} initialized with 'mmeachflush'), see also \code{\link[=physical.ff]{physical}} }
#!   \item{finalizer}{ name of finalizer function called when ff object is \code{\link{remove}d} (default: ff files created in \code{getOptions("fftempdir")} are considered temporary and have default finalizer \code{\link[ff:delete.ff]{delete}}, files created in other locations have default finalizer \code{\link[ff:close.ff]{close}}); available finalizer generics are "close", "delete" and "deleteIfOpen", available methods are \code{\link{close.ff}}, \code{\link{delete.ff}} and \code{\link{deleteIfOpen.ff}}, see also argument 'finonexit' and \code{\link{finalizer}} }
#!   \item{finonexit}{ logical scalar determining whether  and \code{\link{finalize}} is also called when R is closed via \code{\link{q}}, (default TRUE from \code{getOptions("fffinonexit")}) }
#!   \item{FF_RETURN}{ logical scalar or ff object to be used. The default TRUE creates a new ff file. FALSE returns a ram object. Handing over an ff object here uses this or stops if not \code{\link{ffsuitable}} }
#!   \item{BATCHSIZE}{ integer scalar limiting the number of elements to be processed in \code{\link{update.ff}} when length(initdata)>1, default from \code{.Machine$integer.max} }
#!   \item{BATCHBYTES}{ integer scalar limiting the number of bytes to be processed in \code{\link{update.ff}} when length(initdata)>1, default from \code{getOption("ffbatchbytes")}, see also \code{\link{.rambytes}} }
#!   \item{VERBOSE}{ set to TRUE for verbosing in \code{\link{update.ff}} when length(initdata)>1, default FALSE }
#! }
#! \details{
#!  The atomic data is stored in \code{\link{filename}} as a native encoded raw flat file on disk, OS specific limitations of the file system apply.
#!  The number of elements per ff object is limited to the integer indexing, i.e. \code{\link{.Machine}$integer.max}.
#!  Atomic objects created with \command{ff} are \code{\link{is.open}}, a C++ object is ready to access the file via memory-mapping.
#!  Currently the C++ backend provides two caching schemes: 'mmnoflush' let the OS decide when to flash memory mapped pages
#!  and 'mmeachflush' will flush memory mapped pages at each page swap per ff file.
#!  These minimal memory ressources can be released by \code{\link[ff:close.ff]{close}ing} or \code{\link[ff:delete.ff]{delete}ing} the ff file.
#!  ff objects can be \code{\link{save}d} and \code{\link{load}ed} across R sessions. If the ff file still exists in the same location,
#!  it will be \code{\link[ff:open.ff]{open}ed} automatically at the first attempt to access its data. If the ff object is \code{\link{remove}d},
#!  at the next garbage collection (see \code{\link{gc}}) the ff object's \code{\link{finalizer}} is invoked.
#!  Raw data files can be made accessible as an ff object by explicitly given the filename and vmode but no size information (length or dim).
#!  The ff object will open the file and handle the data with respect to the given vmode.
#!  The \code{\link[ff:close.ff]{close}} finalizer will close the ff file, the \code{\link[ff:delete.ff]{delete}} finalizer will delete the ff file.
#!  The default finalizer \code{\link{deleteIfOpen}} will delete open files and do nothing for closed files. If the default finalizer is used,
#!  two actions are needed to protect the ff file against deletion: create the file outside the standard 'fftempdir' and close the ff object before removing it or before quitting R.
#!  When R is exited through \code{\link{q}}, the finalizer will be invoked depending on the 'fffinonexit' option, furthermore the 'fftempdir' is \code{\link{unlink}ed}. \cr
#! }
#! \value{
#!   If (\code{!FF_RETURN}) then a ram object like those generated by \code{\link{vector}}, \code{\link{matrix}}, \code{\link{array}} but with attributes 'vmode', 'physical' and 'virtual' accessible via \code{\link{vmode}}, \code{\link[=physical.ff]{physical}} and \code{\link[=physical.ff]{virtual}}  \cr
#!   If (\code{FF_RETURN}) an object of class 'ff' which is a a list with two components:
#!   \item{physical}{an external pointer of class '\code{ff_pointer}' which carries attributes with copy by reference semantics: changing a physical attribute of a copy changes the original }
#!   \item{virtual}{an empty list which carries attributes with copy by value semantics: changing a virtual attribute of a copy does not change the original }
#! }
#! \section{Physical object component}{
#!   The '\code{ff_pointer}' carries the following 'physical' or readonly attributes, which are accessible via \code{\link[=physical.ff]{physical}}:
#!  \tabular{rl}{
#!   \code{vmode    } \tab see \code{\link{vmode}} \cr
#!   \code{maxlength} \tab see \code{\link{maxlength}} \cr
#!   \code{pattern  } \tab see parameter 'pattern' \cr
#!   \code{filename } \tab see \code{\link{filename}} \cr
#!   \code{pagesize } \tab see parameter 'pagesize' \cr
#!   \code{caching  } \tab see parameter 'caching' \cr
#!   \code{finalizer} \tab see parameter 'finalizer' \cr
#!   \code{finonexit} \tab see parameter 'finonexit' \cr
#!   \code{readonly } \tab see \code{\link{is.readonly}} \cr
#!   \code{class    } \tab The external pointer needs class 'ff\_pointer' to allow method dispatch of finalizers  \cr
#!  }
#! }
#! \section{Virtual object component}{
#!   The 'virtual' component carries the following attributes (some of which might be NULL):
#!  \tabular{rl}{
#!   \code{Length    } \tab see \code{\link{length.ff}} \cr
#!   \code{Levels    } \tab see \code{\link{levels.ff}} \cr
#!   \code{Names     } \tab see \code{\link{names.ff}} \cr
#!   \code{VW        } \tab see \code{\link{vw.ff}} \cr
#!   \code{Dim       } \tab see \code{\link{dim.ff}} \cr
#!   \code{Dimorder  } \tab see \code{\link{dimorder}} \cr
#!   \code{Symmetric } \tab see \code{\link{symmetric.ff}} \cr
#!   \code{Fixdiag   } \tab see \code{\link{fixdiag.ff}} \cr
#!   \code{ramclass  } \tab see \code{\link{ramclass}} \cr
#!   \code{ramattribs} \tab see \code{\link{ramattribs}} \cr
#!  }
#! }
#! \section{Class}{
#!   You should not rely on the internal structure of ff objects or their ram versions. Instead use the accessor functions like \code{\link{vmode}}, \code{\link[=physical.ff]{physical}} and \code{\link[=physical.ff]{virtual}}.
#!   Still it would be wise to avoid attributes AND classes 'vmode', 'physical' and 'virtual' in any other packages.
#!   Note that the 'ff' object's class attribute also has copy-by-value semantics ('virtual').
#!   For the 'ff' object the following class attritibutes are known:
#!  \tabular{rl}{
#!   vector \tab \code{c("ff_vector","ff")} \cr
#!   matrix \tab \code{c("ff_matrix","ff_array","ff")} \cr
#!   array \tab \code{c("ff_array","ff")} \cr
#!   symmetric matrix \tab \code{c("ff_symm","ff")} \cr
#!   distance matrix \tab \code{c("ff_dist","ff_symm","ff")} \cr
#!   reserved for future use \tab \code{c("ff_mixed","ff")} \cr
#!  }
#! }
#! \section{Methods}{
#!  The following methods and functions are available for ff objects:
#!  \tabular{lrll}{
#!   \emph{ Type} \tab  \emph{ Name }  \tab \emph{ Assign }  \tab \emph{Comment}  \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Basic functions}  \cr
#!   function \tab  \code{\link{ff}}                         \tab \emph{ }  \tab constructor for ff and ram objects \cr
#!   generic  \tab  \code{\link[ff:update.ff]{update}}       \tab \emph{ }  \tab updates one ff object with the content of another \cr
#!   generic  \tab  \code{\link{clone}}                      \tab \emph{ }  \tab clones an ff object optionally changing some of its features \cr
#!   method   \tab  \code{\link[ff:print.ff]{print}}         \tab \emph{ }  \tab print ff \cr
#!   method   \tab  \code{\link[ff:str.ff]{str}}             \tab \emph{ }  \tab ff object structure \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Class test and coercion}  \cr
#!   function \tab  \code{\link{is.ff}}                      \tab \emph{ }  \tab check if inherits from ff \cr
#!   generic  \tab  \code{\link{as.ff}}                      \tab \emph{ }  \tab coerce to ff, if not yet \cr
#!   generic  \tab  \code{\link{as.ram}}                     \tab \emph{ }  \tab coerce to ram retaining some of the ff information \cr
#!   generic  \tab  \code{\link[=as.bit.ff]{as.bit}}          \tab \emph{ }  \tab coerce to \code{\link[bit]{bit}} \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Virtual storage mode} \cr
#!   generic  \tab  \code{\link{vmode}}                      \tab \code{<-} \tab get and set virtual mode (setting only for ram, not for ff objects) \cr
#!   generic  \tab  \code{\link{as.vmode}}                   \tab \emph{ }  \tab coerce to vmode (only for ram, not for ff objects) \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Physical attributes}  \cr
#!   function \tab  \code{\link[=physical.ff]{physical}}                   \tab \code{<-} \tab set and get physical attributes \cr
#!   generic  \tab  \code{\link{filename}}                   \tab \emph{<-}  \tab get and set filename \cr
#!   generic  \tab  \code{\link{pattern}}                    \tab \emph{<-}  \tab get pattern and set filename path and prefix via pattern \cr
#!   generic  \tab  \code{\link{maxlength}}                  \tab \emph{ }  \tab get maxlength \cr
#!   generic  \tab  \code{\link{is.sorted}}                  \tab \code{<-} \tab set and get if is marked as sorted \cr
#!   generic  \tab  \code{\link{na.count}}                   \tab \code{<-} \tab set and get NA count, if set to non-NA only swap methods can change and na.count is maintained automatically \cr
#!   generic  \tab  \code{\link{is.readonly}}                \tab \emph{ }   \tab get if is readonly \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }   \tab \bold{Virtual attributes} \cr
#!   function \tab  \code{\link[=physical.ff]{virtual}}                    \tab \code{<-} \tab set and get virtual attributes \cr
#!   method   \tab  \code{\link[ff:length.ff]{length}}       \tab \code{<-} \tab set and get length \cr
#!   method   \tab  \code{\link[ff:dim.ff]{dim}}             \tab \code{<-} \tab set and get dim \cr
#!   generic  \tab  \code{\link{dimorder}}                   \tab \code{<-} \tab set and get the order of dimension interpretation \cr
#!   generic  \tab  \code{\link{vt}}                         \tab \code{}   \tab virtually transpose ff_array \cr
#!   method  \tab   \code{\link[ff:t.ff]{t}}                 \tab \code{}   \tab create transposed clone of ff_array \cr
#!   generic  \tab  \code{\link{vw}}                         \tab \code{<-} \tab set and get virtual windows \cr
#!   method   \tab  \code{\link[ff:names.ff]{names}}         \tab \code{<-} \tab set and get names \cr
#!   method   \tab  \code{\link[ff:dimnames.ff]{dimnames}}   \tab \code{<-} \tab set and get dimnames \cr
#!   generic  \tab  \code{\link{symmetric}}                  \tab \emph{ }   \tab get if is symmetric \cr
#!   generic  \tab  \code{\link{fixdiag}}                    \tab \code{<-} \tab set and get fixed diagonal of symmetric matrix \cr
#!   method   \tab  \code{\link{levels}}                     \tab \code{<-} \tab levels of factor  \cr
#!   generic  \tab  \code{\link{recodeLevels}}               \tab \code{ }  \tab recode a factor to different levels \cr
#!   generic  \tab  \code{\link{sortLevels}}                 \tab \code{ }  \tab sort the levels and recoce a factor \cr
#!   method   \tab  \code{\link{is.factor}}                  \tab \emph{ }  \tab if is factor \cr
#!   method   \tab  \code{\link{is.ordered}}                 \tab \emph{ }  \tab if is ordered (factor) \cr
#!   generic  \tab  \code{\link{ramclass}}                   \tab \code{}   \tab get ramclass \cr
#!   generic  \tab  \code{\link{ramattribs}}                 \tab \code{}   \tab get ramattribs \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Access functions}  \cr
#!   function \tab  \code{\link{get.ff}}                     \tab \emph{ }  \tab get single ff element (currently \code{\link[ff:[[.ff]{[[}} is a shortcut) \cr
#!   function \tab  \code{\link{set.ff}}                     \tab \emph{ }  \tab set single ff element (currently \code{\link[ff:[[<-.ff]{[[<-}} is a shortcut) \cr
#!   function \tab  \code{\link{getset.ff}}                  \tab \emph{ }  \tab set single ff element and get old value in one access operation \cr
#!   function \tab  \code{\link{read.ff}}                    \tab \emph{ }  \tab get vector of contiguous elements \cr
#!   function \tab  \code{\link{write.ff}}                   \tab \emph{ }  \tab set vector of contiguous elements  \cr
#!   function \tab  \code{\link{readwrite.ff}}               \tab \emph{ }  \tab set vector of contiguous elements and get old values in one access operation \cr
#!   method   \tab  \code{\link[ff:[.ff]{[}}                 \tab \emph{ }  \tab get vector of indexed elements, uses HIP, see \code{\link{hi}} \cr
#!   method   \tab  \code{\link[ff:[<-.ff]{[<-}}             \tab \emph{ }  \tab set vector of indexed elements, uses HIP, see \code{\link{hi}} \cr
#!   generic  \tab  \code{\link[ff:swap.ff]{swap}}           \tab \emph{ }  \tab set vector of indexed elements and get old values in one access operation \cr
#!   generic  \tab  \code{\link[ff:add.ff]{add}}             \tab \emph{ }  \tab (almost) unifies '+=' operation for ff and ram objects \cr
#!   generic  \tab  \code{\link[ff:bigsample.ff]{bigsample}} \tab \emph{ }  \tab sample from ff object \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Opening/Closing/Deleting}                                             \cr
#!   generic  \tab  \code{\link{is.open}}                    \tab \emph{ }  \tab check if ff is open \cr
#!   method   \tab  \code{\link[ff:open.ff]{open}}           \tab \emph{ }  \tab open ff object (is done automatically on access) \cr
#!   method   \tab  \code{\link[ff:close.ff]{close}}         \tab \emph{ }  \tab close ff object (releases C++ memory and protects against file deletion if  \code{\link{deleteIfOpen}}) finalizer is used \cr
#!   generic  \tab  \code{\link[ff:delete.ff]{delete}}       \tab \emph{ }  \tab deletes ff file (unconditionally) \cr
#!   generic  \tab  \code{\link{deleteIfOpen}}               \tab \emph{ }  \tab deletes ff file if ff object is open (finalization method) \cr
#!   generic  \tab  \code{\link{finalizer}}                  \tab \emph{<-} \tab get and set finalizer \cr
#!   generic  \tab  \code{\link{finalize}}                   \tab \emph{}   \tab force finalization \cr
#!   \emph{ } \tab  \emph{ }                                 \tab \emph{ }  \tab \bold{Other}                                                     \cr
#!   function \tab  \code{\link{geterror.ff}}                \tab \emph{ }  \tab get error code \cr
#!   function \tab  \code{\link{geterrstr.ff}}               \tab \emph{ }  \tab get error message \cr
#!   }
#! }
#! \section{ff options}{
#!   Through \code{\link{options}} or \code{\link{getOption}} one can change and query global features of the ff package:
#!  \tabular{rll}{
#!   \emph{option}        \tab \emph{description}                           \tab \emph{default} \cr
#!   \code{fftempdir}     \tab default directory for creating ff files      \tab \code{\link{tempdir}} \cr
#!   \code{fffinalizer}   \tab name of default finalizer                    \tab \code{\link{deleteIfOpen}} \cr
#!   \code{fffinonexit}   \tab default for invoking finalizer on exit of R  \tab \code{TRUE} \cr
#!   \code{ffpagesize}    \tab default pagesize                             \tab \code{\link{getdefaultpagesize}} \cr
#!   \code{ffcaching}     \tab caching scheme for the C++ backend           \tab \code{'mmnoflush'} \cr
#!   \code{ffdrop}        \tab default for the \option{drop} parameter in the ff subscript methods  \tab TRUE \cr
#!   \code{ffbatchbytes}  \tab default for the byte limit in batched/chunked processing             \tab \code{\link{memory.limit}() \%/\% 100} \cr
#!  }
#! }
#! \section{OS specific}{
#!  The following table gives an overview of file size limits for common file systems (see \url{http://en.wikipedia.org/wiki/Comparison_of_file_systems} for further details):
#!  \tabular{ll}{
#!   \strong{File System} \tab \strong{File size limit} \cr
#!   FAT16              \tab 2GB \cr
#!   FAT32              \tab 4GB \cr
#!   NTFS               \tab 16GB \cr
#!   ext2/3/4           \tab 16GB to 2TB \cr
#!   ReiserFS           \tab 4GB (up to version 3.4) / 8TB (from version 3.5) \cr
#!   XFS                \tab 8EB \cr
#!   JFS                \tab 4PB \cr
#!   HFS                \tab 2GB \cr
#!   HFS Plus           \tab 16GB \cr
#!   USF1               \tab 4GB to 256TB \cr
#!   USF2               \tab 512GB to 32PB \cr
#!   UDF                \tab 16EB \cr
#!   }
#! }
#! \section{Credits}{
#!  Package Version 1.0
#!  \tabular{ll}{
#!   Daniel Adler       \tab \email{dadler@uni-goettingen.de} \cr
#!                      \tab R package design, C++ generic file vectors, Memory-Mapping, 64-bit Multi-Indexing adapter and Documentation, Platform ports \cr
#!   Oleg Nenadic       \tab \email{onenadi@uni-goettingen.de} \cr
#!                      \tab Index sequence packing, Documentation \cr
#!   Walter Zucchini    \tab \email{wzucchi@uni-goettingen.de} \cr
#!                      \tab Array Indexing, Sampling, Documentation \cr
#!   Christian Gläser   \tab \email{christian\_glaeser@gmx.de} \cr
#!                      \tab Wrapper for biglm package \cr
#!   }
#!  Package Version 2.0
#!  \tabular{ll}{
#!   Jens Oehlschlägel  \tab \email{Jens.Oehlschlaegel@truecluster.com} \cr
#!                      \tab R package redesign; Hybrid Index Preprocessing; transparent object creation and finalization; vmode design; virtualization and hybrid copying; arrays with dimorder and bydim; symmetric matrices; factors and POSIXct; virtual windows and transpose; new generics update, clone, swap, add, as.ff and as.ram; ffapply and collapsing functions. R-coding, C-coding and Rd-documentation. \cr
#!   Daniel Adler       \tab \email{dadler@uni-goettingen.de} \cr
#!                      \tab C++ generic file vectors, vmode implementation and low-level bit-packing/unpacking, arithmetic operations and NA handling, Memory-Mapping and backend caching. C++ coding and platform ports. R-code extensions for opening existing flat files readonly and shared. \cr
#!   }
#! }
#! \note{ Note that the standard finalizers are generic functions, their dispatch to the '\code{ff_pointer}' method happens at finalization time, their 'ff' methods exist for direct calling.
#! }
#! \section{Licence}{Package under GPL-2, included C++ code released by Daniel Adler under the less restrictive ISCL}
#! \seealso{ \code{\link{vector}}, \code{\link{matrix}}, \code{\link{array}}, \code{\link{as.ff}}, \code{\link{as.ram}} }
#! \examples{
#!   message("make sure you understand the following ff options 
#!     before you start using the ff package!!")
#!   oldoptions <- options(fffinalizer="deleteIfOpen", fffinonexit="TRUE", fftempdir=tempdir())
#!   message("an integer vector")
#!   ff(1:12)                  
#!   message("a double vector of length 12")
#!   ff(0, 12)
#!   message("a 2-bit logical vector of length 12 (vmode='boolean' has 1 bit)")
#!   ff(vmode="logical", length=12)
#!   message("an integer matrix 3x4 (standard colwise physical layout)")
#!   ff(1:12, dim=c(3,4))
#!   message("an integer matrix 3x4 (rowwise physical layout, but filled in standard colwise order)")
#!   ff(1:12, dim=c(3,4), dimorder=c(2,1))
#!   message("an integer matrix 3x4 (standard colwise physical layout, but filled in rowwise order
#! aka matrix(, byrow=TRUE))")
#!   ff(1:12, dim=c(3,4), bydim=c(2,1))
#!   gc()
#!   options(oldoptions)
#!
#!   if (ffxtensions()){
#!      message("a 26-dimensional boolean array using 1-bit representation
#!       (file size 8 MB compared to 256 MB int in ram)")
#!      a <- ff(vmode="boolean", dim=rep(2, 26))
#!      dimnames(a) <- dummy.dimnames(a)
#!      rm(a); gc()
#!   }
#!
#!   \dontrun{
#!
#!      message("This 2GB biglm example can take long, you might want to change
#!        the size in order to define a size appropriate for your computer")
#!      require(biglm)
#!
#!      b <- 1000
#!      n <- 100000
#!      k <- 3
#!      memory.size(max = TRUE)
#!      system.time(
#!      x <- ff(vmode="double", dim=c(b*n,k), dimnames=list(NULL, LETTERS[1:k]))
#!      )
#!      memory.size(max = TRUE)
#!      system.time(
#!      ffrowapply({
#!         l <- i2 - i1 + 1
#!         z <- rnorm(l)
#!         x[i1:i2,] <- z + matrix(rnorm(l*k), l, k)
#!      }, X=x, VERBOSE=TRUE, BATCHSIZE=n)
#!      )
#!      memory.size(max = TRUE)
#!
#!      form <- A ~ B + C
#!      first <- TRUE
#!      system.time(
#!      ffrowapply({
#!         if (first){
#!           first <- FALSE
#!           fit <- biglm(form, as.data.frame(x[i1:i2,,drop=FALSE]))
#!         }else
#!           fit <- update(fit, as.data.frame(x[i1:i2,,drop=FALSE]))
#!      }, X=x, VERBOSE=TRUE, BATCHSIZE=n)
#!      )
#!      memory.size(max = TRUE)
#!      first
#!      fit
#!      summary(fit)
#!      rm(x); gc()
#!   }
#! }
#! \keyword{ IO }
#! \keyword{ array }
#! \keyword{ attribute }
#! \keyword{ classes }
#! \keyword{ package }


ff <- function(
  initdata    = NULL
, length      = NULL
, levels      = NULL
, ordered     = NULL
, dim         = NULL
, dimorder    = NULL    # this is the (transparent) storage layout
, bydim       = NULL    # this is the dimorder used to read in the initdata (e.g. use 2:1 to mimic matrix(,byrow=TRUE)) (passed to update)
, symmetric   = FALSE   # allows to define a matrix as symmetric: in this case we need all(diff(dim)==0) and we assume that all subscript-combinations are equal (storage mapped to sorted subscripts, e.g. distance matrices), no dimorder allowed
, fixdiag     = NULL    # symmetric matrix only: value for fixdiag if diagonal is redundant (e.g. fixdiag=0 for dist)
, names       = NULL    # not taken fom initdata
, dimnames    = NULL    # not taken fom initdata
, ramclass    = NULL
, ramattribs  = NULL
, vmode       = NULL    # be default we get the vmode from initdata
, update      = NULL    # set to FALSE to suppress upating ff object with initdata
, pattern     = NULL
, filename    = NULL
, overwrite   = FALSE
, readonly    = FALSE
, pagesize    = NULL    # getOption("ffpagesize")
, caching     = NULL    # getOption("ffcaching")
, finalizer   = NULL    # "delete" for tempfiles and "close" for named files
, finonexit   = NULL    # getOption("fffinonexit")
#, hideclass = FALSE   # experimental, please ignore (would be needed to bypass <-.ff)
, FF_RETURN   = TRUE
, BATCHSIZE   = .Machine$integer.max       # optional batch size restriction in cases (limiting is useful if pulling values from function accesses raw data with more columns that k)
, BATCHBYTES  = getOption("ffbatchbytes")  # batch size restriction in bytes (needed for update() if length(initdata)>1)
, VERBOSE     = FALSE
)
{
  if (!is.logical(FF_RETURN) || length(FF_RETURN)!=1)
      stop("in ff() FF_RETURN must be a logical scalar")

  if (is.null(symmetric)||!is.logical(symmetric)||length(symmetric)!=1)
    stop("symmetric must be FALSE or TRUE")

  if (is.null(vmode)){
    if (is.null(initdata))
      stop("need vmode or initdata")
    vmode <- vmode(initdata)
  } else{
    if (is.na(match(vmode, c("boolean", "logical", "quad", "nibble", "byte", "ubyte", "short", "ushort", "integer", "single", "double", "raw")))){
     stop("vmode '", vmode,"' not implemented")
    }
  }

  if (!ffxtensions()){
    if (!ffsymmxtensions()){
      symm <- function(...).NotYetImplemented()
      if (symmetric){
        message("You are requesting a dual-licence feature that currently is only available to parties who support the development of package ff and friends")
       .NotYetUsed("symmetric", error = TRUE)
      }
      if (!is.null(fixdiag)){
        message("You are requesting a dual-licence feature that currently is only available to parties who support the development of package ff and friends")
       .NotYetUsed("fixdiag", error = TRUE)
      }
    }
    if (!is.na(match(vmode, c("boolean", "quad", "nibble", "byte", "ubyte", "short", "ushort", "single")))){
      message("You are requesting a dual-licence feature that currently is only available to parties who support the development of package ff and friends")
      stop("vmode='", vmode, "' not on CRAN")
    }
  }

  # determine filename and finalizer
  if (is.null(filename)){
    if (is.null(pattern))
      pattern <- "ff"
    filename <- fftempfile(pattern)
  }else{
    if (is.null(pattern))
      pattern <- paste(splitPathFile(filename)$path, "/", sep="")
  }
  # gurantee absolute path
  dfile <- dirname(filename)
  bfile <- basename(filename)
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(dfile)
  dfile <- getwd()
  filename <- file.path(dfile, bfile)
  # fix problem in file.path
  filename <- gsub("/+","/",filename)

  if (is.null(finalizer)){
    finalizer <- getOption("fffinalizer")
    if (is.null(finalizer)){
      if (dfile==getOption("fftempdir"))
        finalizer <- "delete"   # temporary ff object
      else
        finalizer <- "close"    # persistent ff object
    }
  }else{
    finalizer <- match.arg(finalizer, choices=c("deleteIfOpen", "delete", "close"))
  }

  # handle file reuse and set maxlength (0=ram object or reuse of ff, otherwise = size of new ff file in bytes)
  if (FF_RETURN){
    if ( file.exists(filename) ){
      if (overwrite){
        maxlength <- 0L
      }else{
        if (file.access(filename,4) == -1){ # no read access
          stop("read permission denied for file")
        }
        if (file.access(filename,2) == -1){ # no write access
          if (!readonly) {
            readonly <- TRUE
            warning("force read-only access on file")
          }
        }
        filesize <- file.info(filename)$size
        if (is.na(filesize))
          stop("unable to open file")
        if (!is.null(initdata)) {
          stop("bad argument initdata for existing file; initializing existing file is invalid")
        }
        fillength <- filesize/.ffbytes[vmode]
        if (fillength>.Machine$integer.max){
          warning("file contains more than .Machine$integer.max elements of this vmode")
          maxlength <- as.integer(floor(.Machine$integer.max * .ffbytes[vmode] / .rambytes[vmode]) * .rambytes[vmode])
        }else
          maxlength <- as.integer(fillength)
      }
    }else{
      if (readonly)
        stop("allocation of a new 'readonly' flat file vector not supported")
      maxlength <- 0L
    }
  }else{
    maxlength <- 0L
  }

  # handle levels
  if (is.null(levels))
    levels <- levels(initdata)
  if (!is.null(levels)){
    if (vmode=="character")
      vmode <- "integer"
    if (.rammode[vmode] != "integer")
        stop("factors must be .rammode integer")
    if (length(levels)>.vmax[vmode]+.vunsigned[vmode])
      stop("vmode '", vmode, "' can carry max ", .vmax[vmode]+.vunsigned[vmode], " levels")
    if (is.null(ramclass)){
      if (is.null(ordered))
        ordered <- is.ordered(initdata)
      if (ordered)
        ramclass <- c("ordered","factor")
      else
        ramclass <- "factor"
    }
  }
  # handle initdata
  if (!is.null(initdata)){
    if (is.null(ramclass))
      ramclass <- ramclass(initdata)
    if (is.null(ramattribs)){
      ramattribs <- ramattribs(initdata)
    }
  }
  if (!is.atomic(initdata[1]))
    stop("initdata[1] must be atomic")

  # handle or derive length from dim or initdata
  if (is.null(dim)){ # ff_vector
    if (is.null(length)){
      if (maxlength){
        length <- maxlength
      }else{
        if (is.null(initdata))
          stop("need length or initdata")
        else{
          length <- length(initdata)
        }
      }
    }
    if(!is.null(dimorder))
      stop("dimorder must be null with vectors")
    if (symmetric)
      stop("symmetric only allowed with matrices")
    if (!is.null(fixdiag))
      stop("fixed-diagonal only allowed with symmetric matrices")
    length <- as.integer(length)
    ffclass <- c("ff_vector","ff")
  }else{ # ff_array || ff_symm
    dim <- as.integer(dim)
    ndim <- length(dim)
    if (is.null(dimorder))
      dimorder <- seq(length.out=ndim)
    else{
      dimorder <- as.integer(dimorder)
      if (!identical(sort(dimorder), 1:ndim))
        stop("dimorder does not match dimension")
    }
    if (symmetric){ # ff_symm
      if (ndim!=2 || dim[1]!=dim[2])
        stop("symmetric matrices require parameter dim with 2 equal values")
      if (!dimorderStandard(dimorder))
        stop("non-standard dimorder not allowed for symmetric matrices")
      if (!is.null(bydim))
        stop("bydim not allowed for symmetric matrices")
      if (is.null(fixdiag))
        length <- dim[1]*(dim[1]-1L)/2L + dim[1]
      else
        length <- dim[1]*(dim[1]-1L)/2L
      ffclass <- c("ff_symm", "ff")
    }else{
      if(!is.null(fixdiag))
        stop("fixed-diagonal only allowed with symmetric matrices")
      if (!is.null(bydim)){
        bydim <- as.integer(bydim)
        if (!identical(sort(bydim), 1:ndim))
          stop("bydim does not match dimension")
      }
      n <- as.integer(prod(dim))
      if(!is.null(length) && length!=n)
        stop("dim and length don't match")
      else
        length <- n
      if (ndim==2)
        ffclass <- c("ff_matrix", "ff_array","ff")
      else
        ffclass <- c("ff_array","ff")
    }
  }
  if (length<0 || length>.Machine$integer.max)
    stop("length must be between 1 and .Machine$integer.max")

  if (is.null(pagesize))
    pagesize <- getOption("ffpagesize")
  else
    pagesize <- getalignedpagesize(pagesize)
  if (is.null(caching))
    caching <- getOption("ffcaching")
  else
    caching <- match.arg(caching, caching_schemes)
  if (is.null(finonexit))
    finonexit <- getOption("fffinonexit")

	# rr <- logical(1)
	# rr[1] <- readonly[1]
		
  pattr <- list(  # physical attributes
    vmode     = vmode
  , maxlength = if (maxlength){if (length > maxlength) stop("length exceeds file length") else maxlength} else length # the physical file size
  , pattern   = pattern
  , filename  = filename
  , pagesize  = pagesize
  , finalizer = finalizer
  , finonexit = finonexit
  , readonly  = readonly 
  , caching   = caching
  , class     = "ff_pointer"  # class of pointer within ff class needed for finalizer dispatch
  )

  vattr <- list(
    Length      = length      # the current length we use from R
  , Dim         = dim
  , Dimorder    = dimorder
  , Symmetric   = symmetric
  , Fixdiag     = fixdiag
  , Levels      = levels
  , ramclass    = ramclass
  , ramattribs  = ramattribs
  )

  if (is.logical(FF_RETURN) && length(FF_RETURN)==1){
    if (!FF_RETURN){
      if (length*.rambytes[vmode] > getOption("ffbatchbytes"))
        warning("creating large ram object with ", length*.rambytes[vmode], " > ", getOption("ffbatchbytes"))
      if (is.null(initdata)){
        initdata <- vector.vmode(vmode,1)
      }else{
        if (is.null(levels))
          initdata <- as.vmode(initdata[], vmode) # NOTE that [] takes care of dimorder, [i1:i2] would not
        else
          initdata <- as.vmode(ram2ramcode(initdata[], levels), vmode) # NOTE that [] takes care of dimorder, [i1:i2] would not
      }
      ret <- switch(ffclass[1]
      , ff_symm   = { if (!is.null(levels)) stop("symm currently not defined for factors"); symm(initdata, dim=dim, dimnames = dimnames, fixdiag = fixdiag)}
      , ff_vector = {temp <- rep(initdata, length.out=length); names(temp)<-names; temp}
      , ff_matrix = {temp <- vector2array(initdata, dim=dim, dimorder=bydim); dimnames(temp)<- dimnames;temp}
      , ff_array  = {temp <- vector2array(initdata, dim=dim, dimorder=bydim); dimnames(temp)<- dimnames;temp}
      )
      if (!is.null(levels)){
        attr(ret, "levels") <- levels
      }
      if (!is.null(ramattribs))
        attributes(ret) <- c(attributes(ret), ramattribs)
      physical(ret) <- pattr[!is.na(match(names(pattr), ramphysical_includes))]
       virtual(ret) <- vattr[!is.na(match(names(vattr), ramvirtual_includes ))]
      attr(ret, "vmode") <- vmode
      if (!is.null(ramclass))
        class(ret) <- ramclass
      return(ret)
    }
  }else{
    stop("in ff() FF_RETURN must be a logical scalar")
  }

  # stopifnot( file.access(filename,0) || (!file.access(filename,2) && overwrite ) )  # file must not exist OR be writable and overwrite
  # create fast file and return external pointer
  initval <- if (is.null(levels)) {
    if (length(initdata))
      as.vmode(initdata[1], vmode)
    else
      vector.vmode(vmode, 1)
  } else {
    if (length(initdata)){
      if (.vunsigned[vmode])
        as.vmode(match(initdata[1], levels), vmode) - 1L
      else
        as.vmode(match(initdata[1], levels), vmode)
    }else{
      if (.vunsigned[vmode])
        0L
      else
        1L
    }
  }
  ffpointer <- .Call("new"
  , as.character(filename)
  , .ffmode[vmode]
  , initval
  , if (maxlength) 0L else length # already integer
  , pagesize    # already integer
  , readonly
  , caching == "mmeachflush"
  , PACKAGE="ff")

  #if (hideclass){
  #  a[[".class"]] <- a[["class"]]
  #  a[["class"]] <- NULL
  #}
  attributes(ffpointer) <- pattr  # unusual copying semantics: changing an attribute of a 'copy' of an ff object changes the attributes of ALL copies
  reg.finalizer(ffpointer, finalize.ff_pointer, onexit=finonexit)  # for details see ?finalize
  v <- list()
  attributes(v) <- vattr
  ret <- list()
  attributes(ret) <- list(physical=ffpointer, virtual=v, class=ffclass)
  # now this is a legal ff object

  # complete initialization
  if ( length(initdata) && (is.null(update) || update) ){
    if ( (!identical(as.vector(initdata[1]), as.vector(ret[1])) || length(initdata)>1) )
      ret <- update.ff(ret
      , from    = initdata
      , delete  = FALSE
      , bydim   = bydim
      , BATCHSIZE   = BATCHSIZE
      , BATCHBYTES  = BATCHBYTES
      , VERBOSE     = VERBOSE
      )
  }else{
    # note that at this point new ff objects are initialized with 0
    # thus unsigned factor are - without writing - initialized at their first level
    # signed factor need to be initialized with NA in order not to be in an illegal state
    if (!is.null(levels) && !.vunsigned[vmode])
      ret[] <- NA
  }
  if (is.null(dim)){
    if (!is.null(names)) # !hideclass &&
      names(ret) <- names
  }else{
    if (!is.null(dimnames)) # !hideclass &&
      dimnames(ret) <- dimnames
  }
  return(ret)
}



#! \name{update.ff}
#! \alias{update.ff}
#! \alias{update.ffdf}
#! \title{ Update ff content from another object }
#! \description{
#!   \command{update} copies updates one ff object with the content of another object.
#! }
#! \usage{
#! \method{update}{ff}(object, from, delete = FALSE, bydim = NULL, fromdim = NULL
#! , BATCHSIZE = .Machine$integer.max, BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE, \dots)
#! \method{update}{ffdf}(object, from, \dots)
#! }
#! \arguments{
#!   \item{object}{ an ff object to which to update }
#!   \item{from}{ an object from which to uodate  }
#!   \item{delete}{ NA for quick update with file-exchange, TRUE for quick update with deleting the 'from' object after the update, can speed up updating significantly }
#!   \item{bydim}{ how to interpret the content of the object, see \code{\link{ff}} }
#!   \item{fromdim}{ how to interpret the content of the 'from' object, see \code{\link{ff}} }
#!   \item{BATCHSIZE}{ \code{BATCHSIZE} }
#!   \item{BATCHBYTES}{ \code{BATCHBYTES} }
#!   \item{VERBOSE}{ \code{VERBOSE} }
#!   \item{\dots}{ further arguments }
#! }
#! \details{
#!   If the source object \code{is.ff} and not \code{delete=FALSE} then instead of slow copying we - if possible - try to swap and rename the files behind the ff objects.
#!   Quick update requires that the two ff objects are \code{\link{vectorCompatible}},
#!   that both don't use \code{\link{vw}},
#!   that they have identical \code{\link{maxlength}}
#!   and identical \code{\link{levels.ff}}.
#! }
#! \note{
#!   You don't have a guarantee that with \code{delete=TRUE} the 'from' object gets deleted or with \code{delete=NA} the 'from' objects carries the content of 'object'.
#!   Such expectations only turn true if really a quick update was possible.
#! }
#! \value{
#!   An ff object like the input 'object' updated with the content of the 'from' object.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ff}}, \code{\link{clone}}, \code{\link{ffvecapply}}, \code{\link{vectorCompatible}}, \code{\link{filename}} }
#! \examples{
#!   x <- ff(1:100)
#!   y <- ff(-(1:100))
#!   message("You should make it a habit to re-assign the return value 
#! of update although this is not needed currently.")
#!   x <- update(x, from=y)
#!   x
#!   y
#!   x[] <- 1:100
#!   x <- update(x, from=y, delete=NA)
#!   x
#!   y
#!   x <- update(x, from=y, delete=TRUE)
#!   x
#!   y
#!   x
#!   rm(x,y); gc()
#!
#!   \dontrun{
#!     message("timings")
#!     x <- ff(1:10000000)
#!     y <- ff(-(1:10000000))
#!     system.time(update(x, from=y))
#!     system.time(update(y, from=x, delete=NA))
#!     system.time(update(x, from=y, delete=TRUE))
#!     rm(x,y); gc()
#!   }
#!
#! }
#! \keyword{ IO }
#! \keyword{ data }


# update 'x' with content of 'from' using recycling
# if delete=TRUE and from is.ff of same maxlength then we do fast update by renaming behind the scences
# the latter is especially useful when 'transferring' data from temporary ff to target ff, see sort.ff



update.ff <- function(
  object
, from
, delete        = FALSE
, bydim         = NULL
, fromdim       = NULL
, BATCHSIZE     = .Machine$integer.max       # optional batch size restriction in cases (limiting is useful if pulling values from function accesses raw data with more columns that k)
, BATCHBYTES    = getOption("ffbatchbytes")  # batch size restriction in bytes
, VERBOSE       = FALSE
, ...
)
{
  dto <- dim(object)
  dfrom <- dim(from)
  ndto <- length(dto)
  ndfrom <- length(dfrom)
  if (is.null(dfrom) && !is.null(fromdim))
    stop("from has no dimension for fromdim")
  if (is.null(dto) && !is.null(bydim))
    stop("from has no dimension for fromdim")
  if (inherits(from, "ff")){
    doto <- dimorder(object)
    dofrom <- dimorder(from)
    #dotodev <- !is.null(doto) && !identical(doto, sort(doto))
    #dofromdev <- !is.null(dofrom) && !identical(dofrom, sort(dofrom))
    #dodev <- (dotodev || dofromdev) && !identical(doto, dofrom)
    #dovector <- is.null(bydim) && is.null(fromdim) && !dodev
    if ( vmode(object)==vmode(from)
    && vectorCompatible(dim=dto, dim2=dfrom, dimorder=doto, dimorder2=dofrom, bydim=bydim, bydim2=fromdim)
    && is.null(vw(from)) && is.null(vw(object))
    && (is.na(delete) || delete)
    && length(object)==length(from)
    && maxlength(object)==maxlength(from)
    && identical(levels(from), levels(object))
    ){
      # fast update by file rename
      if (is.open(object)){
        reopento <- TRUE
        close(object)
      }else
        reopento <- FALSE
      if (is.open(from)){
        if (is.na(delete))
          reopenfrom <- TRUE
        else
          reopenfrom <- FALSE
        close(from)
      }else
        reopenfrom <- FALSE

      if (is.na(delete)){
        # do exchange
        tmpfilename <- fftempfile("update")
        objfilename <- filename(object)
        fromfilename <- filename(from)
        if(!file.move(objfilename, tmpfilename))
          stop("renaming object file '", objfilename, "' to temp file '", tmpfilename, "' failed in update(...,delete=NA)")
        if(!file.move(fromfilename, objfilename))
          stop("renaming from file '", fromfilename, "' to object file '", objfilename, "' failed in update(...,delete=NA)")
        if(!file.move(tmpfilename, fromfilename))
          stop("renaming temp file '", tmpfilename, "' to from file '", fromfilename, "' failed in update(...,delete=NA)")
      }else{
        # do plug in and delete
        oldnam <- filename(object)
        if(!file.remove(oldnam))
          stop("removing from file '", filename(from), "' failed in update(..., delete=TRUE)")
        if(!file.move(filename(from), oldnam))
          stop("renaming from file '", filename(from), "' to '", oldnam, "' failed in update(...,delete=TRUE)")
      }

      if (reopento)
        open(object)
      if (reopenfrom)
        open(from)
    }else{
      if (is.open(object)){
        recloseto <- FALSE
      }else{
        recloseto <- TRUE
        open(object)
      }
      if (is.open(from)){
        reclosefrom <- FALSE
      }else{
        reclosefrom <- TRUE
        open(from)
      }
      if (ndto && ndto==ndfrom && all(dto==dfrom)){
        ndim <- ndto
        args <- rep(alist(a = ), ndim)
        argsfrom <- c(args, alist(bydim=fromdim))
        argsto <- c(args, alist(bydim=bydim))
        ret <- ffapply(
          , EXPR = {
            argsfrom[1:ndim] <- argsto[1:ndim] <- lapply(1:ndim, function(i)substitute(b1:b2, list(b1=i1[i], b2=i2[i])))  # substitute avoids unpacking the sequences !!
            temp <- do.call("[", c(list(from), argsfrom, list(drop = FALSE)))
            object <- do.call("[<-", c(list(object), argsto, list(value=temp)))
          }
        , X = object
        , MARGIN  = 1:ndim
        , BATCHSIZE   = BATCHSIZE
        , BATCHBYTES  = BATCHBYTES
        , VERBOSE     = VERBOSE
        )
      }else{
        nto <- length(object)
        nfrom <- length(from)
        if (nto%%nfrom)
          warning("length(object) not a multiple of length(from) in update.ff(object, from, ...)")

        i1 <- i2 <- 0L  # dummy assignment to shut up R CMD CHECK NOTE about no visible binding
        if (nfrom<nto){
          ffvecapply(
            EXPR = object[i1:i2] <- repfromto(from, i1, i2)
          , X = object
          , BATCHSIZE   = BATCHSIZE
          , BATCHBYTES  = BATCHBYTES
          , VERBOSE     = VERBOSE
          )
        }else{
          ffvecapply(
            EXPR = object[i1:i2] <- from[i1:i2]
          , X = object
          , BATCHSIZE   = BATCHSIZE
          , BATCHBYTES  = BATCHBYTES
          , VERBOSE     = VERBOSE
          )
        }
      }

      if (reclosefrom)
        close(from)
      if (recloseto)
        close(object)
    }
  }else{
    # if is.ram(from) we assign all in one chunk and let [<-.ff recycle
    # dimorder(from) is always standard
    # dimorder(to) and bydim is handled by [<-.ff
    # xx we only need to care about fromdim
    if (!is.null(fromdim)){
      fromLevels <- levels(from)
      if (is.null(fromLevels)){
        from <- array2vector(from, dim=dfrom, dimorder=fromdim)
      }else{
        fromramclass <- ramclass(from)
        from <- array2vector(from, dim=dfrom, dimorder=fromdim)
        attr(from, "levels") <- fromLevels
        class(from) <- fromramclass
      }
    }
    if (is.open(object)){
      recloseto <- FALSE
    }else{
      recloseto <- TRUE
      open(object)
    }
    if (is.null(bydim))
      object[] <- from
    else
      object[,bydim=bydim] <- from
    if (recloseto)
      close(object)
  } # end is.ram(from)
  object
}


#! \name{clone}
#! \alias{clone}
#! \alias{clone.ff}
#! \alias{clone.list}
#! \alias{clone.default}
#! \title{ Cloning ff and ram objects }
#! \description{
#!   \command{clone} physically duplicates ff (and ram) objects and can additionally change some features, e.g. length.
#! }
#! \usage{
#! clone(x, \dots)
#! \method{clone}{ff}(x
#! , initdata = x
#! , length = NULL
#! , levels = NULL
#! , ordered = NULL
#! , dim = NULL
#! , dimorder = NULL
#! , bydim = NULL
#! , symmetric = NULL
#! , fixdiag = NULL
#! , names = NULL
#! , dimnames = NULL
#! , ramclass = NULL
#! , ramattribs = NULL
#! , vmode = NULL
#! , update  = NULL
#! , pattern = NULL
#! , filename = NULL
#! , overwrite = FALSE
#! , pagesize = NULL
#! , caching = NULL
#! , finalizer = NULL
#! , finonexit = NULL
#! , FF_RETURN = NULL
#! , BATCHSIZE = .Machine$integer.max
#! , BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE
#! , \dots)
#! \method{clone}{list}(x, \dots)
#! \method{clone}{default}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ \code{x} }
#!   \item{initdata}{ scalar or vector of the \code{\link{.vimplemented}} \code{\link{vmode}}s, recycled if needed, default 0, see also \code{\link{as.vmode}} and \code{\link{vector.vmode}} }
#!   \item{length}{ optional vector \code{\link{length}} of the object (default: derive from 'initdata' or 'dim'), see \code{\link{length.ff}} }
#!   \item{levels}{ optional character vector of levels if (in this case initdata must be composed of these) (default: derive from initdata) }
#!   \item{ordered}{ indicate whether the levels are ordered (TRUE) or non-ordered factor (FALSE, default) }
#!   \item{dim}{ optional array \code{\link{dim}}, see \code{\link{dim.ff}} and \code{\link{array}} }
#!   \item{dimorder}{ physical layout (default 1:length(dim)), see \code{\link{dimorder}} and \code{\link{aperm}} }
#!   \item{bydim}{ dimorder by which to interpret the 'initdata', generalization of the 'byrow' paramter in \code{\link{matrix}} }
#!   \item{symmetric}{ extended feature: TRUE creates symmetric matrix (default FALSE) %, see \code{\link{symm}}, \code{\link{ff_symm}}, \code{\link{ff_dist}}
#!        }
#!   \item{fixdiag}{ extended feature: non-NULL scalar requires fixed diagonal for symmetric matrix (default NULL is free diagonal) }
#!   \item{names}{ see \code{\link{names}} }
#!   \item{dimnames}{ NOT taken from initdata, see \code{\link{dimnames}} }
#!   \item{ramclass}{ class attribute attached when moving all or parts of this ff into ram, see \code{\link{ramclass}} }
#!   \item{ramattribs}{ additional attributes attached when moving all or parts of this ff into ram, see \code{\link{ramattribs}} }
#!   \item{vmode}{ virtual storage mode (default: derive from 'initdata'), see \code{\link{vmode}} and \code{\link{as.vmode}} }
#!   \item{update}{ set to FALSE to avoid updating with 'initdata' (default TRUE) (used by \code{\link{ffdf}}) }
#!   \item{pattern}{ root pattern for automatic ff filename creation (default "ff"), see also \code{\link[=physical.ff]{physical}} }
#!   \item{filename}{ ff \code{\link{filename}} (default tmpfile with 'pattern' prefix), see also \code{\link[=physical.ff]{physical}} }
#!   \item{overwrite}{ set to TRUE to allow overwriting existing files (default FALSE) }
#!   \item{pagesize}{ pagesize in bytes for the memory mapping (default from getOptions("ffpagesize") initialized by \code{\link{getdefaultpagesize}}), see also \code{\link[=physical.ff]{physical}} }
#!   \item{caching}{ caching scheme for the backend, currently 'mmnoflush' or 'mmeachflush' (flush mmpages at each swap, default from getOptions("ffcaching") initialized with 'memorymap'), see also \code{\link[=physical.ff]{physical}} }
#!   \item{finalizer}{ name of finalizer function called when ff object is \code{\link{remove}d}, (default "deleteIfOpen" from getOptions("fffinalizer"))), standard finalizers are \code{\link{close.ff}}, \code{\link{delete.ff}} and \code{\link{deleteIfOpen.ff}}, see also \code{\link[base]{reg.finalizer}} }
#!   \item{finonexit}{ logical scalar determining whether finalizer is also called when R is closed via \code{\link{q}}, (default TRUE from getOptions("fffinonexit")) }
#!   \item{FF_RETURN}{ logical scalar or ff object to be used. The default NULL creates a ff or ram clone, TRUE returns a ff clone, FALSE returns a ram clone. Handing over an ff object here uses this or stops if not \code{\link{ffsuitable}} }
#!   \item{BATCHSIZE}{ integer scalar limiting the number of elements to be processed in \code{\link{update.ff}} when length(initdata)>1, default from getOption("ffbatchsize") }
#!   \item{BATCHBYTES}{ integer scalar limiting the number of bytes to be processed in \code{\link{update.ff}} when length(initdata)>1, default from getOption("ffbatchbytes"), see also \code{\link{.rambytes}} }
#!   \item{VERBOSE}{ set to TRUE for verbosing in \code{\link{update.ff}} when length(initdata)>1, default FALSE }
#!   \item{\dots}{ further arguments to the generic }
#! }
#! \details{
#!   \command{clone} is generic. \command{clone.ff} is the workhorse behind \code{\link{as.ram}} and \code{\link{as.ff}}.
#!   For creating the desired object it calls \code{\link{ff}} which calls \code{\link{update}} for initialization.
#! }
#! \value{
#!   an ff or ram object
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ff}}, \code{\link{update}}, \code{\link{as.ram}}, \code{\link{as.ff}} }
#! \examples{
#!   x <- ff(letters, levels=letters)
#!   y <- clone(x, length=52)
#!   rm(x,y); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }


# clone can change the length of the cloned object, but not create subsetted versions
# increased length is filled by recycling
clone.ff <- function(
  x
, initdata    = x
, length      = NULL
, levels      = NULL
, ordered     = NULL
, dim         = NULL
, dimorder    = NULL
, bydim       = NULL
, symmetric   = NULL
, fixdiag     = NULL
, names       = NULL    # not taken fom initdata
, dimnames    = NULL
, ramclass    = NULL
, ramattribs  = NULL
, vmode       = NULL
, update      = NULL    # set to TRUE to suppress upating ff object with initdata
, pattern     = NULL
, filename    = NULL
, overwrite   = FALSE
, pagesize    = NULL
, caching     = NULL
, finalizer   = NULL
, finonexit   = NULL
#, hideclass  = FALSE
, FF_RETURN   = NULL
, BATCHSIZE   = .Machine$integer.max       # optional batch size restriction in cases (limiting is useful if pulling values from function accesses raw data with more columns that k)
, BATCHBYTES  = getOption("ffbatchbytes")  # batch size restriction in bytes
, VERBOSE     = FALSE
, ... # dummy to keep R CMD check quiet
)
{
  # BTW: we tried fast cloning via file.copy() but that was slower, only system(copy...) would save 33% time (too wacky)
  if (is.null(vmode) && !is.null(x))
    vmode <- vmode(x)

  if (is.null(levels))
    levels <- levels(x)
  if (is.null(ordered))
    ordered <- is.ordered(x)
  if (is.null(names))
    names <- names(x)
  if (is.null(ramclass))
    ramclass <- ramclass(x)
  if (is.null(ramattribs))
    ramattribs <- ramattribs(x)

  oldlen <- length(x)
  if (is.null(dim) && !is.null(length)){  # not dim but length given
    dimorder <- NULL
    symmetric <- FALSE
    fixdiag <- NULL
    dimnames <- NULL
  }else{ # dim given or none of dim/length
    if (is.null(dim))
      dim <- dim(x)
    if (is.null(dim)){
      if (is.null(length))
        length <- length(x)
    }else length <- as.integer(prod(dim))
    if (is.null(dimorder)){
      dimorder <- virtual(x)$Dimorder  # NOT: dimorder(x) because this would not restore the dimorder from ram objects
    }
    if (is.null(symmetric))
      symmetric <- symmetric(x)
    if (is.null(fixdiag))
      fixdiag <- fixdiag(x)
    if (identical(dim, dim(x))){
      if (is.null(dimnames))
        dimnames <- dimnames(x)
    }else
      dimnames <- NULL
  }

  if (is.null(FF_RETURN))
    FF_RETURN <- is.ff(x)
  if ( !(is.logical(FF_RETURN) && length(FF_RETURN)==1) )
    stop("in clone() FF_RETURN must be a logical scalar")

  physical <- physical(x)
  if (is.null(pattern)){
    pattern <- physical$pattern
    if (is.null(pattern))
      pattern <- "clone"
  }
  if (is.null(filename) && !is.ff(x))
    filename <- physical$filename
  if (is.null(pagesize))
    pagesize <- physical$pagesize
  if (is.null(caching))
    caching <- physical$caching
  if (is.null(finonexit))
    finonexit <- physical$finonexit
  if (is.null(finonexit))
    finonexit <- physical$finonexit

  # don't use "<-" operator with ff argument in order to avoid recursion (if anyone defines method <-.ff as cloning)
  assign("ret", ff(
    initdata    = initdata
  , length      = length
  , levels      = levels
  , ordered     = ordered
  , dim         = dim
  , dimorder    = dimorder
  , bydim       = bydim
  , symmetric   = symmetric
  , fixdiag     = fixdiag
  , names       = names
  , dimnames    = dimnames
  , ramclass    = ramclass
  , ramattribs  = ramattribs
  , vmode       = vmode
  , update      = update
  , pattern     = pattern
  , filename    = filename
  , readonly    = FALSE
  , overwrite   = overwrite
  , pagesize    = pagesize
  , caching     = caching
  , finalizer   = finalizer
  , finonexit   = finonexit
  #, hideclass  = hideclass
  , FF_RETURN   = FF_RETURN
  , BATCHSIZE   = BATCHSIZE
  , BATCHBYTES  = BATCHBYTES
  , VERBOSE     = VERBOSE
  ))
  newlen <- length(ret)
  nam <- names(x)
  if (is.null(dim) && !is.null(nam)){
    if (newlen==oldlen)
      names(ret) <- names(x)
    else if (newlen>oldlen)
      names(ret) <- c(names(x), (oldlen+1L):newlen)
    else
      names(ret) <- names(x)[1:newlen]
  }
  if (!is.null(physical$na.count)){
    if (newlen==oldlen)
      na.count(x) <- physical$na.count
    else{
      na.count(x) <- NULL
      warning("cloning removed attribute 'na.count'")
   }
  }
  if (!is.null(physical$is.sorted)){
    if (newlen<=oldlen && is.null(dim(x)) && is.null(dim(ret)))
      is.sorted(ret) <- physical$is.sorted
    else{
      is.sorted(ret) <- FALSE
      warning("cloning set 'is.sorted' to FALSE")
    }
  }
  ret
}

clone.default <- function(x
, ... # passed to clone
){
  if (is.atomic(x)){
    if (length(x))
      x[1] <- x[1]  # force a copy around COPY ON MODIFY
    x
  }else{
    stop("clone not defined for type")
  }
}

clone.list <- function(x
, ... # passed to clone
){
  lapply(x, clone, ...)
}


#! \name{finalizer}
#! \Rdversion{1.1}
#! \alias{finalizer}
#! \alias{finalizer<-}
#! \alias{finalizer.ff}
#! \alias{finalizer<-.ff}
#! \title{
#!   Get and set finalizer (name)
#! }
#! \description{
#!   The generic \code{finalizer} allows to get the current finalizer. The generic \code{finalizer<-} allows to set the current finalizer or to change an existing finalizer (but not to remove a finalizer).
#! }
#! \usage{
#! finalizer(x, ...)
#! finalizer(x, ...) <- value
#! \method{finalizer}{ff}(x, ...)
#! \method{finalizer}{ff}(x, ...) <- value
#! }
#! \arguments{
#!   \item{x}{an \code{\link{ff}} object}
#!   \item{value}{the name of the new finalizer}
#!   \item{\dots}{ignored}
#! }
#! \details{
#!   If an \code{\link{ff}}  object is created a finalizer is assigned, it has the task to free ressources no longer needed, for example remove the ff file or free the C++ RAM associated with an open ff file.
#!   The assigned finalizer depends on the location of the ff file:
#!   if the file is created in \code{getOption(fftempdir)} it is considered considered temporary and has default finalizer \code{\link[ff:delete.ff]{delete}},
#!   files created in other locations have default finalizer \code{\link[ff:close.ff]{close}}.
#!   The user can override this either by setting \code{options("fffinalizer")} or by using argument \code{finalizer} when creating single \code{ff} objects.
#!   Available finalizer generics are "close", "delete" and "deleteIfOpen", available methods are \code{\link{close.ff}}, \code{\link{delete.ff}} and \code{\link{deleteIfOpen.ff}}.
#!   \cr
#!   In order to be able to change the finalizer before finalization, the finalizer is NOT directly passed to R's finalization mechanism \code{\link[base]{reg.finalizer}} (an active finalizer can never be changed other than be executed).
#!   Instead the NAME of the desired finalizer is stored in the ff object and \code{\link{finalize.ff_pointer}} is passed to \code{reg.finalizer}.
#!   \code{finalize.ff_pointer} will at finalization-time determine the desired finalizer and call it.
#!   \cr
#!   There are two possible triggers for execution \code{finalize.ff_pointer}:
#!   \enumerate{
#!     \item the garbage collection \code{\link{gc}} following removal \code{\link{rm}} of the ff object
#!     \item closing R if \code{finonexit} was \code{TRUE} at ff creation-time, determined by \code{options("fffinonexit")} and ff argument \code{finonexit}
#!   }
#!   Furthermore there are two possible triggers for calling the finalizer
#!   \enumerate{
#!     \item an explicit call to \code{\link{finalize}}
#!     \item an explicit call to one of the finalizers \code{\link[ff:close.ff]{close}}, \code{\link{delete}} and \code{\link{deleteIfOpen}}
#!   }
#!   The user can define custom finalizers by creating a generic function like \code{\link{delete}}, a ff_pointer method like \code{\link{delete.ff_pointer}} and a ff method for manual calls like \code{\link{delete.ff}}. The user then is responsible to take care of two things
#!   \enumerate{
#!     \item adequate freeing of ressources
#!     \item proper maintenance of the finalizer name in the ff object via \code{\link[=physical.ff]{physical}$finalizer}
#!   }
#!   \code{is.null(finalizer(ff))} indicates NO active finalizer, i.e. no pending execution of \code{finalize.ff_pointer} lurking around after call of \code{reg.finalizer}.
#!   This requires that
#!   \enumerate{
#!     \item the \code{ff_pointer} method sets the finalizer name to \code{NULL}
#!     \item the \code{ff} may change a non-NULL finalizer name to a different name but not change it to NULL
#!   }
#! }
#! \value{
#!   \code{finalizer} returns the name of the active finalizer or \code{NULL} if no finalizer is active. \cr
#!   \code{finalizer<-} returns the changed ff object (reassignment of this return value not needed to keep the change).
#!   If there was no pending call to \code{\link{finalize.ff_pointer}} (\code{is.null(finalizer(ff))}), \code{finalizer<-} will create one by calling \code{reg.finalizer} with the current setting of \code{\link[=physical.ff]{physical}$finonexit}.
#! }
#! \note{
#!   You can not assign NULL to an active finalizer using \code{finalizer<-} because this would not stop R's finalization mechanism and would carry the risk of assiging MULTIPLE finalization tasks.
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{ff}}, \code{\link{finalize}}, \code{\link[base]{reg.finalizer}}
#! }
#! \examples{
#!   x <- ff(1:12, pattern="./finalizerdemo")
#!   fnam <- filename(x)
#!   finalizer(x)
#!   finalizer(x) <- "delete"
#!   finalizer(x)
#!   rm(x)
#!   file.exists(fnam)
#!   gc()
#!   file.exists(fnam)
#! }
#! \keyword{ IO }

finalizer.ff <- function(x, ...){
  attr(attr(x, "physical"), "finalizer")
}

"finalizer<-.ff" <- function(x, ..., value){
  physical <- attr(x, "physical")
  oldfin <- attr(physical, "finalizer")
  if (is.null(value)){
    if (!is.null(oldfin))
      stop("Active finalizer '", oldfin, "'can only be REPLACED by a new name of a finalizer function, but NOT de-activated")
  }else{
    attr(attr(x, "physical"), "finalizer") <- value
    if (is.null(oldfin)){
      reg.finalizer(physical, finalize.ff_pointer, onexit=attr(physical, "finonexit"))
    }
  }
  x
}

#! \name{finalize}
#! \Rdversion{1.1}
#! \alias{finalize}
#! \alias{finalize.ff_pointer}
#! \alias{finalize.ff}
#! \alias{finalize.ffdf}
#! \title{
#!   Call finalizer
#! }
#! \description{
#!   This calls the currently assigned finalizer, either via R's finalization mechanism or manually.
#! }
#! \usage{
#! finalize(x, ...)
#! \method{finalize}{ff_pointer}(x, ...)
#! \method{finalize}{ff}(x, ...)
#! \method{finalize}{ffdf}(x, ...)
#! }
#! \arguments{
#!   \item{x}{ either an \code{\link{ff}} or \code{\link{ffdf}} object or an \code{ff_pointer}, see details }
#!   \item{\dots}{ currently ignored }
#! }
#! \details{
#!   The \code{finalize.ff_pointer} method is called from R after it had been passed to \code{\link[base]{reg.finalizer}}. It will set the finalizer name to \code{NULL} and call the finalizer.
#!   \cr
#!   The \code{finalize} generic can be called manually on \code{\link{ff}} or \code{\link{ffdf}} objects. It will call the finalizer but not touch the finalizer name.
#!   \cr
#!   For more details see \code{\link{finalizer}}
#! }
#! \note{
#!   \code{finalize.ff_pointer} MUST NEVER be called manually - neither directly nor by calling the generic on an ff_pointer (could erroneously signal that there is no pending finalization lurking around)
#! }
#! \value{
#!   returns whatever the called finalizer returns, for ffdf a list with the finalization returns of each physical component is returned.
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{finalizer}}
#! }
#! \examples{
#!   x <- ff(1:12, pattern="./finalizerdemo")
#!   fnam <- filename(x)
#!   finalizer(x)
#!   is.open(x)
#!   file.exists(fnam)
#!
#!   finalize(x)
#!
#!   finalizer(x)
#!   is.open(x)
#!   file.exists(fnam)
#!
#!   delete(x)
#!   finalizer(x)
#!   is.open(x)
#!   file.exists(fnam)
#!
#!   rm(x)
#!   gc()
#! }
#! \keyword{ IO }

finalize.ff_pointer <- function(
  x     # ff_pointer
, ...   # ignored
){
  #message("R is finalizing" , attr(x, "filename"), "")
  fin <- attr(x, "finalizer")
  if (is.null(fin))
    TRUE
  else{
    attr(x, "finalizer") <- NULL
    do.call(fin, list(x))
  }
}

finalize.ff <- function(
x       # ff object
, ...   # passed to finalizer
){
  physical <- attr(x, "physical")
  fin <- attr(physical, "finalizer")
  if (is.null(fin))
    TRUE
  else
    do.call(fin, c(list(x), list(...)))
}

finalize.ffdf <- function(
x       # ff object
, ...   # passed to finalizer
){
    p <- .subset2(x, "physical")
    ret <- lapply(p, finalize, ...)
    rnam <- .subset2(x, "row.names")
    if (is.ff(rnam))
        ret <- c(row.names=finalize(rnam, ...), ret)
    ret
}



#! \name{open.ff}
#! \alias{open.ff}
#! \alias{open.ffdf}
#! \title{ Opening an ff file }
#! \description{
#!   \command{open.ff} opens an ff file, optionally marking it readonly and optionally specifying a caching scheme.
#! }
#! \usage{
#!  \method{open}{ff}(con, readonly = FALSE, pagesize = NULL, caching = NULL, assert = FALSE, \dots)
#!  \method{open}{ffdf}(con, readonly = FALSE, pagesize = NULL, caching = NULL, assert = FALSE, \dots)
#! }
#! \arguments{
#!   \item{con}{ an \code{\link{ff}} or \code{\link{ffdf}} object }
#!   \item{readonly}{ \code{readonly} }
#!   \item{pagesize}{ number of bytes to use as pagesize or NULL to take the pagesize stored in the \code{\link[=physical.ff]{physical}} attribute of the ff object, see \code{\link{getalignedpagesize}} }
#!   \item{caching}{ one of 'mmnoflush' or 'mmeachflush', see \code{\link{ff}} }
#!   \item{assert}{ setting this to TRUE will give a message if the ff was not open already }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!   ff objects will be opened automatically when accessing their content and the file is still closed.
#!   Opening ffdf objects will open all of their \code{\link[=physical.ffdf]{physical}} components including their \code{\link[=row.names.ffdf]{row.names}} if they are \code{\link{is.ff}}
#! }
#! \value{
#!   TRUE if object could be opened, FALSE if it was opened already (or NA if not all components of an ffdf returned FALSE or TRUE on opening)
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ff}}, \code{\link{close.ff}}, \code{\link{delete}}, \code{\link{deleteIfOpen}}, \code{\link{getalignedpagesize}} }
#! \examples{
#!   x <- ff(1:12)
#!   close(x)
#!   is.open(x)
#!   open(x)
#!   is.open(x)
#!   close(x)
#!   is.open(x)
#!   x[]
#!   is.open(x)
#!   y <- x
#!   close(y)
#!   is.open(x)
#!   rm(x,y); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

open.ff <- function(con
, readonly  = FALSE
, pagesize = NULL
, caching = NULL
, assert = FALSE
, ... # dummy to keep R CMD check quiet
)
{
  readonly <- as.logical(readonly)
  physical <- attr(con, "physical")
  if (is.open(con)){
    if (attr(physical, "readonly")!=readonly && !assert)
      stop(paste("ff is already open with readonly=", attr(physical, "readonly"), sep=""))
    return(FALSE)
  }else if (assert){
		message("opening ff ", filename(con))
	}
  filename <- attr(physical, "filename")
  stopifnot(file.access(filename,0)==0 )
  if (!readonly && file.access(filename,2)){
    message("opening ff ", filename, " readonly")
    readonly <- TRUE
  }
  if (!is.null(pagesize)){
    attr(attr(con, "physical"), "pagesize") <- getalignedpagesize(pagesize) # C-code currently reads pagesize directly from the attribute
  }
  if (is.null(caching))
    caching <- attr(physical, "caching")
  else
    caching <- match.arg(caching, caching_schemes)
  if (.Call("open", physical, .ffmode[attr(physical, "vmode")], readonly, caching=="mmeachflush", PACKAGE="ff")){
    attr(attr(con, "physical"), "caching") <- caching
    if (is.null(attr(physical, "finalizer"))){
      attr(attr(con, "physical"), "finalizer") <- "close"
      reg.finalizer(physical, finalize.ff_pointer, onexit=attr(physical, "finonexit"))
    }
    return(TRUE)
  }else{
		stop("failed opening ff ", filename(con), "because ", geterrstr.ff(con))
    return(FALSE)
  }
}



#! \name{close.ff}
#! \alias{close.ff}
#! \alias{close.ffdf}
#! \alias{close.ff_pointer}
#! \title{ Closing ff files }
#! \description{
#!   Close frees the Memory Mapping resources and closes the ff file without deleting the file data.
#! }
#! \usage{
#! \method{close}{ff}(con, \dots)
#! \method{close}{ffdf}(con, \dots)
#! \method{close}{ff_pointer}(con, \dots)
#! }
#! \arguments{
#!   \item{con}{ an open ff object }
#!   \item{\dots}{ \code{\dots} }
#! }
#! \details{
#!   The \code{ff_pointer} method is not intended for manual use, it is used at finalizer dispatch time.
#!   Closing ffdf objects will close all of their \code{\link[=physical.ffdf]{physical}} components including their \code{\link[=row.names.ffdf]{row.names}} if they are \code{\link{is.ff}}
#! }
#! \value{
#!   TRUE if the file could be closed, FALSE if it was closed already (or NA if not all components of an ffdf returned FALSE or TRUE on closing)
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ff}}, \code{\link{open.ff}}, \code{\link{delete}}, \code{\link{deleteIfOpen}} }
#! \examples{
#!   x <- ff(1:12)
#!   close(x)
#!   x
#!   open(x)
#!   x
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

# version to which the finalizer dispatches
close.ff_pointer <- function(con
, ... # dummy to keep R CMD check quiet
)
{
  .Call("delete", con, PACKAGE="ff")  # this is intentionally 'delete' which releases (almost) all ressources (residual C++ RAM ist released when rm(ff) removes the external pointer (said Daniel Adler, 15.11.2007)
}
# version for manual use
close.ff <- function(con
, ... # dummy to keep R CMD check quiet
)
{
  # we do NOT touch the finalizer name if called manually (could be a temporary ff file with a 'delete' finalizer) !!
  .Call("delete", attr(con, "physical"), PACKAGE="ff")  # this is intentionally 'delete' which releases (almost) all ressources (residual C++ RAM ist released when rm(ff) removes the external pointer (said Daniel Adler, 15.11.2007)
}


#! \name{delete}
#! \alias{delete}
#! \alias{delete.ff}
#! \alias{delete.ffdf}
#! \alias{delete.ff_pointer}
#! \alias{delete.default}
#! \alias{deleteIfOpen}
#! \alias{deleteIfOpen.ff}
#! \alias{deleteIfOpen.ff_pointer}
#! \title{ Deleting the file behind an ff object }
#! \description{
#!   The generic \command{delete} deletes the content of an object without removing the object itself.
#!   The generic \command{deleteIfOpen} does the same, but only if \code{\link{is.open}} returns TRUE.
#! }
#! \usage{
#! delete(x, \dots)
#! deleteIfOpen(x, \dots)
#! \method{delete}{ff}(x, \dots)
#! \method{delete}{ffdf}(x, \dots)
#! \method{delete}{ff_pointer}(x, \dots)
#! \method{delete}{default}(x, \dots)
#! \method{deleteIfOpen}{ff}(x, \dots)
#! \method{deleteIfOpen}{ff_pointer}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an ff or ram object }
#!   \item{\dots}{ further arguments (not used) }
#! }
#! \details{
#!   The proper sequence to fully delete an ff object is: \code{delete(x);rm(x)}, where \command{delete.ff} frees the Memory Mapping resources and deletes the ff file,
#!   leaving intact the R-side object including its \code{\link{class}}, \code{\link[=physical.ff]{physical}} and \code{\link[=physical.ff]{virtual}} attributes.
#!   The default method is a compatibility function doing something similar with ram objects: by assiging an empty list to the name of the ram object to the parent frame
#!   we destroy the content of the object, leaving an empty stub that prevents raising an error if the parent frame calls the \code{delete(x);rm(x)} sequence. \cr
#!   The \command{deleteIfOpen} does the same as \command{delete} but protects closed ff objects from deletion, it is mainly intended for use through a finalizer, as are the \code{ff_pointer} methods.
#! }
#! \value{
#!   \command{delete} returns TRUE if the/all ff files could be removed and FALSE otherwise. \cr
#!   \command{deleteIfOpen} returns TRUE if the/all ff files could be removed, FALSE if not and NA if the ff object was open.
#! }
#! \author{ Jens Oehlschlägel }
#! \note{
#!   Deletion of ff files can be triggerd automatically via three routes:
#!   \enumerate{
#!     \item if an ff object with a 'delete' finalizer is removed
#!     \item if an ff object was created with \code{fffinonexit=TRUE} the finalizer is also called when R shuts down.
#!     \item if an ff object was created in \code{getOption("fftempdir")}, it will be unlinked together with the fftempdir \code{\link[base:ns-hooks]{.onUnload}}
#!   }
#!   Thus in order to retain an ff file, one has to create it elsewhere than in fftempdir with a finalizer that does not destroy the file (by default files outside fftempdir get a 'close' finalizer') i.e. one of the following:
#!   \enumerate{
#!     \item name the file AND use \code{fffinalizer="close"}
#!     \item name the file AND use \code{fffinalizer="deleteIfOpen"} AND close the ff object before leaving R
#!     \item name the file AND use \code{fffinalizer="delete"} AND use \code{fffinonexit=FALSE}
#!   }
#! }
#! \seealso{ \code{\link{ff}}, \code{\link{close.ff}}, \code{\link{open.ff}}, \code{\link[base]{reg.finalizer}} }
#! \examples{
#!   message('create the ff file outside getOption("fftempir"), 
#!     it will have default finalizer "close", so you need to delete it explicitely')
#!   x <- ff(1:12, pattern="./ffexample")
#!   delete(x)
#!   rm(x)
#! }
#! \keyword{ IO }
#! \keyword{ data }

# version to which the finalizer dispatches
delete.ff_pointer <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  .Call("delete", x, PACKAGE="ff")
  all(file.remove(filename(x)))        # filename() might be a vector of files in the future
}
# version for manual use
delete.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  physical <- attr(x, "physical")
  # make sure 'delete' is not called a second time by assigning a harmless dummy finalizer name which signals that there is still a call to finalize.ff_pointer lurking around
  if (!is.null(attr(physical, "finalizer")))
    attr(physical, "finalizer") <- "close"
  .Call("delete", physical, PACKAGE="ff")
  all(file.remove(attr(physical, "filename"))) # filename() might be a vector of files in the future
}
delete.default <- function(x
, ... # dummy to keep R CMD check quiet
){
  assign(deparse(substitute(x)), list(), parent.frame())  # delete memory associated with x but leave stub (because often we have delete(x);rm(x) in the parent frame)
}


# version to which the finalizer dispatches
deleteIfOpen.ff_pointer <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  if (is.open(x)){
    .Call("delete", x, PACKAGE="ff")
    all(file.remove(filename(x)))   # filename() might be a vector of files in the future
  }else{
    NA
  }
}


# version for manual use
deleteIfOpen.ff <- function(x
, ... # dummy to keep R CMD check quiet
)
{
  if (is.open(x)){
    delete(x)
  }else{
    NA
  }
}


if (FALSE){
  deleteIfOpen.default <- function(x
  , ... # dummy to keep R CMD check quiet
  )
  {
    message("--- Here deleteIfOpen.default ---")
    message("--- x ---")
    print(x)
    str(x)
    message("--- ... ---")
    print(list(...))
    str(list(...))
    message("------")
  }
}



# --- ff read / write / readwrite ----------------------------------------------------------

# the getset / readwrite / swap functions are efficient for read/write at once (they do maintain na.count() if na.count() is activated )
# the get,set / read,write / [,[<- functions are simplified clones of these, the write versions stop if na.count() is activated

#! \name{getset.ff}
#! \alias{getset.ff}
#! \alias{get.ff}
#! \alias{set.ff}
#! \title{ Reading and writing vectors of values (low-level) }
#! \description{
#!   The three functions \command{get.ff}, \command{set.ff} and \command{getset.ff} provide the simplest interface to access an ff file: getting and setting vector of values identified by positive subscripts
#! }
#! \usage{
#! get.ff(x, i)
#! set.ff(x, i, value, add = FALSE)
#! getset.ff(x, i, value, add = FALSE)
#! }
#! \arguments{
#!   \item{x}{ an ff object }
#!   \item{i}{ an index position within the ff file }
#!   \item{value}{ the value to write to position i }
#!   \item{add}{ TRUE if the value should rather increment than overwrite at the index position }
#! }
#! \details{
#!   \command{getset.ff} combines the effects of \command{get.ff} and \command{set.ff} in a single operation: it retrieves the old value at position \code{i} before changing it.
#!   \command{getset.ff} will maintain \code{\link{na.count}}.
#! }
#! \value{
#!   \command{get.ff} returns a vector, \command{set.ff} returns the 'changed' ff object (like all assignment functions do) and \command{getset.ff} returns the value at the subscript positions.
#!   More precisely \code{getset.ff(x, i, value, add=FALSE)} returns the old values at the subscript positions \code{i} while \code{getset.ff(x, i, value, add=TRUE)} returns the incremented values at the subscript positions.
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ \command{get.ff}, \command{set.ff} and \command{getset.ff} are low level functions that do not support \code{ramclass} and \code{ramattribs} and thus will not give the expected result with \code{factor} and \code{POSIXct} }
#! \seealso{ \code{\link{readwrite.ff}} for low-level access to contiguous chunks and \code{\link{[.ff}} for high-level access }
#! \examples{
#!  x <- ff(0, length=12)
#!  get.ff(x, 3L)
#!  set.ff(x, 3L, 1)
#!  x
#!  set.ff(x, 3L, 1, add=TRUE)
#!  x
#!  getset.ff(x, 3L, 1, add=TRUE)
#!  getset.ff(x, 3L, 1)
#!  x
#!  rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }


getset.ff <- function(x, i, value, add=FALSE)
{
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")
  if (is.double(i))
    i <- as.integer(i)
  if (!is.integer(i) || i<1 || i>length(x)) stop("illegal index")
  if(!is.null(vw(x))) stop("please use '[' to access ff with vw")

  nc <- na.count(x)
  if (!is.na(nc))
    new.nc <- is.na(value)
  vmode <- vmode(x)
  if (add)
    ret <- .Call("addgetset_vec", .ffmode[vmode], attr(x, "physical"), i, length(i), as.vmode(value, vmode), PACKAGE="ff")
  else
    ret <- .Call("getset_vec", .ffmode[vmode], attr(x, "physical"), i, length(i), as.vmode(value, vmode), PACKAGE="ff")
  if (!is.na(nc)){
    old.nc <- is.na(ret)
    na.count(x) <- nc - old.nc + new.nc
  }
  ret
}

get.ff   <- function(x, i)
{
	open(x, assert=TRUE)
  if (is.double(i))
    i <- as.integer(i)
  if (!is.integer(i) || i<1 || i>length(x)) stop("illegal index")
  if(!is.null(vw(x))) stop("please use '[' to access ff with vw")

  .Call("get_vec", .ffmode[vmode(x)], attr(x, "physical"), i, length(i), PACKAGE="ff")
}

set.ff   <- function(x, i, value, add=FALSE)
{
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")
  if (is.double(i))
    i <- as.integer(i)
  if (!is.integer(i) || i<1 || i>length(x)) stop("illegal index")
  if(!is.null(vw(x))) stop("please use '[' to access ff with vw")

  if(!is.null(physical(x)$na.count)) stop("use readwrite.ff instead to maintain na.count (or deactivate na.count(x)<-NULL)")
  vmode <- vmode(x)
  if (add)
    attr(x, "physical") <- .Call("addset_vec", .ffmode[vmode], attr(x, "physical"), i, length(i), as.vmode(value, vmode), PACKAGE="ff")
  else
    attr(x, "physical") <- .Call("set_vec", .ffmode[vmode], attr(x, "physical"), i, length(i), as.vmode(value, vmode), PACKAGE="ff")
  x
}

"[[.ff" <- function(x, i){
  if (length(i)!=1L)
    stop("i must have length 1")
  set.ff(x=x, i=i)
}

"[[<-.ff" <- function(x, i, add=FALSE, value){
  if (length(i)!=1L)
    stop("i must have length 1")
  if (length(value)!=1L)
    stop("value must have length 1")
  set.ff(x=x, i=i, value=value, add=add)
}


#! \name{readwrite.ff}
#! \alias{read.ff}
#! \alias{write.ff}
#! \alias{readwrite.ff}
#! \title{ Reading and writing vectors (low-level) }
#! \description{
#!   Simpe low-level interface for reading and writing vectors from ff files.
#! }
#! \usage{
#! read.ff(x, i, n)
#! write.ff(x, i, value, add = FALSE)
#! readwrite.ff(x, i, value, add = FALSE)
#! }
#! \arguments{
#!   \item{x}{ an ff object }
#!   \item{i}{ a start position in the ff file }
#!   \item{n}{ number of elements to read }
#!   \item{value}{ vector of elements to write }
#!   \item{add}{ TRUE if the values should rather increment than overwrite at the target positions }
#! }
#! \details{
#!   \command{readwrite.ff} combines the effects of \command{read.ff} and \command{write.ff} in a single operation: it retrieves the old values starting from position \code{i} before changing them.
#!   \command{getset.ff} will maintain \code{\link{na.count}}.
#! }
#! \value{
#!   \command{read.ff} returns a vector of values, \command{write.ff} returns the 'changed' ff object (like all assignment functions do) and \command{readwrite.ff} returns the values at the target position.
#!   More precisely \code{readwrite.ff(x, i, value, add=FALSE)} returns the old values at the position \code{i} while \code{readwrite.ff(x, i, value, add=TRUE)} returns the incremented values of \code{x}.
#! }
#! \author{ Jens Oehlschlägel }
#! \note{ \command{read.ff}, \command{write.ff} and \command{readwrite.ff} are low level functions that do not support \code{ramclass} and \code{ramattribs} and thus will not give the expected result with \code{factor} and \code{POSIXct} }
#! \seealso{ \code{\link{getset.ff}} for low-level scalar access and \code{\link{[.ff}} for high-level access }
#! \examples{
#!   x <- ff(0, length=12)
#!   read.ff(x, 3, 6)
#!   write.ff(x, 3, rep(1, 6))
#!   x
#!   write.ff(x, 3, rep(1, 6), add=TRUE)
#!   x
#!   readwrite.ff(x, 3, rep(1, 6), add=TRUE)
#!   readwrite.ff(x, 3, rep(1, 6))
#!   x
#!   rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }


readwrite.ff <- function(x, i, value, add=FALSE)
{
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")
  if(!is.null(vw(x))) stop("please use '[' to access ff with vw")

  n <- length(value)
  stopifnot( 0 < i && i+n-1 <= length(x) )

  nc <- na.count(x)
  if (!is.na(nc))
    new.nc <- sum(is.na(value))
  vmode <- vmode(x)
  if (add){
    ret <- .Call("addgetset_contiguous", .ffmode[vmode], attr(x, "physical"), as.integer(i), n, as.vmode(value, vmode), PACKAGE="ff")
  }else{
    ret <- .Call("getset_contiguous", .ffmode[vmode], attr(x, "physical"), as.integer(i), n, as.vmode(value, vmode), PACKAGE="ff")
  }
  if (!is.na(nc)){
    old.nc <- sum(is.na(ret))
    na.count(x) <- nc - old.nc + new.nc
  }
  ret
}

read.ff <- function(x, i, n)
{
	open(x, assert=TRUE)
  if(!is.null(vw(x))) stop("please use '[' to access ff with vw")

  stopifnot( 0 < i && i+n-1 <= length(x) )
  .Call("get_contiguous", .ffmode[vmode(x)], attr(x, "physical"), as.integer(i), as.integer(n), PACKAGE="ff")
}

write.ff <- function(x, i, value, add=FALSE)
{
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")
  if(!is.null(vw(x))) stop("please use '[' to access ff with vw")

  if(!is.null(physical(x)$na.count)) stop("use readwrite.ff instead to maintain na.count (or deactivate na.count(x)<-NULL)")

  n <- length(value)
  stopifnot( 0 < i && i+n-1 <= length(x) )

  vmode <- vmode(x)
  if (add)
    attr(x, "physical") <- .Call("addset_contiguous", .ffmode[vmode], attr(x, "physical"), as.integer(i), n, as.vmode(value, vmode), PACKAGE="ff")
  else
    attr(x, "physical") <- .Call("set_contiguous", .ffmode[vmode], attr(x, "physical"), as.integer(i), n, as.vmode(value, vmode), PACKAGE="ff")
  x
}




#! \name{swap}
#! \alias{swap}
#! \alias{swap.ff}
#! \alias{swap.ff_array}
#! \alias{swap.default}
#! \title{ Reading and writing in one operation (high-level) }
#! \description{
#!   The generic \command{swap} combines \code{x[i]} and \code{x[i] <- value} in a single operation.
#! }
#! \usage{
#! swap(x, value, \dots)
#! \method{swap}{ff}(x, value, i, add = FALSE, pack = FALSE, \dots)
#! \method{swap}{ff_array}(x, value, \dots, bydim = NULL, drop = getOption("ffdrop"), add = FALSE, pack = FALSE)
#! \method{swap}{default}(x, value, \dots, add = FALSE)
#! }
#! \arguments{
#!   \item{x}{ a ff or ram object }
#!   \item{value}{ the new values to write, possibly recycled, see \code{\link{[.ff}} }
#!   \item{i}{ index information, see \code{\link{[.ff}} }
#!   \item{\dots}{ missing OR up to length(dim(x)) index expressions OR (ff only) \code{\link{hi}} objects }
#!   \item{drop}{ logical scalar indicating whether array dimensions shall be dropped }
#!   \item{bydim}{ how to interpret vector to array data, see \code{\link{[.ff}} }
#!   \item{add}{ TRUE if the values should rather increment than overwrite at the target positions, see \code{\link{readwrite.ff}} }
#!   \item{pack}{ FALSE to prevent rle-packing in hybrid index preprocessing, see \code{\link{as.hi}} }
#! }
#! \details{
#!   \preformatted{
#!   y <- swap(x, value, i, add=FALSE, ...)
#!
#!   is a shorter and more efficient version of
#!
#!   y <- x[i, add=FALSE, ...]
#!   x[i, add=FALSE, ...] <- value
#!
#!   and
#!
#!   y <- swap(x, value, i, add=TRUE, ...)
#!
#!   is a shorter and more efficient version of
#!
#!   y <- x[i, add=TRUE, ...]
#!   y <- y + value
#!   x[i, add=FALSE, ...] <- y
#!   }
#! }
#! \value{
#!   Values at the target positions.
#!   More precisely \code{swap(x, value, i, add=FALSE)} returns the old values at the position \code{i} while \code{swap(x, value, i, add=TRUE)} returns the incremented values of \code{x}.
#! }
#! \author{ Jens Oehlschlägel }
#! \note{
#!   Note that \code{swap.default} changes the object in its parent frame and thus violates R's usual functional programming logic.
#!   When using \code{add=TRUE}, duplicated index positions should be avoided, because ff and ram objects behave differently:
#!   \preformatted{
#!   swap.ff(x, 1, c(3,3), add=TRUE)
#!   # will increment x at position 3 TWICE by 1, while
#!   swap.default(x, 1, c(3,3), add=TRUE)
#!   # will increment x at position 3 just ONCE by 1
#!   }
#! }
#! \seealso{ \code{\link{[.ff}}, \code{\link{add}}, \code{\link{readwrite.ff}}, \code{\link{getset.ff}}, \code{\link{LimWarn}} }
#! \examples{
#!   x <- ff("a", levels=letters, length=52)
#!   y <- swap(x, "b", sample(length(x), 26))
#!   x
#!   y
#!   rm(x,y); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ array }


#! \name{Extract.ff}
#! \alias{Extract.ff}
#! \alias{[.ff}
#! \alias{[<-.ff}
#! \alias{[.ff_array}
#! \alias{[<-.ff_array}
#! \alias{[[.ff}
#! \alias{[[<-.ff}
#! \title{ Reading and writing vectors and arrays (high-level) }
#! \description{
#!   These are the main methods for reading and writing data from ff files.
#! }
#! \usage{
#! \method{[}{ff}(x, i, pack = FALSE)
#! \method{[}{ff}(x, i, add = FALSE, pack = FALSE) <- value
#! \method{[}{ff_array}(x, \dots, bydim = NULL, drop = getOption("ffdrop"), pack = FALSE)
#! \method{[}{ff_array}(x, \dots, bydim = NULL, add = FALSE, pack = FALSE) <- value
#! \method{[[}{ff}(x, i)
#! \method{[[}{ff}(x, i, add = FALSE) <- value
#! }
#! \arguments{
#!   \item{x}{ an ff object }
#!   \item{i}{ missing OR a single index expression OR a \code{\link{hi}} object }
#!   \item{\dots}{ missing OR up to length(dim(x)) index expressions OR \code{\link{hi}} objects }
#!   \item{drop}{ logical scalar indicating whether array dimensions shall be dropped }
#!   \item{bydim}{ the dimorder which shall be used in interpreting vector to/from array data }
#!   \item{pack}{ FALSE to prevent rle-packing in hybrid index preprocessing, see \code{\link{as.hi}} }
#!   \item{value}{ the values to be assigned, possibly recycled }
#!   \item{add}{ TRUE if the values should rather increment than overwrite at the target positions, see \code{\link{readwrite.ff}} }
#! }
#! \details{
#!   The single square bracket operators \command{[} and \command{[<-} are the workhorses for accessing the content of an ff object.
#!   They support \code{ff_vector} and \code{ff_array} access (\code{\link{dim.ff}}), they respect virtual windows (\code{\link{vw}}),
#!   \code{\link{names.ff}} and \code{\link{dimnames.ff}} and retain \code{\link{ramclass}} and \code{\link{ramattribs}}
#!   and thus support \code{\link{POSIXct}} and \code{\link{factor}}, see \code{\link{levels.ff}}.
#!   \cr
#!   The functionality of \command{[} and \command{[<-} cn be combined into one efficient operation, see \code{\link{swap}}.
#!   \cr
#!   The double square bracket operator \command{[[} is a shortcut for \code{\link{get.ff}}
#!   resp. \code{\link{set.ff}}, however, you should not rely on this for the future, see \code{\link{LimWarn}}. For programming please prefer \command{[}.
#! }
#! \section{Index expressions}{
#!   \code{x <- ff(1:12, dim=c(3,4), dimnames=list(letters[1:3], NULL))}
#!   \tabular{rll}{
#!   \emph{allowed expression}    \tab -- \tab \emph{\code{example}}                \cr
#!    positive integers           \tab    \tab \code{x[ 1 ,1]}                      \cr
#!    negative integers           \tab    \tab \code{x[ -(2:12) ]}                  \cr
#!    logical                     \tab    \tab \code{x[ c(TRUE, FALSE, FALSE) ,1]}  \cr
#!    character                   \tab    \tab \code{x[ "a" ,1]}                    \cr
#!    integer matrices            \tab    \tab \code{x[ rbind(c(1,1)) ]}            \cr
#!    hybrid index                \tab    \tab \code{x[ hi ,1]}                     \cr
#!   \emph{disallowed expression} \tab -- \tab \emph{\code{example}}                \cr
#!    zeros                       \tab    \tab \code{x[ 0 ]}                        \cr
#!    NAs                         \tab    \tab \code{x[ NA ]}                       \cr
#!   }
#! }
#! \section{Dimorder and bydim}{
#!   Arrays in R have always standard \code{\link{dimorder} 1:length(dim(x))} while ff allows to store an array in a different dimorder.
#!   Using nonstandard dimorder (see \code{\link{dimorderStandard}}) can speed up certain access operations: while matrix \code{dimorder=c(1,2)} -- column-major order --
#!   allows fast extraction of columns, \code{dimorder=c(2,1)} allows fast extraction of rows.
#!   \cr
#!   While the dimorder -- being an attribute of an \code{ff_array} -- controls how the vector in an ff file is interpreted,
#!   the \code{bydim} argument to the extractor functions controls, how assigment vector values
#!   in \command{[<-} are translated to the array and how the array is translated to a vector in \command{[} subscripting.
#!   Note that \code{bydim=c(2,1)} corresponds to \code{matrix(..., byrow=TRUE)}.
#! }
#! \section{Multiple vector interpretation in arrays }{
#!   In case of non-standard dimorder (see \code{\link{dimorderStandard}})
#!   the vector sequence of array elements in R and in the ff file differs.
#!   To access array elements in file order, you can use \code{\link{getset.ff}}, \code{\link{readwrite.ff}}
#!   or copy the ff object and set \code{dim(ff)<-NULL} to get a vector view into the ff object
#!   (using \code{[} dispatches the vector method \code{\link{[.ff}}).
#!   To access the array elements in R standard dimorder you simply use \code{[} which dispatches
#!   to \code{\link{[.ff_array}}. Note that in this case \code{\link{as.hi}} will unpack the complete index, see next section.
#! }
#! \section{RAM expansion of index expressions}{
#!   Some index expressions do not consume RAM due to the \code{\link{hi}} representation,
#!   for example \code{1:n} will almost consume no RAM hoewever large n.
#!   However, some index expressions are expanded and require to \code{\link{maxindex}(i) * .rambytes["integer"]} bytes,
#!   either because the sorted sequence of index positions cannot be rle-packed efficiently
#!   or because \code{\link{hiparse}} cannot yet parse such expression and falls back to evaluating/expanding the index expression.
#!   If the index positions are not sorted, the index will be expanded and a second vector is needed to store the information for re-ordering,
#!   thus the index requires \code{2 * \link{maxindex}(i) * .rambytes["integer"]} bytes.
#! }
#! \section{RAM expansion when recycling assigment values}{
#!   Some assignment expressions do not consume RAM for recycling, for example \code{x[1:n] <- 1:k}
#!   will not consume RAM hoewever large n compared to k, when x has standard \code{\link{dimorder}}.
#!   However, if \code{length(value)>1}, assignment expressions with non-ascending index positions trigger recycling the value R-side to the full index length.
#!   This will happen if \code{\link{dimorder}} does not match parameter \code{bydim} or if the index is not sorted ascending.
#! }
#! \value{
#!   The read operators \command{[} and \command{[[} return data from the ff object,
#!   possibly decorated with \code{\link[ff:names.ff]{names}}, \code{\link[ff:dim.ff]{dim}},
#!   \code{\link[ff:dimnames.ff]{dimnames}} and further attributes and classes (see \code{\link{ramclass}}, \code{\link{ramattribs}})
#!   \cr
#!   The write operators \command{[<-} and \command{[[<-} return the 'modified' ff object (like all assignment operators do).
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ff}}, \code{\link{swap}}, \code{\link{add}}, \code{\link{readwrite.ff}}, \code{\link{LimWarn}} }
#! \examples{
#!    message("look at different dimorders")
#!    x <- ff(1:12, dim=c(3,4), dimorder=c(1,2))
#!    x[]
#!    as.vector(x[])
#!    x[1:12]
#!    x <- ff(1:12, dim=c(3,4), dimorder=c(2,1))
#!    x[]
#!    as.vector(x[])
#!    message("Beware (might be changed)")
#!    x[1:12]
#!
#!    message("look at different bydim")
#!    matrix(1:12, nrow=3, ncol=4, byrow=FALSE)
#!    x <- ff(1:12, dim=c(3,4), bydim=c(1,2))
#!    x
#!    matrix(1:12, nrow=3, ncol=4, byrow=TRUE)
#!    x <- ff(1:12, dim=c(3,4), bydim=c(2,1))
#!    x
#!    x[,, bydim=c(2,1)]
#!    as.vector(x[,, bydim=c(2,1)])
#!    message("even consistent interpretation of vectors in assignments")
#!    x[,, bydim=c(1,2)] <- x[,, bydim=c(1,2)]
#!    x
#!    x[,, bydim=c(2,1)] <- x[,, bydim=c(2,1)]
#!    x
#!    rm(x); gc()
#!
#!   \dontrun{
#!    message("some performance implications of different dimorders")
#!    n <- 100
#!    m <- 100000
#!    a <- ff(1L,dim=c(n,m))
#!    b <- ff(1L,dim=c(n,m), dimorder=2:1)
#!    system.time(lapply(1:n, function(i)sum(a[i,])))
#!    system.time(lapply(1:n, function(i)sum(b[i,])))
#!    system.time(lapply(1:n, function(i){i<-(i-1)*(m/n)+1; sum(a[,i:(i+m/n-1)])}))
#!    system.time(lapply(1:n, function(i){i<-(i-1)*(m/n)+1; sum(b[,i:(i+m/n-1)])}))
#!
#!    n <- 100
#!    a <- ff(1L,dim=c(n,n,n,n))
#!    b <- ff(1L,dim=c(n,n,n,n), dimorder=4:1)
#!    system.time(lapply(1:n, function(i)sum(a[i,,,])))
#!    system.time(lapply(1:n, function(i)sum(a[,i,,])))
#!    system.time(lapply(1:n, function(i)sum(a[,,i,])))
#!    system.time(lapply(1:n, function(i)sum(a[,,,i])))
#!    system.time(lapply(1:n, function(i)sum(b[i,,,])))
#!    system.time(lapply(1:n, function(i)sum(b[,i,,])))
#!    system.time(lapply(1:n, function(i)sum(b[,,i,])))
#!    system.time(lapply(1:n, function(i)sum(b[,,,i])))
#!
#!    n <- 100
#!    m <- 100000
#!    a <- ff(1L,dim=c(n,m))
#!    b <- ff(1L,dim=c(n,m), dimorder=2:1)
#!    system.time(ffrowapply(sum(a[i1:i2,]), a, RETURN=TRUE, CFUN="csum", BATCHBYTES=16104816\%/\%20))
#!    system.time(ffcolapply(sum(a[,i1:i2]), a, RETURN=TRUE, CFUN="csum", BATCHBYTES=16104816\%/\%20))
#!    system.time(ffrowapply(sum(b[i1:i2,]), b, RETURN=TRUE, CFUN="csum", BATCHBYTES=16104816\%/\%20))
#!    system.time(ffcolapply(sum(b[,i1:i2]), b, RETURN=TRUE, CFUN="csum", BATCHBYTES=16104816\%/\%20))
#!
#!    rm(a,b); gc()
#!   }
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ array }


swap.ff <- function(
  x
, value
, i
, add   = FALSE
, pack  = FALSE
, ... # dummy to keep R CMD check quiet
){
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")

  vmode <- vmode(x)
  vw <- vw(x)
  fflen <- length(x)
  fflev <- levels(x)
  ffnam <- names(x)
  if (fflen){
    if (missing(i)){
      index <- hi(from=1, to=fflen, maxindex=fflen, vw=vw, pack=pack)
      nreturn <- fflen
    }else if(is.ff(i)){
      stop("ff subscripts not implemented for swap.ff")
    }else{
      index <- as.hi(substitute(i), maxindex=fflen, vw=vw, pack=pack, names=ffnam, envir=parent.frame())
      nreturn <- poslength(index)
    }
  }else{
    nreturn <- 0:0
  }
  if (nreturn){
    nc <- na.count(x)
    ixre <- !is.null(index$ix) || index$re
    if (is.null(fflev))
      value <- as.vmode(value, vmode)
    else
      value <- ram2ffcode(value, fflev, vmode)
    nvalue <- length(value)
    if (nvalue>1){
      if (nreturn!=nvalue){
        rb <- nreturn%%nvalue
        if (!is.na(nc))
          new.nc <- nreturn%/%nvalue * sum(is.na(value)) + if (rb) sum(is.na(value[1:rb])) else 0L
        if (rb)
          warning("number of elements to replace is not multiple of values for replacement")
        if (ixre){
          # if possible we try to recycle on the C-side, but in case of ix or re we need to recycle already here
          value <- rep(value, length=nreturn)
          nvalue <- nreturn
        }
      }else{
        if (!is.na(nc))
          new.nc <- sum(is.na(value))
      }
      if (!is.null(index$ix))
        value <- value[index$ix]
      else if (index$re)
        value <- rev(value)
    }else{
      if(nvalue<1)
        stop("no value for replacement")
      if (!is.na(nc))
        new.nc <- nreturn * is.na(value)
    }
    if (add)
      ret <- .Call("addgetset_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
    else
      ret <- .Call("getset_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
    ret <- unsort.hi(ret, index)
    if (!is.na(nc)){
      old.nc <- sum(is.na(ret))
      na.count(x) <- nc - old.nc + new.nc
    }
    if (!is.null(ffnam))
      setattr(ret, "names", ffnam[as.integer(index)])  #names(ret) <- ffnam[as.integer(index)]
  }else{
    ret <- vector(mode=.rammode[vmode], length=0)
  }
  if (!is.null(fflev)){
    if (.vunsigned[vmode])
      ret <- ret + 1L
    setattr(ret, "levels", fflev) #levels(ret) <- fflev
  }
  ramattribs <- attr(attr(x, "virtual"), "ramattribs")
  if (!is.null(ramattribs)){
    setattributes(ret, ramattribs)  #attributes(ret) <- c(attributes(ret), ramattribs)
  }
  setattr(ret, "class",attr(attr(x, "virtual"), "ramclass"))  #class(ret) <- attr(attr(x, "virtual"), "ramclass")
  .vset[[vmode]](ret)
  return(ret)
}


"[.ff" <- function(
  x
, i
, pack = FALSE
)
{

  vmode <- vmode(x)
  vw <- vw(x)
  fflen <- length(x)
  fflev <- levels(x)
  ffnam <- names(x)
  if (fflen){
    if (missing(i)){
      simple <- TRUE
      nreturn <- fflen
    }else if(is.ff(i)){
      return(ffindexget(x, i))
    }else{
      simple <- FALSE
      # note: do not ask inherits(index, "hi") directly, because this would evaluate and could create large objects in RAM
      # calling as.hi(substitute()) will first try to parse, this will fail on hi input and then only evaluate hi
      # this is a bit indirect but can substantially save memory
      index <- as.hi(substitute(i), maxindex=fflen, vw=vw, pack=pack, names=ffnam, envir=parent.frame())
      nreturn <- poslength(index)
    }
  }else{
    nreturn <- 0:0
  }
  if (nreturn){
		open(x, assert=TRUE)
    if (simple){
      ret <- .Call("get_contiguous", .ffmode[vmode], attr(x, "physical"), if (is.null(vw)) 1L else vw[1]+1L, nreturn, PACKAGE="ff")
      if (!is.null(ffnam))
        setattr(ret, "names", ffnam)  #names(ret) <- ffnam
    }else{
      ret <- .Call("get_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, PACKAGE="ff")
      ret <- unsort.hi(ret,index)
      if (!is.null(ffnam))
        setattr(ret, "names", ffnam[as.integer(index)])  #names(ret) <- ffnam[as.integer(index)]
    }
  }else{
    ret <- vector(mode=.rammode[vmode], length=0)
  }
  if (!is.null(fflev)){
    if (.vunsigned[vmode])
      ret <- ret + 1L
    setattr(ret, "levels", fflev)  #levels(ret) <- fflev
  }
  ramattribs <- attr(attr(x, "virtual"), "ramattribs")
  if (!is.null(ramattribs)){
    setattributes(ret, ramattribs) #attributes(ret) <- c(attributes(ret), ramattribs)
  }
  setattr(ret, "class",attr(attr(x, "virtual"), "ramclass"))  #class(ret) <- attr(attr(x, "virtual"), "ramclass")
  .vset[[vmode]](ret)

  return(ret)
}



"[<-.ff" <- function(
  x
, i
, add   = FALSE
, pack  = FALSE
, value
)
{
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")
  physical <- physical(x)
  if(!is.null(physical$na.count)) stop("use swap instead to maintain na.count (or deactivate na.count(x)<-NULL for assigning)")
  if(!is.null(physical$is.sorted)) stop("deactivate is.sorted(x)<-NULL for assigning)")

  vmode <- vmode(x)
  vw <- vw(x)
  fflen <- length(x)
  fflev <- levels(x)
  ffnam <- names(x)
  if (fflen){
    if (missing(i)){
      index <- hi(from=1, to=fflen, maxindex=fflen, vw=vw, pack=pack)
      nreturn <- fflen
    }else if(is.ff(i)){
      return(ffindexset(x, i, value))
    }else{
      index <- as.hi(substitute(i), maxindex=fflen, vw=vw, pack=pack, names=ffnam, envir=parent.frame())
      nreturn <- poslength(index)
    }
  }else{
    nreturn <- 0:0
  }
  if (nreturn){
    ixre <- !is.null(index$ix) || index$re
    vmode <- vmode(x)
    if (is.null(fflev))
      value <- as.vmode(value, vmode)
    else
      value <- ram2ffcode(value, fflev, vmode)
    nvalue <- length(value)
    if (nvalue>1){
      if (nreturn!=nvalue){
        if (nreturn%%nvalue)
          warning("number of elements to replace is not multiple of values for replacement")
        if (ixre){
          # if possible we try to recycle on the C-side, but in case of ix or re we need to recycle already here
          value <- rep(value, length=nreturn)
          nvalue <- nreturn
        }
      }
      if (!is.null(index$ix))
        value <- value[index$ix]
      else if (index$re)
        value <- rev(value)
    }else if(nvalue<1){
      stop("no value for replacement")
    }
    if (add)
      attr(x, "physical") <- .Call("addset_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
    else
      attr(x, "physical") <- .Call("set_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
  }
  return(x)
}


swap.ff_array <- function(
  x
, value
, ...
, bydim = NULL
, drop  = getOption("ffdrop")
, add   = FALSE
, pack  = FALSE
)
{
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")

  vmode <- vmode(x)
  vw <- vw(x)
  fflen <- length(x)
  fflev <- levels(x)
  ffdim <- dim(x)
  ffdimnam <- dimnames(x)
  ffdimord <- dimorder(x)
  ndim <- length(ffdim)

  arguments <- match.call(expand.dots=FALSE)[["..."]]
  narguments <- length(arguments)

  if (narguments==0 || (narguments==1 && is.language(arguments[[1]]) && arguments[[1]]=="")){  # ad[]
    index <- lapply(1:ndim, function(d){
      hi(from=1, to=ffdim[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack)
    })
    indexdim <- ffdim
    nreturn <- fflen
    simple <- FALSE
    ixre <- FALSE
    drop <- FALSE  # R always had: no subsetting no drop
  }else if (narguments==1){                                     # ad[cbind()] ad[1:n] 1d[1:n]
    index <- as.hi(arguments[[1]], maxindex=fflen, vw=vw, dim=ffdim, dimorder=ffdimord, pack=pack, envir=parent.frame())
    nreturn <- poslength(index)
    ixre <- !is.null(index$ix) || index$re
    if (!is.null(index$dim) || ndim>1){ # ad[cbind()] ad[1:n]
      indexdim <- NULL
      simple <- TRUE
    }else{ # 1d[1:n]
      index <- list(index)
      indexdim <- nreturn
      simple <- FALSE
    }
  }else if (narguments==ndim){
    noselect <- TRUE
    ixre <- FALSE
    envir <- parent.frame()
    index <- lapply(1:ndim, function(d){
      if (is.language(arguments[[d]]) && arguments[[d]]==""){
        hi(from=1, to=ffdim[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack)
      }else{
        noselect <<- FALSE
        lret <- as.hi(arguments[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack, names=ffdimnam[[d]], envir=envir)
        if (!is.null(lret$ix) || lret$re)
          ixre <<- TRUE
        lret
      }
    })
    indexdim <- sapply(index,poslength)
    nreturn <- as.integer(prod(indexdim))
    simple <- FALSE
    if (noselect)
      drop <- FALSE
    rm(noselect)
  }else{
    stop("wrong number of dimensions")
  }


  if (fflen && nreturn){
    nc <- na.count(x)
    if (is.null(fflev))
      value <- as.vmode(value, vmode)
    else
      value <- ram2ffcode(value, fflev, vmode)
    nvalue <- length(value)
    if (nvalue>1){
      if (!dimorderStandard(bydim)){
        value <- vector2array(value, dim=indexdim, dimorder=bydim)
      }
      if (nreturn!=nvalue){
        rb <- nreturn%%nvalue
        if (!is.na(nc))
          new.nc <- nreturn%/%nvalue * sum(is.na(value)) + if (rb) sum(is.na(value[1:rb])) else 0L
        if (rb)
          warning("number of elements to replace is not multiple of values for replacement")
        if (ixre){
          # if possible we try to recycle on the C-side, but in case of ix or re we need to recycle already here
          value <- rep(value, length=nreturn)
          nvalue <- nreturn
        }
      }else{
        if (!is.na(nc))
          new.nc <- sum(is.na(value))
      }
      if (simple){
        if (!is.null(index$ix))
          value <- value[index$ix]
        else if (index$re)
          value <- rev(value)
      }else{
        if (ixre){
          value <- array(value, dim=indexdim)
          ix <- lapply(index, function(i){
                  if (is.null(i$ix)){
                    if (i$re)
                      rev(1:length(i))
                    else
                      1:length(i)
                  }else{
                    i$ix
                  }
                })
          value <- do.call("[", c(list(value), ix))
        }else{
          ix <- NULL
        }
      }
    }else{
      if(nvalue<1)
        stop("no value for replacement")
      if (!is.na(nc))
        new.nc <- nreturn * is.na(value)
    }

    if (simple){
      if (add)
        ret <- .Call("addgetset_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
      else
        ret <- .Call("getset_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
      ret <- unsort.hi(ret, index)
    }else{
      if (!is.null(vw))
        ffdim <- as.integer(colSums(vw))
      cumffdim <- c(1L,cumprod(ffdim[ffdimord])[-ndim])
      cumindexdim <- cumprod(indexdim)
      nreturn <- cumindexdim[ndim]
      cumindexdim <- c(1,cumindexdim[-ndim])[ffdimord]

      if (add)
        ret <- .Call("addgetset_array", .ffmode[vmode], attr(x, "physical"), index[ffdimord], as.integer(indexdim[ffdimord]), as.integer(ffdim[ffdimord]), as.integer(ndim), as.integer(nreturn), as.integer(cumindexdim), as.integer(cumffdim), value, PACKAGE="ff")
      else
        ret <- .Call("getset_array", .ffmode[vmode], attr(x, "physical"), index[ffdimord], as.integer(indexdim[ffdimord]), as.integer(ffdim[ffdimord]), as.integer(ndim), as.integer(nreturn), as.integer(cumindexdim), as.integer(cumffdim), value, PACKAGE="ff")
      #ret <- array(ret, dim=indexdim)
      setattr(ret, "dim", indexdim)   #dim(ret) <- indexdim
      if (nvalue==1)
        ret <- unsort.ahi(ret, index, ixre)
      else
        ret <- unsort.ahi(ret, index, ixre, ix)
      do.bydim <- !dimorderStandard(bydim)
      if (do.bydim){
        ret <- array2vector(ret, dim=indexdim, dimorder=bydim)
        setattr(ret, "dim", indexdim[bydim])   #dim(ret) <- indexdim[bydim]
      }
      if (!is.null(ffdimnam)){
        dimnam <- 1:ndim
        setattr(dimnam, "names", names(ffdimnam))  #names(dimnam) <- names(ffdimnam)
        dimnam <- lapply(dimnam, function(d){ffdimnam[[d]][as.integer(index[[d]])]})
        if (do.bydim)
          setattr(ret, "dimnames",  dimnam[bydim]) #dimnames(ret) <- dimnam[bydim]
        else
          setattr(ret, "dimnames",  dimnam) #dimnames(ret) <- dimnam
      }
      if (drop){
        ret <- drop(ret)
        if (length(dim(ret))==1)
          ret <- as.vector(ret)
      }
    }
    if (!is.na(nc)){
      old.nc <- sum(is.na(ret))
      na.count(x) <- nc - old.nc + new.nc
    }
  }else{
    if (drop)
      ret <- vector(mode=.rammode[vmode], length=0)
    else{
      ret <- vector(mode=.rammode[vmode], length=0)
      setattr(ret, "dim", indexdim) #dim(ret) <- indexdim
    }
  }
  if (!is.null(fflev)){
    if (.vunsigned[vmode])
      ret <- ret + 1L
    setattr(ret, "levels", fflev) #levels(ret) <- fflev
  }
  ramattribs <- attr(attr(x, "virtual"), "ramattribs")
  if (!is.null(ramattribs)){
    setattributes(ret, ramattribs)  #attributes(ret) <- c(attributes(ret), ramattribs)
  }
  setattr(ret, "class",attr(attr(x, "virtual"), "ramclass"))  #class(ret) <- attr(attr(x, "virtual"), "ramclass")
  .vset[[vmode]](ret)
  return(ret)
}


"[.ff_array" <- function(
  x
, ...
, bydim = NULL
, drop  = getOption("ffdrop")
, pack  = FALSE
)
{
	open(x, assert=TRUE)
  vmode <- vmode(x)
  vw <- vw(x)
  fflen <- length(x)
  fflev <- levels(x)
  ffdim <- dim(x)
  ffdimord <- dimorder(x)
  ffdimnam <- dimnames(x)
  ndim <- length(ffdim)
  arguments <- match.call(expand.dots=FALSE)[["..."]]
  narguments <- length(arguments)
  if (narguments==0 || (narguments==1 && is.language(arguments[[1]]) && arguments[[1]]=="")){  # ad[]
    index <- lapply(1:ndim, function(d){
      hi(from=1, to=ffdim[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack)
    })
    indexdim <- ffdim
    nreturn <- fflen
    simple <- FALSE
    ixre <- FALSE
    drop <- FALSE  # R always had: no subsetting no drop
  }else if (narguments==1){                                     # ad[cbind()] ad[1:n] 1d[1:n]
    index <- as.hi(arguments[[1]], maxindex=fflen, vw=vw, dim=ffdim, dimorder=ffdimord, pack=pack, envir=parent.frame())
    nreturn <- poslength(index)
    ixre <- !is.null(index$ix) || index$re
    if (!is.null(index$dim) || ndim>1){ # ad[cbind()] ad[1:n]
      indexdim <- NULL
      simple <- TRUE
    }else{ # 1d[1:n]
      index <- list(index)
      indexdim <- nreturn
      simple <- FALSE
    }
  }else if (narguments==ndim){
    noselect <- TRUE
    ixre <- FALSE
    envir <- parent.frame()
    index <- lapply(1:ndim, function(d){
      if (is.language(arguments[[d]]) && arguments[[d]]==""){
        hi(from=1, to=ffdim[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack)
      }else{
        noselect <<- FALSE
        lret <- as.hi(arguments[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack, names=ffdimnam[[d]], envir=envir)
        if (!is.null(lret$ix) || lret$re)
          ixre <<- TRUE
        lret
      }
    })
    indexdim <- sapply(index,poslength)
    nreturn <- as.integer(prod(indexdim))
    simple <- FALSE
    if (noselect)
      drop <- FALSE
    rm(noselect)
  }else{
    stop("wrong number of dimensions")
  }
  if (fflen && nreturn){
    if (simple){
      ret <- .Call("get_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, PACKAGE="ff")
      ret <- unsort.hi(ret, index)
    }else{
      if (!is.null(vw))
        ffdim <- as.integer(colSums(vw))
      cumffdim <- c(1,cumprod(ffdim[ffdimord])[-ndim])
      cumindexdim <- cumprod(indexdim)
      nreturn <- cumindexdim[ndim]
      cumindexdim <- c(1,cumindexdim[-ndim])[ffdimord]

      ret <- .Call("get_array", .ffmode[vmode], attr(x, "physical"), index[ffdimord], as.integer(indexdim[ffdimord]), as.integer(ffdim[ffdimord]), as.integer(ndim), as.integer(nreturn), as.integer(cumindexdim), as.integer(cumffdim), PACKAGE="ff")
      #ret <- array(ret, dim=indexdim)
      setattr(ret, "dim", indexdim) #dim(ret) <- indexdim
      ret <- unsort.ahi(ret, index, ixre)
      do.bydim <- !dimorderStandard(bydim)
      if (do.bydim){
        ret <- array2vector(ret, dim=indexdim, dimorder=bydim)
        setattr(ret, "dim", indexdim[bydim]) #dim(ret) <- indexdim[bydim]
      }
      if (!is.null(ffdimnam)){
        dimnam <- 1:ndim
        #names(dimnam) <- names(ffdimnam)
        setattr(dimnam, "names", names(ffdimnam))
        dimnam <- lapply(dimnam, function(d){ffdimnam[[d]][as.integer(index[[d]])]})
        if (do.bydim)
          setattr(ret, "dimnames", dimnam[bydim]) #dimnames(ret) <- dimnam[bydim]
        else
          setattr(ret, "dimnames", dimnam) #dimnames(ret) <- dimnam
      }
      if (drop){
        ret <- drop(ret)
        if (length(dim(ret))==1)
          ret <- as.vector(ret)
      }

    }
  }else{
    if (drop)
      ret <- vector(mode=.rammode[vmode], length=0)
    else{
      ret <- vector(mode=.rammode[vmode], length=0)
      setattr(ret, "dim", indexdim) #dim(ret) <- indexdim
    }
  }
  if (!is.null(fflev)){
    if (.vunsigned[vmode])
      ret <- ret + 1L
    setattr(ret, "levels", fflev) #levels(ret) <- fflev
  }
  ramattribs <- attr(attr(x, "virtual"), "ramattribs")
  if (!is.null(ramattribs)){
    setattributes(ret, ramattribs)  #attributes(ret) <- c(attributes(ret), ramattribs)
  }
  setattr(ret, "class",attr(attr(x, "virtual"), "ramclass"))  #class(ret) <- attr(attr(x, "virtual"), "ramclass")
  .vset[[vmode]](ret)
  return(ret)
}


"[<-.ff_array" <- function(
  x
, ...
, bydim = NULL
, add   = FALSE
, pack  = FALSE
, value
)
{
	open(x, assert=TRUE)
  if( is.readonly(x) ) stop("ff is readonly")
  physical <- physical(x)
  if(!is.null(physical$na.count)) stop("use swap instead to maintain na.count (or deactivate na.count(x)<-NULL for assigning)")
  if(!is.null(physical$is.sorted)) stop("deactivate is.sorted(x)<-NULL for assigning)")

  vmode <- vmode(x)
  vw <- vw(x)
  fflen <- length(x)
  fflev <- levels(x)
  ffdim <- dim(x)
  ffdimord <- dimorder(x)
  ffdimnam <- dimnames(x)
  ndim <- length(ffdim)

  arguments <- match.call(expand.dots=FALSE)[["..."]]
  narguments <- length(arguments)
  if (narguments==0 || (narguments==1 && is.language(arguments[[1]]) && arguments[[1]]=="")){  # ad[]
    index <- lapply(1:ndim, function(d){
      hi(from=1, to=ffdim[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack)
    })
    indexdim <- ffdim
    nreturn <- fflen
    simple <- FALSE
    ixre <- FALSE
  }else if (narguments==1){                                     # ad[cbind()] ad[1:n] 1d[1:n]
    index <- as.hi(arguments[[1]], maxindex=fflen, vw=vw, dim=ffdim, dimorder=ffdimord, pack=pack, envir=parent.frame())
    nreturn <- poslength(index)
    ixre <- !is.null(index$ix) || index$re
    if (!is.null(index$dim) || ndim>1){ # ad[cbind()] ad[1:n]
      indexdim <- NULL
      simple <- TRUE
    }else{ # 1d[1:n]
      index <- list(index)
      indexdim <- nreturn
      simple <- FALSE
    }
  }else if (narguments==ndim){
    ixre <- FALSE
    envir <- parent.frame()
    index <- lapply(1:ndim, function(d){
      if (is.language(arguments[[d]]) && arguments[[d]]==""){
        hi(from=1, to=ffdim[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack)
      }else{
        lret <- as.hi(arguments[[d]], maxindex=ffdim[[d]], vw=vw[,d], pack=pack, names=ffdimnam[[d]], envir=envir)
        if (!is.null(lret$ix) || lret$re)
          ixre <<- TRUE
        lret
      }
    })
    indexdim <- sapply(index,poslength)
    nreturn <- as.integer(prod(indexdim))
    simple <- FALSE
  }else{
    stop("wrong number of dimensions")
  }

  if (fflen && nreturn){
    if (is.null(fflev))
      value <- as.vmode(value, vmode)
    else
      value <- ram2ffcode(value, fflev, vmode)
    nvalue <- length(value)
    if (nvalue>1){
      if (!dimorderStandard(bydim)){
        value <- vector2array(value, dim=indexdim, dimorder=bydim)
      }
      if (nreturn!=nvalue){
        if (nreturn%%nvalue)
          warning("number of elements to replace is not multiple of values for replacement")
        if (ixre){
          # if possible we try to recycle on the C-side, but in case of ix or re we need to recycle already here
          value <- rep(value, length=nreturn)
          nvalue <- nreturn
        }
      }
      if (simple){
        if (!is.null(index$ix))
          value <- value[index$ix]
        else if (index$re)
          value <- rev(value)
      }else{
        if (ixre){
          value <- array(value, dim=indexdim)
          ix <- lapply(index, function(i){
                  if (is.null(i$ix)){
                    if (i$re)
                      rev(1:length(i))
                    else
                      1:length(i)
                  }else{
                    i$ix
                  }
                })
          value <- do.call("[", c(list(value), ix))
        }else{
          ix <- NULL
        }
      }
    }else if(nvalue<1){
      stop("no value for replacement")
    }

    if (simple){
      if (add)
        attr(x, "physical") <- .Call("addset_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
      else
        attr(x, "physical") <- .Call("set_vector", .ffmode[vmode], attr(x, "physical"), index, nreturn, value, PACKAGE="ff")
    }else{
      if (!is.null(vw))
        ffdim <- as.integer(colSums(vw))
      cumffdim <- c(1,cumprod(ffdim[ffdimord])[-ndim])
      cumindexdim <- cumprod(indexdim)
      nreturn <- cumindexdim[ndim]
      cumindexdim <- c(1,cumindexdim[-ndim])[ffdimord]

      if (add)
        attr(x, "physical") <- .Call("addset_array", .ffmode[vmode], attr(x, "physical"), index[ffdimord], as.integer(indexdim[ffdimord]), as.integer(ffdim[ffdimord]), as.integer(ndim), as.integer(nreturn), as.integer(cumindexdim), as.integer(cumffdim), value, PACKAGE="ff")
      else
        attr(x, "physical") <- .Call("set_array", .ffmode[vmode], attr(x, "physical"), index[ffdimord], as.integer(indexdim[ffdimord]), as.integer(ffdim[ffdimord]), as.integer(ndim), as.integer(nreturn), as.integer(cumindexdim), as.integer(cumffdim), value, PACKAGE="ff")
    }
  }

  return(x)
}

#! \name{add}
#! \alias{add}
#! \alias{add.ff}
#! \alias{add.default}
#! \title{ Incrementing an ff or ram object }
#! \description{
#!   Yet another assignment interface in order to allow to formulate \code{x[index,...,add=TRUE]<-value}
#!   in a way which works transparently, not only for ff, but also for ram objects: \code{add(x, value, index, ...)}.
#! }
#! \usage{
#! add(x, \dots)
#! \method{add}{ff}(x, value, \dots)
#! \method{add}{default}(x, value, \dots)
#! }
#! \arguments{
#!   \item{x}{ an ff or ram object }
#!   \item{value}{ the amount to increment, possibly recylcled }
#!   \item{\dots}{ further arguments -- especially index information -- passed to \code{\link{[<-}} or \code{\link{[<-.ff}} }
#! }
#! \value{
#!   invisible()
#! }
#! \author{ Jens Oehlschlägel }
#! \note{
#!   Note that \code{add.default} changes the object in its parent frame and thus violates R's usual functional programming logic.
#!   Duplicated index positions should be avoided, because ff and ram objects behave differently:
#!   \preformatted{
#!   add.ff(x, 1, c(3,3))
#!   # will increment x at position 3 TWICE by 1, while
#!   add.default(x, 1, c(3,3))
#!   # will increment x at position 3 just ONCE by 1
#!   }
#! }
#! \seealso{ \code{\link{swap}}, \code{\link{[.ff}}, \code{\link{LimWarn}} }
#! \examples{
#!    message("incrementing parts of a vector")
#!    x <- ff(0, length=12)
#!    y <- rep(0, 12)
#!    add(x, 1, 1:6)
#!    add(y, 1, 1:6)
#!    x
#!    y
#!
#!    message("incrementing parts of a matrix")
#!    x <- ff(0, dim=3:4)
#!    y <- array(0, dim=3:4)
#!    add(x, 1, 1:2, 1:2)
#!    add(y, 1, 1:2, 1:2)
#!    x
#!    y
#!
#!    message("BEWARE that ff and ram methods differ in treatment of duplicated index positions")
#!    add(x, 1, c(3,3))
#!    add(y, 1, c(3,3))
#!    x
#!    y
#!
#!    rm(x); gc()
#! }
#! \keyword{ IO }
#! \keyword{ data }

# for compatibility with the comatibility generic 'add' which unifies += (almost) for ff and ram
add.ff <- function(
  x
, value
, ... # dummy to keep R CMD check quiet
){
  cl <- match.call(expand.dots=TRUE)
  valueterm <- cl$value
  cl[[1]] <- as.symbol("[")
  cl$value <- NULL
  cl$add <- TRUE
  assigncall <- call("<-", cl, valueterm)
  eval.parent(assigncall)
  invisible()
}
add.default <- function(
  x
, value
, ... # dummy to keep R CMD check quiet
){
  cl <- match.call(expand.dots=TRUE)
  valueterm <- cl$value
  cl[[1]] <- as.symbol("[")
  names(cl)[2] <- ""
  cl$value <- NULL
  cl$add <- NULL
  assigncall <- call("<-", cl, call("+", cl, valueterm))
  eval.parent(assigncall)
  invisible()
}



# Attention: changes its argument x
# add=FALSE writes value and returns oldvalue (BEFORE overwriting)
# add=TRUE  writes and returns oldvalue+value (AFTER overwriting)
swap.default <- function(
  x
, value
, ... # dummy to keep R CMD check quiet
, add=FALSE
){
  cl <- match.call(expand.dots=TRUE)
  valueterm <- cl$value
  cl[[1]] <- as.symbol("[")
  names(cl)[2] <- ""
  cl$value <- NULL
  cl$add <- NULL
  if (add){
    assigncall <- call("<-", cl, call("+", cl, valueterm))
    eval.parent(assigncall)
    ret <- eval.parent(cl)
  }else{
    assigncall <- call("<-", cl, value=valueterm)
    ret <- eval.parent(cl)
    eval.parent(assigncall)
  }
  ret
}


#! \name{LimWarn}
#! \alias{LimWarn}
#! \title{ ff Limitations and Warnings }
#! \description{
#!   This help page lists the currently known limitations of package ff,
#!   as well as differences between ff and ram methods.
#! }
#! \section{Automatic file removal}{
#!   Remind that not giving parameter \code{ff(filename=)} will result in a temporary file in \code{fftempdir} with 'delete' finalizer,
#!   while giving parameter \code{ff(filename=)} will result in a permanent file with 'close' finalizer.
#!   Do avoid setting \code{setwd(getOption("fftempdir"))}!
#!   Make sure you really understand the implications of automatic unlinking of getOption("fftempdir") \code{\link{.onUnload}},
#!   of finalizer choice and of finalizing behaviour at the end of R sessions as defaulted in getOption("fffinonexit").
#!   \bold{Otherwise you might experience 'unexpected' losses of files and data.}
#! }
#! \section{Size of objects}{
#!   Currently ff objects cannot have length zero and are limited to \code{.Machine$integer.max} elements. We have not yet ported the R code to support 64bit double indices (in essence 52 bits integer) although the C++ back-end has been prepared for this.
#!   Furthermore filesize limitations of the OS apply, see \code{\link{ff}}.
#! }
#! \section{Side effects}{
#!   In contrast to standard R expressions, ff expressions violate the functional programming logic and are called for their side effects.
#!   This is also true for ram compatibility functions \code{\link{swap.default}}, and \code{\link{add.default}}.
#! }
#! \section{Hybrid copying semantics}{
#!   If you modify a copy of an ff object, changes of data (\code{\link[ff:[<-.ff]{[<-}}) and of \code{\link[=physical.ff]{physical}} attributes
#!   will be shared, but changes in \code{\link[=physical.ff]{virtual}} and class attributes will not.
#! }
#! \section{Limits of compatibility between ff and ram objects}{
#!   If it's not too big, you can move an ff object completely into R's RAM through \code{\link{as.ram}}.
#!   However, you should watch out for three limitations:
#!   \enumerate{
#!     \item Ram objects don't have hybrid copying semantics; changes to a copy of a ram object will never change the original ram object
#!     \item Assigning values to a ram object can easily upgrade to a higher \code{\link{storage.mode}}. This will create conflicts with the
#!           \code{\link{vmode}} of the ram object, which goes undetected until you try to write back to disk through \code{\link{as.ff}}.
#!     \item Writing back to disk with \code{\link{as.ff}} under the same filename requires that the original ff object has been deleted
#!           (or at least closed if you specify parameter \code{overwrite=TRUE}).
#!   }
#! }
#! \section{Index expressions}{
#!   ff index expressions do not allow zeros and NAs, see see \code{\link{[.ff}} and see \code{\link{as.hi}}
#! }
#! \section{Availablility of bydim parameter}{
#!   Parameter \code{bydim} is only available in ff access methods, see \code{\link{[.ff}}
#! }
#! \section{Availablility of add parameter}{
#!   Parameter \code{add} is only available in ff access methods, see \code{\link{[.ff}}
#! }
#! \section{Compatibility of swap and add}{
#!   If index expressions contain duplicated positions, the ff and ram methods for \code{\link{swap}}
#!   and \code{\link{add}} will behave differently, see \code{\link{swap}}.
#! }
#! \section{Definition of [[ and [[<-}{
#!   You should consider the behaviour of \code{\link{[[.ff}} and \code{\link{[[<-.ff}} as undefined and not use them in programming.
#!   Currently they are shortcuts to \code{\link{get.ff}} and \code{\link{set.ff}}, which unlike \code{\link{[.ff}} and \code{\link{[<-.ff}}
#!   do not support \code{\link{factor}} and \code{\link{POSIXct}}, nor \code{\link{dimorder}} or virtual windows \code{\link{vw}}.
#!   In contrast to the standard methods, \code{\link{[[.ff}} and \code{\link{[[<-.ff}} only accepts positive integer index positions.
#!   The definition of \code{\link{[[.ff}} and \code{\link{[[<-.ff}} may be changed in the future.
#! }
#! \section{Multiple vector interpretation in arrays }{
#!   R objects have always standard \code{\link{dimorder} 1:length(dim)}.
#!   In case of non-standard dimorder (see \code{\link{dimorderStandard}})
#!   the vector sequence of array elements in R and in the ff file differs.
#!   To access array elements in file order, you can use \code{\link{getset.ff}}, \code{\link{readwrite.ff}}
#!   or copy the ff object and set \code{dim(ff)<-NULL} to get a vector view into the ff object
#!   (using \code{[} dispatches the vector method \code{\link{[.ff}}).
#!   To access the array elements in R standard dimorder you simply use \code{[} which dispatches
#!   to \code{\link{[.ff_array}}. Note that in this case \code{\link{as.hi}} will unpack the complete index, see next section.
#! }
#! \section{RAM expansion of index expressions}{
#!   Some index expressions do not consume RAM due to the \code{\link{hi}} representation.
#!   For example \code{1:n} will almost consume no RAM however large n.
#!   However, some index expressions are expanded and require to \code{\link{maxindex}(i) * .rambytes["integer"]} bytes,
#!   either because the sorted sequence of index positions cannot be rle-packed efficiently
#!   or because \code{\link{hiparse}} cannot yet parse such expression and falls back to evaluating/expanding the index expression.
#!   If the index positions are not sorted, the index will be expanded and a second vector is needed to store the information for re-ordering,
#!   thus the index requires \code{2 * \link{maxindex}(i) * .rambytes["integer"]} bytes.
#! }
#! \section{RAM expansion when recycling assigment values}{
#!   Some assignment expressions do not consume RAM for recycling. For example \code{x[1:n] <- 1:k}
#!   will not consume RAM however large is n compared to k, when x has standard \code{\link{dimorder}}.
#!   However, if \code{length(value)>1}, assignment expressions with non-ascending index positions trigger recycling the value R-side to the full index length.
#!   This will happen if \code{\link{dimorder}} does not match parameter \code{bydim} or if the index is not sorted in ascending order.
#! }
#! \section{Byteorder imcompatibility}{
#!   Note that ff files cannot been transferred between systems with different byteorder.
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ package }

