# reading/writing csv files into/from ffdf objects
# (c) 2009 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2009-09-09
# Last changed: 2009-09-27

# source("d:/mwp/eanalysis/ff/R/ffcsv.R")


# 3 columns with 78 mio records (integer, factor, factor)
# read.table: read.csv 5 min + write.ffdf 1 min
# clone.ffdf: read+write: 1 min
# write.table: read.ffdf 0.8 min + write.csv 9.5 min

# 7 columns with 78 mio records ("boolean"    "byte"  "single"   "short"   "short"  "single"  "single")
# create ffdf 102.7 sec
# write.table: ffdf-read=100.45sec  csv-write=2553.57sec  TOTAL=2654.02sec
# read.table: csv-read=2769.26sec  ffdf-write=115.63sec  TOTAL=2884.89sec
# clone.ffdf: read+write: 175 min


#readratio  2769.26/100.45 = 27.5
#writeratio 2553.57/115.63 = 22

# ff read 1.3s / 1mio / 7 col = 0.2s / mio / col
# ff write 1.5s / 1mio / 7 col = 0.21s / mio / col
# csv write 35.5s / 1mio / 7 col = 5 sec / 1mio / col
# csv read 32.7s / 1mio / 7 col = 4.7 sec / 1mio / col



# NOTE that this colClass implementation works only because accidentally the last position of the oldClasses is needed
# for c("ordered","factor") read.table does not want "ordered"
# for c("POSIXt","POSIXct") "POSIXct" is in the last position only due to "historical error" (Gabor Grothendieck, r-help, 26.9.2009)

colClass <- function(x)
UseMethod("colClass")

if(getRversion() < "2.12.0"){
  colClass.default <- function(x){
    cl <- class(x)
    cl[length(cl)]
  }
}else{
  colClass.default <- function(x){
    cl <- class(x)
    if(inherits(x, "POSIXct")) "POSIXct" else cl[length(cl)]
  }
}


colClass.ff <- function(x){
  if (length(x))
    x <- x[1]
  else
    x <- x[]
  NextMethod()
}



#! \name{read.table.ffdf}
#! \Rdversion{1.1}
#! \alias{read.table.ffdf}
#! \alias{read.csv.ffdf}
#! \alias{read.csv2.ffdf}
#! \alias{read.delim.ffdf}
#! \alias{read.delim2.ffdf}
#! \title{
#!   Importing csv files into ff data.frames
#! }
#! \description{
#!   Function \code{read.table.ffdf} reads separated flat files into \code{\link{ffdf}} objects, very much like (and using) \code{\link{read.table}}.
#!   It can also work with any convenience wrappers like \code{\link{read.csv}} and provides its own convenience wrapper (e.g. \code{read.csv.ffdf}) for R's usual wrappers.
#! }
#! \usage{
#! read.table.ffdf(
#!   x = NULL
#! , file, fileEncoding = ""
#! , nrows = -1, first.rows = NULL, next.rows = NULL
#! , levels = NULL, appendLevels = TRUE
#! , FUN = "read.table", ...
#! , transFUN = NULL
#! , asffdf_args = list()
#! , BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE
#! )
#! read.csv.ffdf(...)
#! read.csv2.ffdf(...)
#! read.delim.ffdf(...)
#! read.delim2.ffdf(...)
#! }
#! \arguments{
#!   \item{x}{
#! NULL or an optional \code{\link{ffdf}} object to which the read records are appended.
#! If this is provided, it defines crucial features that are otherwise determnined during the 'first' chunk of reading:
#! \code{\link[=vmode.ffdf]{vmode}s}, \code{\link[=dimnames.ffdf]{colnames}}, \code{colClasses}, sequence of predefined \code{\link[=levels.ff]{levels}}.
#! In order to also read the first chunk into such predefined ffdf, an \code{x} with 1 row is treated special: instead of appending the first row will be overwritten.
#! This is necessary because we cannot provide \code{x} with zero rows (we cannot create ff vectors with zero elements).
#! }
#!   \item{file}{
#!     the name of the file which the data are to be read from.
#!     Each row of the table appears as one line of the file.  If it does
#!     not contain an \emph{absolute} path, the file name is
#!     \emph{relative} to the current working directory,
#!     \code{\link{getwd}()}. Tilde-expansion is performed where supported.
#!
#!     Alternatively, \code{file} can be a readable text-mode
#!     \code{\link{connection}} (which will be opened for reading if
#!     necessary, and if so \code{\link{close}}d (and hence destroyed) at
#!     the end of the function call).
#! }
#!   \item{fileEncoding}{
#!     character string: if non-empty declares the
#!     encoding used on a file (not a connection) so the character data can
#!     be re-encoded.  See \code{\link{file}}.
#! }
#!   \item{nrows}{
#!   integer: the maximum number of rows to read in (includes first.rows in case a 'first' chunk is read)
#!   Negative and other invalid values are ignored.
#! }
#!   \item{first.rows}{
#!   integer: number of rows to be read in the first chunk, see details. Default is the value given at \code{next.rows} or \code{1e3} otherwise.
#!   Ignored if \code{x} is given.
#! }
#!   \item{next.rows}{
#!   integer: number of rows to be read in further chunks, see details.
#!   By default calculated as \code{BATCHBYTES \%/\% sum(.rambytes[\link{vmode}(x)])}
#! }
#!   \item{levels}{
#!   NULL or an optional list, each element named with col.names of factor columns specifies the \code{\link[=levels.ff]{levels}}
#!   Ignored if \code{x} is given.
#! }
#!   \item{appendLevels}{
#!   logical.
#!   A vector of permissions to expand \code{\link[=levels.ff]{levels}} for factor columns.
#!   Recycled as necessary, or if the logical vector is named, unspecified values are taken to be \code{TRUE}.
#!   Ignored during processing of the 'first' chunk
#! }
#!   \item{FUN}{
#!   character: name of a function that is called for reading each chunk, see \code{\link{read.table}}, \code{\link{read.csv}}, etc.
#! }
#!   \item{\dots}{
#!   further arguments, passed to \code{FUN} in \code{read.table.ffdf}, or passed to \code{read.table.ffdf} in the convenience wrappers
#! }
#!   \item{transFUN}{
#!   NULL or a function that is called on each data.frame chunk after reading with \code{FUN} and before further processing (for filtering, transformations etc.)
#! }
#!   \item{asffdf_args}{
#!   further arguments passed to \code{\link{as.ffdf}} when converting the \code{\link{data.frame}} of the first chunk to \code{\link{ffdf}}.
#!   Ignored if \code{x} is given.
#! }
#!   \item{BATCHBYTES}{
#!   integer: bytes allowed for the size of the \code{\link{data.frame}} storing the result of reading one chunk. Default \code{getOption("ffbatchbytes")}.
#! }
#!   \item{VERBOSE}{
#!   logical: TRUE to verbose timings for each processed chunk (default FALSE)
#! }
#! }
#! \details{
#!   \code{read.table.ffdf} has been designed to read very large (many rows) separated flatfiles in row-chunks
#!   and store the result in a \code{\link{ffdf}} object on disk, but quickly accessible via \code{\link{ff}} techniques.
#!   \cr
#!   The first chunk is read with a default of 1000 rows, for subsequent chunks the number of rows is calculated to not require more RAM than \code{getOption("ffbatchbytes")}.
#!   The following could be indications to change the parameter \code{first.rows}:
#!   \enumerate{
#!     \item set \code{first.rows=-1} to read the complete file in one go (requires enough RAM)
#!     \item set \code{first.rows} to a smaller number if the pre-allocation of RAM for the first chunk with parameter \code{nrows} in \code{\link{read.table}} is too large, i.e. with many columns on machine with little RAM.
#!     \item set \code{first.rows} to a larger number if you expect better factor level ordering (factor levels are sorted in the first chunk, but not at subsequent chunks, however, factor level ordering can be fixed later, see below).
#!   }
#!   By default the \code{\link{ffdf}} object is created on the fly at the end of reading the 'first' chunk, see argument \code{first.rows}.
#!   The creation of the \code{\link{ffdf}} object is done via \code{\link{as.ffdf}} and can be finetuned by passing argument \code{asffdf_args}.
#!   Even more control is possible by passing in a \code{\link{ffdf}} object as argument \code{x} to which the read records are appended.
#!   \cr
#!   \code{read.table.ffdf} has been designed to behave as much like \code{\link{read.table}} as possible. Hoever, note the following differences:
#!   \enumerate{
#!     \item Arguments 'colClasses' and 'col.names' are now enforced also during 'next.rows' chunks.
#!           For example giving \code{colClasses=NA} will force that no colClasses are derived from the \code{first.rows} respective from the \code{\link{ffdf}} object in parameter \code{x}.
#!     \item colClass 'ordered' is allowed and will create an \code{\link{ordered}} factor
#!     \item character vector are not supported, character data must be read as one of the following colClasses: 'Date', 'POSIXct', 'factor, 'ordered'.
#!           By default character columns are read as factors.
#!           Accordingly arguments 'as.is' and 'stringsAsFactors' are not allowed.
#!     \item the sequence of \code{\link{levels.ff}} from chunked reading can depend on chunk size: by default new levels found on a chunk are appended to the levels found in previous chunks, no attempt is made to sort and recode the levels during chunked processing, levels can be sorted and recoded most efficiently \emph{after} all records have been read using \code{\link{sortLevels}}.
#!     \item the default for argument 'comment.char' is \code{""} even for those FUN that have a different default. However, explicit specification of 'comment.char' will have priority.
#!   }
#! }
#! \note{
#!     Note that using the 'skip' argument still requires to read the file from beginning in order to count the lines to be skipped.
#!     If you first read part of the file in order to understand its structure and then want to continue,
#!     a more efficient solution that using 'skip' is opening a \code{\link{file}} \code{\link{connection}} and pass that to argument 'file'.
#!     \code{read.table.ffdf} does the same in order to skip efficiently over previously read chunks.
#! }
#! \value{
#!   An \code{\link{ffdf}} object. If created during the 'first' chunk pass, it will have one \code{\link[=physical.ffdf]{physical}} component per \code{\link[=virtual.ffdf]{virtual}} column.
#! }
#! \author{
#!   Jens Oehlschlägel, Christophe Dutang
#! }
#! \seealso{
#!   \code{\link{write.table.ffdf}}, \code{\link{read.table}}, \code{\link{ffdf}}
#! }
#! \examples{
#!     x <- data.frame(log=rep(c(FALSE, TRUE), length.out=26), int=1:26, dbl=1:26 + 0.1
#! , fac=factor(letters), ord=ordered(LETTERS)
#! , dct=Sys.time()+1:26, dat=seq(as.Date("1910/1/1"), length.out=26, by=1))
#!     x <- x[c(13:1, 13:1),]
#!     csvfile <- tempPathFile(path=getOption("fftempdir"), extension="csv")
#!     write.csv(x, file=csvfile, row.names=FALSE)
#!
#!     y <- read.csv(file=csvfile, header=TRUE)
#!     y
#!     cat("Read csv with header\n")
#!     ffx <- read.csv.ffdf(file=csvfile, header=TRUE)
#!     ffx
#!     sapply(ffx[,], class)
#!
#!     message("NOTE that read.table fails for ordered factors, this is fixed in read.table.ffdf")
#!     try(read.csv(file=csvfile, header=TRUE, colClasses=c(ord="ordered")))
#!     # TODO could fix this with the following two commands (Gabor Grothendieck) 
#!     # but does not know what bad side-effects this could have
#!      #setOldClass("ordered")
#!      #setAs("character", "ordered", function(from) ordered(from))
#!     y <- read.csv(file=csvfile, header=TRUE, colClasses=c(dct="POSIXct", dat="Date"))
#!     y
#!
#!     ffx <- read.csv.ffdf(file=csvfile, header=TRUE
#! , colClasses=c(ord="ordered", dct="POSIXct", dat="Date"))
#!     ffx
#!     sapply(ffx[,], class)
#!
#!     message("NOTE that reading in chunks can change the sequence of levels and thus the coding")
#!     message("(Sorting levels during chunked reading can be too expensive)")
#!     ffx <- read.csv.ffdf(file=csvfile, header=TRUE
#! , colClasses=c(ord="ordered", dct="POSIXct", dat="Date"), first.rows=6, next.rows=10
#! , VERBOSE=TRUE)
#!     y <- ffx$fac[]
#!     print(levels(y))
#!     data.frame(values=as.character(y), codes=as.integer(y))
#!
#!     message("If we don't know the levels we can sort then after reading")
#!     message("(Will rewrite all factor codes)")
#!     message("NOTE that you MUST assign the return value of sortLevels()")
#!     ffx <- sortLevels(ffx)
#!     y <- ffx$fac[]
#!     print(levels(y))
#!     data.frame(values=as.character(y), codes=as.integer(y))
#!
#!     message("If we KNOW the levels we can fix levels upfront")
#!     ffx <- read.csv.ffdf(file=csvfile, header=TRUE
#! , colClasses=c(ord="ordered", dct="POSIXct", dat="Date")
#! , first.rows=6, next.rows=10, levels=list(fac=letters, ord=LETTERS))
#!     y <- ffx$fac[]
#!     print(levels(y))
#!     data.frame(values=as.character(y), codes=as.integer(y))
#!
#!     message("Or we inspect a sufficiently large chunk of data and use those")
#!     ffx1 <- read.csv.ffdf(file=csvfile, header=TRUE
#! , colClasses=c(ord="ordered", dct="POSIXct", dat="Date"), nrows=13)
#!     ffx <- read.csv.ffdf(x=ffx1, file=csvfile, header=FALSE, skip=1+nrow(ffx1), VERBOSE=TRUE)
#!
#!     message("We can check for unexpected factor levels, say we only allowed a:l")
#!     ffx <- read.csv.ffdf(file=csvfile, header=TRUE
#! , colClasses=c(ord="ordered", dct="POSIXct", dat="Date")
#! , levels=list(fac=letters[1:12], ord=LETTERS[1:12]), appendLevels=FALSE)
#!     sapply(colnames(ffx), function(i)sum(is.na(ffx[[i]][])))
#!
#!     message("We can fine-tune the creation of the ffdf:")
#!     message("- let's create the ff files outside of fftempdir")
#!     message("- let's reduce required disk space and thus file.system cache RAM")
#!     vmode(ffx)
#!     message("By default we had record size")
#!     sum(.ffbytes[vmode(ffx)])
#!
#!     ffy <- read.csv.ffdf(file=csvfile, header=TRUE
#!     , colClasses=c(ord="ordered", dct="POSIXct", dat="Date")
#!     , asffdf_args=list(
#!         vmode = c(log="boolean", int="byte", dbl="single", fac="nibble"
#!                 , ord="nibble", dct="single", dat="single")
#!       , col_args=list(pattern = "./csv")  # create in getwd() with prefix csv
#!       )
#!     )
#!     vmode(ffy)
#!     message("This recordsize is more than 50\% reduced")
#!     sum(.ffbytes[vmode(ffy)])
#!
#!     message("Don't forget to wrap-up files that are not in fftempdir")
#!     delete(ffy); rm(ffy)
#!     message("It's a good habit to also wrap-up temporary stuff (or at least know how this is done)")
#!     rm(ffx); gc()
#!
#!     fwffile <- tempfile()
#!
#!     cat(file=fwffile, "123456", "987654", sep="\n")
#!     x <- read.fwf(fwffile, widths=c(1,2,3))    #> 1 23 456 \ 9 87 654
#!     y <- read.table.ffdf(file=fwffile, FUN="read.fwf", widths=c(1,2,3))
#!     stopifnot(identical(x, y[,]))
#!     x <- read.fwf(fwffile, widths=c(1,-2,3))   #> 1 456 \ 9 654
#!     y <- read.table.ffdf(file=fwffile, FUN="read.fwf", widths=c(1,-2,3))
#!     stopifnot(identical(x, y[,]))
#!     unlink(fwffile)
#!
#!     cat(file=fwffile, "123", "987654", sep="\n")
#!     x <- read.fwf(fwffile, widths=c(1,0, 2,3))    #> 1 NA 23 NA \ 9 NA 87 654
#!     y <- read.table.ffdf(file=fwffile, FUN="read.fwf", widths=c(1,0, 2,3))
#!     stopifnot(identical(x, y[,]))
#!     unlink(fwffile)
#!
#!     cat(file=fwffile, "123456", "987654", sep="\n")
#!     x <- read.fwf(fwffile, widths=list(c(1,0, 2,3), c(2,2,2))) #> 1 NA 23 456 98 76 54
#!     y <- read.table.ffdf(file=fwffile, FUN="read.fwf", widths=list(c(1,0, 2,3), c(2,2,2)))
#!     stopifnot(identical(x, y[,]))
#!     unlink(fwffile)
#!
#!     \dontshow{
#!        x <- read.csv(file=csvfile, header=TRUE)
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE)
#!        stopifnot(identical(x, y[,]))
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE, nrows=13)
#!        stopifnot(identical(x[1:13,], y[,]))
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE, first.rows=12)
#!        y <- sortLevels(y)
#!        stopifnot(identical(x, y[,]))
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE, nrows=13, first.rows=12)
#!        y <- sortLevels(y)
#!        stopifnot(identical(x[1:13,], y[,]))
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE, nrows=12, first.rows=12)
#!        y <- sortLevels(y)
#!        stopifnot(!identical(x[1:12,], y[,]))
#!        stopifnot(identical(as.character(as.matrix(x[1:12,])), as.character(as.matrix(y[,]))))
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE, nrows=11, first.rows=12)
#!        y <- sortLevels(y)
#!        stopifnot(!identical(x[1:11,], y[,]))
#!        stopifnot(identical(as.character(as.matrix(x[1:11,])), as.character(as.matrix(y[,]))))
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE, first.rows=-1)
#!        stopifnot(identical(x, y[,]))
#!
#!        y <- read.csv.ffdf(file=csvfile, header=TRUE, nrows=13, first.rows=12, appendLevels=c(ord=FALSE))
#!        stopifnot(is.na(y$ord[13]) && !is.na(y$fac[13]))
#!     }
#!
#!     unlink(csvfile)
#! }
#! \keyword{IO}
#! \keyword{file}
#! \keyword{connection}


read.table.ffdf <- function(
  x = NULL
, file
, fileEncoding = ""
, nrows       = -1
, first.rows  = NULL
, next.rows   = NULL
, levels      = NULL
, appendLevels = TRUE  # levels are expanded (but not sorted)
, FUN = "read.table"
, ...    # further arguments pmatched against args(read.table)
, transFUN = NULL
, asffdf_args = list()    # further arguments passed to as.ffdf() (ignored if 'x' gives an ffdf object )
, BATCHBYTES = getOption("ffbatchbytes")
, VERBOSE = FALSE
){

  append <- !is.null(x)
  if (append && !inherits(x, "ffdf"))
    stop("only ffdf objects can be used for appending (and skipping the first.row chunk)")

  if (VERBOSE){
    read.time <- 0
    write.time <- 0
  }
  rt.args <- list(...)
  lnam <- names(rt.args)
  if (FUN=="read.fwf")
    fnam <- names(as.list(args(read.fwf)))
  else
    fnam <- names(as.list(args(read.table)))
  i <- pmatch(lnam, fnam)
  if (any(is.na(i)))
    stop("unkown arguments: ", paste(lnam[is.na(i)], sep="", collapse=","))
  names(rt.args) <- fnam[i]

  if (!is.null(rt.args$as.is))
    stop("argument 'as.is' not allowed, all columns must be read as factors or ordered factors")
  if (!is.null(rt.args$stringsAsFactors))
    stop("argument 'stringsAsFactors' not allowed, all columns must be read as factors or ordered factors")


  # { begin 1:1 from read.table
  if (is.character(file)) {
      file <- if (nzchar(fileEncoding))
          file(file, "r", encoding = fileEncoding)
      else file(file, "r")
      on.exit(close(file))
  }
  if (!inherits(file, "connection"))
      stop("'file' must be a character string or connection")
  if (!isOpen(file, "r")) {
      open(file, "r")
      on.exit(close(file))
  }
  # } end 1:1 from read.table

  rt.args$file <- file
  rt.args$fileEncoding <- fileEncoding

  if (is.null(rt.args$comment.char))
    rt.args$comment.char <- ""  # faster than read.table default "#"

  nrows <- as.integer(nrows)
  N <- 0L

  no.colClasses <- is.na(match("colClasses", names(rt.args)))
  no.col.names <- is.na(match("col.names", names(rt.args)))

  if (!append){

    if (VERBOSE){
      cat("read.table.ffdf ", N+1L,"..", sep="")
      read.start <- proc.time()[3]
    }


    if (no.colClasses)
      rt.args$colClasses <- as.character(NA)
    # xx fix a failure of read.table to properly handle "ordered"
    colClasses <- rt.args$colClasses
    rt.args$colClasses[!is.na(colClasses) & colClasses=="ordered"] <- "factor"

    if (is.null(first.rows)){
      if (is.null(next.rows))
        first.rows <- 1e3L
      else
        first.rows <- as.integer(next.rows)
    }else
      first.rows <- as.integer(first.rows)
    if (nrows>=0L && nrows<first.rows)
      first.rows <- nrows

    if (first.rows==1)
      stop("first.rows must not be 1")

    rt.args$nrows <- first.rows

    dat <- do.call(FUN, rt.args)

    # do this already here
    if (no.col.names){
      rt.args$col.names <- colnames(dat)
      no.col.names <- FALSE
    }
    if (no.colClasses){
      rt.args$colClasses <- sapply(seq_len(ncol(dat)), function(i)colClass(dat[[i]]))
      no.colClasses <- FALSE
    }

    n.orig <- nrow(dat)
    if (!is.null(transFUN))
      dat <- transFUN(dat)

    need.next <- n.orig==first.rows
    n <- nrow(dat)
    N <- n

    if (!is.null(levels)){
      cnam <- colnames(dat)
      lnam <- names(levels)
      if (!is.list(levels) || is.null(lnam) || any(is.na(match(lnam,cnam))))
        stop("levels must be a list with names matching column names of the first data.frame read")
      for (i in lnam){
        dat[[i]] <- recodeLevels(dat[[i]], levels[[i]])
      }
    }

    if (VERBOSE){
      write.start <- proc.time()[3]
      read.time <- read.time + (write.start - read.start)
      cat(N, " (", n, ")  csv-read=", round(write.start-read.start, 3), "sec", sep="")
    }

    x <- do.call("as.ffdf", c(list(dat), asffdf_args))

    # now fix ordered
    colClasses <- repnam(colClasses, colnames(x), default=NA)
    i.fix <- seq_len(ncol(dat))[!is.na(match(colClasses, "ordered"))]
    for (i in i.fix)
      virtual(x[[i]])$ramclass <- c("ordered","factor")

    if (VERBOSE){
      write.stop <- proc.time()[3]
      write.time <- write.time + (write.stop - write.start)
      cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
    }

    rt.args$skip <- 0L
    rt.args$header <- FALSE
  }

  if (append || need.next){

    k <- ncol(x)
    if (no.col.names){
      rt.args$col.names <- colnames(x)
    }
    if (no.colClasses){
      rt.args$colClasses <- sapply(seq_len(ncol(x)), function(i)colClass(x[[i]]))
    }

    if (is.null(next.rows)){
      recordsize <- sum(.rambytes[vmode(x)])
      next.rows <- BATCHBYTES %/% recordsize
      if (next.rows<1L){
        next.rows <- 1L
        warning("single record does not fit into BATCHBYTES")
      }
    }else{
      next.rows <- as.integer(next.rows)
    }
    rt.args$nrows <- next.rows

    appendLevels <- repnam(appendLevels, colnames(x), default=TRUE)
    if(any(appendLevels)){
      i.fac <- seq_len(k)
      i.fac <- i.fac[appendLevels & sapply(i.fac, function(i)is.factor(x[[i]]))]
    }

    while(TRUE){
      if (nrows>=0L && N+next.rows > nrows)
        rt.args$nrows <- nrows - N
      if (rt.args$nrows<1L){
        break
      }

      if (VERBOSE){
        cat("read.table.ffdf ", N+1L,"..", sep="")
        read.start <- proc.time()[3]
      }

      dat <- do.call(FUN, rt.args)
      n.orig <- nrow(dat)
      if (!is.null(transFUN))
        dat <- transFUN(dat)

      n <- nrow(dat)
      N <- N + n
      if (VERBOSE){
        write.start <- proc.time()[3]
        read.time <- read.time + (write.start - read.start)
        cat(N, " (", n, ")  csv-read=", round(write.start-read.start, 3), "sec", sep="")
      }
      if (n<1L){
        if (VERBOSE)
          cat("\n")
        break
      }

      if(any(appendLevels))
      for (i in i.fac){
        lev <- unique(c(levels(x[[i]]),levels(dat[[i]])))  # we save a call to the more general appendLevels() here
        levels(x[[i]]) <- lev
        dat[[i]] <- recodeLevels(dat[[i]], lev)
      }

      nff <- nrow(x)
      if (nff==1){  # since we cannot create a ffdf with zero rows, we allow to use nrow()==1 instead
        nrow(x) <- n
        x[,] <- dat
      }else{
        nrow(x) <- nff + n
        i <- hi(nff+1L, nff+n)
        x[i,] <- dat
      }

      if (VERBOSE){
        write.stop <- proc.time()[3]
        write.time <- write.time + (write.stop - write.start)
        cat(" ffdf-write=", round(write.stop-write.start, 3), "sec\n", sep="")
      }

      if (n.orig<rt.args$nrows){
        break
      }

      rt.args$skip <- 0L
      rt.args$header <- FALSE
    }

  }

  if (VERBOSE){
    cat(" csv-read=", round(read.time, 3), "sec  ffdf-write=", round(write.time, 3), "sec  TOTAL=", round(read.time+write.time, 3), "sec\n", sep="")
  }

  return(x)
}



#! \name{write.table.ffdf}
#! \Rdversion{1.1}
#! \alias{write.table.ffdf}
#! \alias{write.csv.ffdf}
#! \alias{write.csv2.ffdf}
#! \alias{write.csv}
#! \alias{write.csv2}
#! \title{
#!   Exporting csv files from ff data.frames
#! }
#! \description{
#!   Function \code{write.table.ffdf} writes a \code{\link{ffdf}} object to a separated flat file, very much like (and using) \code{\link{write.table}}.
#!   It can also work with any convenience wrappers like \code{\link{write.csv}} and provides its own convenience wrapper (e.g. \code{write.csv.ffdf}) for R's usual wrappers.
#! }
#! \usage{
#! write.table.ffdf(x = NULL
#! , file, append = FALSE
#! , nrows = -1, first.rows = NULL, next.rows = NULL
#! , FUN = "write.table", ...
#! , transFUN = NULL
#! , BATCHBYTES = getOption("ffbatchbytes")
#! , VERBOSE = FALSE
#! )
#! write.csv.ffdf(...)
#! write.csv2.ffdf(...)
#! write.csv(...)
#! write.csv2(...)
#! }
#! \arguments{
#!   \item{x}{
#!   a \code{\link{ffdf}} object which to export to the separated file
#! }
#!   \item{file}{
#!   either a character string naming a file or a connection
#!   open for writing.  \code{""} indicates output to the console.
#! }
#!   \item{append}{
#!   logical. Only relevant if \code{file} is a character
#!   string.  If \code{TRUE}, the output is appended to the
#!   file.  If \code{FALSE}, any existing file of the name is destroyed.
#! }
#!   \item{nrows}{
#!   integer: the maximum number of rows to write in (includes first.rows in case a 'first' chunk is read)
#!   Negative and other invalid values are ignored.
#! }
#! \item{first.rows}{
#!   the number of rows to write with the first chunk (default: next.rows)
#! }
#!   \item{next.rows}{
#!   integer: number of rows to write in further chunks, see details.
#!   By default calculated as \code{BATCHBYTES \%/\% sum(.rambytes[\link{vmode}(x)])}
#! }
#!   \item{FUN}{
#!   character: name of a function that is called for writing each chunk, see \code{\link{write.table}}, \code{\link{write.csv}}, etc.
#! }
#!   \item{\dots}{
#!   further arguments, passed to \code{FUN} in \code{write.table.ffdf}, or passed to \code{write.table.ffdf} in the convenience wrappers
#! }
#!   \item{transFUN}{
#!   NULL or a function that is called on each data.frame chunk before writing with \code{FUN} (for filtering, transformations etc.)
#! }
#!   \item{BATCHBYTES}{
#!   integer: bytes allowed for the size of the \code{\link{data.frame}} storing the result of reading one chunk. Default \code{getOption("ffbatchbytes")}.
#! }
#!   \item{VERBOSE}{
#!   logical: TRUE to verbose timings for each processed chunk (default FALSE)
#! }
#! }
#! \details{
#!   \code{write.table.ffdf} has been designed to export very large \code{\link{ffdf}} objects to separated flatfiles in chunks.
#!   The first chunk is potentially written with col.names. Further chunks are appended.
#!   \cr
#!   \code{write.table.ffdf} has been designed to behave as much like \code{\link{write.table}} as possible. However, note the following differences:
#!   \enumerate{
#!     \item by default \code{\link[=dimnames.ffdf]{row.names}} are only written if the \code{\link{ffdf}} has row.names.
#!   }
#! }
#! \value{
#!   \code{\link{invisible}}
#! }
#! \note{
#!   \code{\link[utils]{write.csv}} and \code{\link[utils]{write.csv2}} have been fixed in order to suppress \code{col.names} if \code{append=TRUE} is passed.
#!   Note also that \code{write.table.ffdf} passes \code{col.names=FALSE} for all chunks following the first chunk - but not so for \code{FUN="write.csv"} and \code{FUN="write.csv2"} .
#! }
#! \author{
#!   Jens Oehlschlägel, Christophe Dutang
#! }
#! \seealso{
#!   \code{\link{read.table.ffdf}}, \code{\link{write.table}}, \code{\link{ffdf}}
#! }
#! \examples{
#!    x <- data.frame(log=rep(c(FALSE, TRUE), length.out=26), int=1:26, dbl=1:26 + 0.1
#! , fac=factor(letters), ord=ordered(LETTERS), dct=Sys.time()+1:26
#! , dat=seq(as.Date("1910/1/1"), length.out=26, by=1))
#!    ffx <- as.ffdf(x)
#!
#!    csvfile <- tempPathFile(path=getOption("fftempdir"), extension="csv")
#!
#!    write.csv.ffdf(ffx, file=csvfile)
#!    write.csv.ffdf(ffx, file=csvfile, append=TRUE)
#!
#!    ffy <- read.csv.ffdf(file=csvfile, header=TRUE
#! , colClasses=c(ord="ordered", dct="POSIXct", dat="Date"))
#!
#!    rm(ffx, ffy); gc()
#!    unlink(csvfile)
#!
#!  \dontrun{
#!   # Attention, this takes very long
#!   vmodes <- c(log="boolean", int="byte", dbl="single"
#!, fac="short", ord="short", dct="single", dat="single")
#!
#!   message("create a ffdf with 7 columns and 78 mio rows")
#!   system.time({
#!     x <- data.frame(log=rep(c(FALSE, TRUE), length.out=26), int=1:26, dbl=1:26 + 0.1
#! , fac=factor(letters), ord=ordered(LETTERS), dct=Sys.time()+1:26
#! , dat=seq(as.Date("1910/1/1"), length.out=26, by=1))
#!     x <- do.call("rbind", rep(list(x), 10))
#!     x <- do.call("rbind", rep(list(x), 10))
#!     x <- do.call("rbind", rep(list(x), 10))
#!     x <- do.call("rbind", rep(list(x), 10))
#!     ffx <- as.ffdf(x, vmode = vmodes)
#!     for (i in 2:300){
#!       message(i, "\n")
#!       last <- nrow(ffx) + nrow(x)
#!       first <- last - nrow(x) + 1L
#!       nrow(ffx) <- last
#!       ffx[first:last,] <- x
#!     }
#!   })
#!
#!
#!   csvfile <- tempPathFile(path=getOption("fftempdir"), extension="csv")
#!
#!   write.csv.ffdf(ffx, file=csvfile, VERBOSE=TRUE)
#!   ffy <- read.csv.ffdf(file=csvfile, header=TRUE
#! , colClasses=c(ord="ordered", dct="POSIXct", dat="Date")
#! , asffdf_args=list(vmode = vmodes), VERBOSE=TRUE)
#!
#!   rm(ffx, ffy); gc()
#!   unlink(csvfile)
#!  }
#! }
#! \keyword{IO}
#! \keyword{file}
#! \keyword{connection}


write.table.ffdf <- function(
  x = NULL
, file
, append  = FALSE
, nrows       = -1
, first.rows  = NULL
, next.rows   = NULL
, FUN = "write.table"
, ...    # further arguments pmatched against args(write.table)
, transFUN = NULL
, BATCHBYTES = getOption("ffbatchbytes")
, VERBOSE = FALSE
){

  stopifnot(inherits(x, "ffdf"))

  if (nrows<0)
    nrows <- nrow(x)
  else
    nrows <- min(as.integer(nrows), nrow(x))

  if (nrows<1L)
    return(invisible())

  if (is.null(next.rows)){
    recordsize <- sum(.rambytes[vmode(x)])
    next.rows <- BATCHBYTES %/% recordsize
    if (next.rows<1L){
      next.rows <- 1L
      warning("single record does not fit into BATCHBYTES")
    }
  }else{
    next.rows <- as.integer(next.rows)
  }
  if (is.null(first.rows)){
    if (is.null(next.rows))
      first.rows <- 1e3L
    else
      first.rows <- next.rows
  }else
    first.rows <- as.integer(first.rows)
  if (nrows<first.rows)
    first.rows <- nrows

  if (VERBOSE){
    read.time <- 0
    write.time <- 0
  }
  wt.args <- list(...)
  lnam <- names(wt.args)
  fnam <- names(as.list(do.call("args", list(FUN))))
  i <- pmatch(lnam, fnam)
  if (any(is.na(i)))
    stop("unkown arguments: ", paste(lnam[is.na(i)], sep="", collapse=","))
  names(wt.args) <- fnam[i]

  wt.args$append <- append

  if (is.null(wt.args$row.names))
    wt.args$row.names <- !is.null(rownames(x))

  # { begin 1:1 from write.table
  if (file == "")
      file <- stdout()
  else if (is.character(file)) {
      file <- file(file, ifelse(append, "a", "w"))
      on.exit(close(file))
  }
  else if (!isOpen(file, "w")) {
      open(file, "w")
      on.exit(close(file))
  }
  if (!inherits(file, "connection"))
      stop("'file' must be a character string or connection")
  # } end 1:1 from write.table

  wt.args$file <- file

  N <- 0L

  # begin first.rows

  if (VERBOSE){
    cat("write.table.ffdf ", N+1L,"..", sep="")
    read.start <- proc.time()[3]
  }

  i <- hi(1L,first.rows)
  y <- x[i,,drop=FALSE]
  n <- first.rows
  N <- n

  if (VERBOSE){
    write.start <- proc.time()[3]
    read.time <- read.time + (write.start - read.start)
      cat(N, " (", n, ", ", round(N/nrows*100, 1), "%)  ffdf-read=", round(write.start-read.start, 3), "sec", sep="")
  }

  if (is.null(transFUN))
    dat <- do.call(FUN, c(list(y), wt.args))
  else
    dat <- do.call(FUN, c(list(transFUN(y)), wt.args))

  if (VERBOSE){
    write.stop <- proc.time()[3]
    write.time <- write.time + (write.stop - write.start)
    cat(" csv-write=", round(write.stop-write.start, 3), "sec\n", sep="")
  }

  # end first.rows

  # avoid col.names arg in R's wrapper functions
  if (!(FUN %in% c("write.csv","write.csv2")))
    wt.args$col.names <- FALSE

  wt.args$append <- TRUE

  if (N==first.rows)
  while(TRUE){
    if (N+next.rows > nrows)
      next.rows <- nrows - N
    if (next.rows<1L)
      break

    if (VERBOSE){
      cat("write.table.ffdf ", N+1L,"..", sep="")
      read.start <- proc.time()[3]
    }

    i <- hi(N+1L,N+next.rows)
    y <- x[i,,drop=FALSE]
    n <- next.rows
    N <- N + n

    if (VERBOSE){
      write.start <- proc.time()[3]
      read.time <- read.time + (write.start - read.start)
      cat(N, " (", n, ", ", round(N/nrows*100, 1), "%)  ffdf-read=", round(write.start-read.start, 3), "sec", sep="")
    }

    dat <- do.call(FUN, c(list(y), wt.args))

    if (VERBOSE){
      write.stop <- proc.time()[3]
      write.time <- write.time + (write.stop - write.start)
      cat(" csv-write=", round(write.stop-write.start, 3), "sec\n", sep="")
    }
  }

  if (VERBOSE){
    cat(" ffdf-read=", round(read.time, 3), "sec  csv-write=", round(write.time, 3), "sec  TOTAL=", round(read.time+write.time, 3), "sec\n", sep="")
  }

  return(invisible())
}


read.csv.ffdf <- function(...)
  read.table.ffdf(FUN="read.csv", ...)


read.csv2.ffdf <- function(...)
  read.table.ffdf(FUN="read.csv2", ...)

read.delim.ffdf <- function(...)
  read.table.ffdf(FUN="read.delim", ...)

read.delim2.ffdf <- function(...)
  read.table.ffdf(FUN="read.delim2", ...)

packageStartupMessage("we fix write.csv and write.csv2 such that:")
packageStartupMessage("- match.call identifies abbreviated arguments (matching against write.table)")
packageStartupMessage("- we allow col.names=FALSE in combination with append=TRUE")
write.csv <-
function (...)
{
    Call <- match.call(write.table, expand.dots = TRUE)
    for (argname in c("col.names", "sep", "dec", "qmethod")) if (!is.null(Call[[argname]]))
        warning(gettextf("attempt to set '%s' ignored", argname),
            domain = NA)
    rn <- eval.parent(Call$row.names)
    ap <- eval.parent(Call$append)
    Call$col.names <- if (is.logical(ap) && ap) FALSE else {if (is.logical(rn) && !rn) TRUE else NA}
    Call$sep <- ","
    Call$dec <- "."
    Call$qmethod <- "double"
    Call[[1L]] <- as.name("write.table")
    eval.parent(Call)
}
write.csv2 <-
function (...)
{
    Call <- match.call(write.table, expand.dots = TRUE)
    for (argname in c("col.names", "sep", "dec", "qmethod")) if (!is.null(Call[[argname]]))
        warning(gettextf("attempt to set '%s' ignored", argname),
            domain = NA)
    rn <- eval.parent(Call$row.names)
    ap <- eval.parent(Call$append)
    Call$col.names <- if (is.logical(ap) && ap) FALSE else {if (is.logical(rn) && !rn) TRUE else NA}
    Call$sep <- ";"
    Call$dec <- ","
    Call$qmethod <- "double"
    Call[[1L]] <- as.name("write.table")
    eval.parent(Call)
}

write.csv.ffdf <- function(...)
  write.table.ffdf(FUN="write.csv", ...)

write.csv2.ffdf <- function(...)
  write.table.ffdf(FUN="write.csv2", ...)


