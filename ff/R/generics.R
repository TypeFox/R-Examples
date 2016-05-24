# New generic functions for ff (and R.ff)
# (c) 2007 Jens Oehlschlägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-08-24
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/generics.R")

# coerce to hybrid index (part of the HIP concept)
as.hi <- function(x, ...){
  UseMethod("as.hi")
}

# physical and virtual features of an object
if (!exists("physical"))
  physical <- function(x)UseMethod("physical")
if (!exists("physical<-"))
  "physical<-" <- function(x, value)UseMethod("physical<-")
if (!exists("virtual"))
  virtual <- function(x)UseMethod("virtual")
if (!exists("virtual<-"))
  "virtual<-" <- function(x, value)UseMethod("virtual<-")

# complementing 'length' generic
if (!exists("maxindex"))
  maxindex  <- function(x, ...)UseMethod("maxindex")  # max possible index length (required for negative indices)
if (!exists("poslength"))
  poslength <- function(x, ...)UseMethod("poslength") # length of selected elements (even for negative indices)
maxlength <- function(x, ...)UseMethod("maxlength") # max physical length
#maxwidth <- function (x, ...)UseMethod("maxwidth")  # max width of fixed char (fffc)



# further physical and readonly attributes
filename    <- function(x, ...)UseMethod("filename")    # (physical attribute, part of the HCS concept)
"filename<-" <- function(x, ..., value)UseMethod("filename<-")    # (physical attribute, part of the HCS concept)
pattern    <- function(x, ...)UseMethod("pattern")    # (physical attribute, part of the HCS concept)
"pattern<-" <- function(x, ..., value)UseMethod("pattern<-")    # (physical attribute, part of the HCS concept)
is.readonly <- function(x, ...)UseMethod("is.readonly") # (physical attribute, part of the HCS concept)
is.open     <- function(x, ...)UseMethod("is.open")     # readonly
pagesize    <- function(x, ...)UseMethod("pagesize")     # readonly

# querying and coercing vmode (physical attribute, part of the HCS concept)
    vmode <- function(x, ...)UseMethod("vmode")
 as.vmode <- function(x, ...)UseMethod("as.vmode")
"vmode<-" <- function(x, value)UseMethod("vmode<-")

# finalization
finalizer   <- function(x, ...)UseMethod("finalizer")               # get name of active finalizer
"finalizer<-" <- function(x, ..., value)UseMethod("finalizer<-")    # change and activate finalizer (but not: de-activate)
finalize  <- function(x, ...)UseMethod("finalize")                  # call active finalizer (if any)


# generics complementing 'open' and 'close' in ff file handling
delete     <- function(x, ...)UseMethod("delete")
deleteIfOpen <- function(x, ...)UseMethod("deleteIfOpen")

# complements 'ff' and 'update': cloning and convenience wrappers for swap in/out of cache
clone  <- function(x, ...)UseMethod("clone")
as.ff  <- function(x, ...)UseMethod("as.ff")
as.ram <- function(x, ...)UseMethod("as.ram")

# super-classes
#as.fffc <- function(x, ...)UseMethod("as.fffc")
as.ffdf <- function(x, ...)UseMethod("as.ffdf")


# querying and setting dimorder (virtual attribute, part of the HCS concept)
 dimorder    <- function(x, ...)UseMethod("dimorder")
"dimorder<-" <- function(x, ..., value)UseMethod("dimorder<-")

# querying and setting virtual windows (virtual attribute, part of the HCS concept)
 vw    <- function(x, ...)UseMethod("vw")
"vw<-" <- function(x, ..., value)UseMethod("vw<-")

# querying is.factor and is.ordered
if (!exists("is.factor.default")){
  packageStartupMessage("creating generic 'is.factor'")
  is.factor.default <- get("is.factor", mode = "function") #pos = NULL,
  is.factor <- function (x)
    UseMethod("is.factor")
}
if (!exists("is.ordered.default")){
  packageStartupMessage("creating generic 'is.ordered'")
  is.ordered.default <- get("is.ordered", mode = "function") #pos = NULL,
  is.ordered <- function (x)
    UseMethod("is.ordered")
}

ramclass <- function (x, ...)
  UseMethod("ramclass")
ramattribs <- function (x, ...)
  UseMethod("ramattribs")

# fffactor
recodeLevels <- function(x, lev)
  UseMethod("recodeLevels")

sortLevels <- function(x)
  UseMethod("sortLevels")



# -- special access methods --

# add(x, value)  <==>  x[,add=TRUE] <- value  <==>  x += value
add <- function(x, ...)UseMethod("add")

# swap: combines existing generics '[' and '[<-' and is needed e.g. for maintaining NA counts
swap <- function(x, value, ...)UseMethod("swap")


# -- for future use in ff ------------------------------------------------------------

# check for 'hard' symmetric attribute as contrasted with R's isSymmetric() that checks for symmetry on the fly via all.equal()
symmetric <- function(x, ...)UseMethod("symmetric")

# read and write fixed diagonal in symmetric matrices
fixdiag <- function(x, ...)UseMethod("fixdiag")
"fixdiag<-" <- function(x, ..., value)UseMethod("fixdiag<-")


# -- R.ff stuff currently already in ff

# large sampling
bigsample <- function(x, ...)UseMethod("bigsample")

# virtual matrix transpose
vt <- function(x, ...)UseMethod("vt")


# -- not yet generics ----------------------------------------------------------------

#! \name{nrowAssign}
#! \Rdversion{1.1}
#! \alias{nrow<-}
#! \alias{ncol<-}
#! \title{
#!   Assigning the number of rows or columns
#! }
#! \description{
#!   Function \code{nrow<-} assigns \code{\link[base]{dim}} with a new number of rows. \cr
#!   Function \code{ncol<-} assigns \code{\link[base]{dim}} with a new number of columns.
#! }
#! \usage{
#! nrow(x) <- value
#! ncol(x) <- value
#! }
#! \arguments{
#!   \item{x}{ a object that has \code{\link[base]{dim}} AND can be assigned ONE new dimension }
#!   \item{value}{ the new size of the assigned dimension }
#! }
#! \details{
#!   Currently only asssigning new rows to \code{\link{ffdf}} is supported.
#!   The new ffdf rows are not initialized (usually become zero).
#!   NOTE that
#! }
#! \value{
#!   The object with a modified dimension
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{ffdf}}, \code{\link{dim.ffdf}}
#! }
#! \examples{
#!   a <- as.ff(1:26)
#!   b <- as.ff(factor(letters)) # vmode="integer"
#!   c <- as.ff(factor(letters), vmode="ubyte")
#!   df <- ffdf(a,b,c)
#!   nrow(df) <- 2*26
#!   df
#!   message("NOTE that the new rows have silently the first level 'a' for UNSIGNED vmodes")
#!   message("NOTE that the new rows have an illegal factor level <0> for SIGNED vmodes")
#!   message("It is your responsibility to put meaningful content here")
#!   message("As an example we replace the illegal zeros by NA")
#!   df$b[27:52] <- NA
#!   df
#!
#!   rm(a,b,c,df); gc()
#! }
#! \keyword{ array }


"nrow<-" <- function(x, value){
  d <- dim(x)
  if (is.null(d) || length(d)!=2)
    stop("not a two-dimensional array")
  dim(x) <- c(as.integer(value), d[[2]])
  x
}

"ncol<-" <- function(x, value){
  d <- dim(x)
  if (is.null(d) || length(d)!=2)
    stop("not a two-dimensional array")
  dim(x) <- c(d[[1]], as.integer(value))
  x
}


