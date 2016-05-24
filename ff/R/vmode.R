# vmode = virtual mode = attribute that both, ff and R RAM objects have
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-10-15
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/vmode.R")

#! \name{vmode}
#! \alias{vmode}
#! \alias{vmode.default}
#! \alias{vmode.ff}
#! \alias{vmode<-}
#! \alias{vmode<-.default}
#! \alias{vmode<-.ff}
#! \alias{.vmode}
#! \alias{.vunsigned}
#! \alias{.vvalues}
#! \alias{.vimplemented}
#! \alias{.rammode}
#! \alias{.ffmode}
#! \alias{.vmin}
#! \alias{.vmax}
#! \alias{.vNA}
#! \alias{.rambytes}
#! \alias{.ffbytes}
#! \alias{.vcoerceable}
#! \alias{regtest.vmode}
#! \title{ Virtual storage mode }
#! \description{
#!   Function \command{vmode} returns virtual storage modes of 'ram' or 'ff' objects, the generic \command{vmode<-} sets the vmode of ram objects (vmode of ff objects cannot be changed).
#! }
#! \usage{
#! vmode(x, ...)
#! vmode(x) <- value
#! \method{vmode}{default}(x, \dots)
#! \method{vmode}{ff}(x, \dots)
#! \method{vmode}{default}(x) <- value
#! \method{vmode}{ff}(x) <- value
#!  regtest.vmode()
#! }
#! \arguments{
#!   \item{x}{ any object }
#!   \item{value}{ a vmode from .vmode }
#!   \item{\dots}{ The \code{\dots} don't have a function yet, they are only defined to keep the generic flexible. }
#! }
#! \details{
#!  \command{vmode} is generic with default and ff methods. The following meta data vectors can be queried by .vmode or .ffmode:
#!  \tabular{rl}{
#!   \code{.vmode}         \tab virtual mode \cr
#!   \code{.vunsigned}     \tab TRUE if unsigned vmode \cr
#!   \code{.vvalues}       \tab number of possible values (incl. NA) \cr
#!   \code{.vimplemented}  \tab TRUE if this vmode is available in ff (initialized \code{\link{.onLoad}} and stored in \code{\link{globalenv}} ) \cr
#!   \code{.rammode}       \tab storage mode of this vmode \cr
#!   \code{.ffmode}        \tab integer used to code the vmode in C-code \cr
#!   \code{.vvalues}       \tab number of possible integers incl. NA in this vmode (or NA for other vmodes) \cr
#!   \code{.vmin}          \tab min integer in this vmode (or NA for other vmodes) \cr
#!   \code{.vmax}          \tab max integer in this vmode (or NA for other vmodes) \cr
#!   \code{.vNA}           \tab NA or 0 if no NA for this vmode \cr
#!   \code{.rambytes}      \tab bytes needed in ram \cr
#!   \code{.ffbytes}       \tab bytes needed by ff on disk \cr
#!   \code{.vcoerceable}   \tab list of vectors with those vmodes that can absorb this vmode \cr
#!  }
#!  the following functions relate to vmode:
#!  \tabular{rl}{
#!   \code{\link{vector.vmode}}   \tab creating (ram) vector of some vmode \cr
#!   \code{\link{as.vmode}}       \tab generic for coercing to some vmode (dropping other attributes) \cr
#!   \code{vmode<-}               \tab generic for coercing to some vmode (keeping other attributes) \cr
#!   \code{\link{maxffmode}}      \tab determine lowest \code{.ffmode} that can absorb all input vmodes without information loss \cr
#!  }
#!  some of those call the vmode-specific functions:
#!  \tabular{lll}{
#!   \strong{creation}        \tab \strong{coercion}           \tab  \strong{vmode description} \cr
#!   \code{\link{boolean}}    \tab \code{\link{as.boolean}}    \tab  1 bit logical without NA \cr
#!   \code{\link{logical}}    \tab \code{\link{as.logical}}    \tab  2 bit logical with NA \cr
#!   \code{\link{quad}}       \tab \code{\link{as.quad}}       \tab  2 bit unsigned integer without NA \cr
#!   \code{\link{nibble}}     \tab \code{\link{as.nibble}}     \tab  4 bit unsigned integer without NA \cr
#!   \code{\link{byte}}       \tab \code{\link{as.byte}}       \tab  8 bit signed integer with NA      \cr
#!   \code{\link{ubyte}}      \tab \code{\link{as.ubyte}}      \tab  8 bit unsigned integer without NA \cr
#!   \code{\link{short}}      \tab \code{\link{as.short}}      \tab 16 bit signed integer with NA      \cr
#!   \code{\link{ushort}}     \tab \code{\link{as.ushort}}     \tab 16 bit unsigned integer without NA \cr
#!   \code{\link{integer}}    \tab \code{\link{as.integer}}    \tab 32 bit signed integer with NA      \cr
#!   \code{\link{single}}     \tab \code{\link{as.single}}     \tab 32 bit float \cr
#!   \code{\link{double}}     \tab \code{\link{as.double}}     \tab 64 bit float \cr
#!   \code{\link{complex}}    \tab \code{\link{as.complex}}    \tab 2x64 bit float \cr
#!   \code{\link{raw}}        \tab \code{\link{as.raw}}        \tab 8 bit unsigned char \cr
#!   \code{\link{character}}  \tab \code{\link{as.character}}  \tab character \cr
#!  }
#! }
#! \value{
#!   \code{vmode} returns a character scalar from \code{.vmode} or "NULL" for NULL \cr
#!   \code{rambytes} returns a vector of byte counts required by each of the vmodes
#! }
#! \note{ \command{regtest.vmode} checks correctness of some vmode features
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{ff}}, \code{\link{storage.mode}}, \code{\link{mode}} }
#! \examples{
#!  data.frame(.vmode=.vmode, .vimplemented=.vimplemented, .rammode=.rammode, .ffmode=.ffmode
#! , .vmin=.vmin, .vmax=.vmax, .vNA=.vNA, .rambytes=.rambytes, .ffbytes=.ffbytes)
#!   vmode(1)
#!   vmode(1L)
#!   .vcoerceable[["byte"]]
#!   .vcoerceable[["ubyte"]]
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }


#! \name{vector.vmode}
#! \alias{vector.vmode}
#! \alias{vector.vmode.default}
#! \alias{vector.vmode.ff}
#! \alias{boolean}
#! \alias{quad}
#! \alias{nibble}
#! \alias{byte}
#! \alias{ubyte}
#! \alias{short}
#! \alias{ushort}
#! \title{ Create vector of virtual mode }
#! \description{
#!   \command{vector.vmode} creates a vector of a given vmode and length
#! }
#! \usage{
#! vector.vmode(vmode = "logical", length = 0)
#! boolean(length = 0)
#! quad(length = 0)
#! nibble(length = 0)
#! byte(length = 0)
#! ubyte(length = 0)
#! short(length = 0)
#! ushort(length = 0)
#! }
#! \arguments{
#!   \item{vmode}{ virtual mode }
#!   \item{length}{ desired length }
#! }
#! \details{
#!  Function \command{vector.vmode} creates the vector in one of the usual \code{\link{storage.mode}s} (see \code{\link{.rammode}}) but flags them with an additional attribute 'vmode' if necessary.
#!  The creators can also be used directly:
#!  \tabular{rl}{
#!   \code{boolean}    \tab  1 bit logical without NA \cr
#!   \code{logical}    \tab  2 bit logical with NA \cr
#!   \code{quad}       \tab  2 bit unsigned integer without NA \cr
#!   \code{nibble}     \tab  4 bit unsigned integer without NA \cr
#!   \code{byte}       \tab  8 bit signed integer with NA      \cr
#!   \code{ubyte}      \tab  8 bit unsigned integer without NA \cr
#!   \code{short}      \tab 16 bit signed integer with NA      \cr
#!   \code{ushort}     \tab 16 bit unsigned integer without NA \cr
#!   \code{integer}    \tab 32 bit signed integer with NA       \cr
#!   \code{single}     \tab 32 bit float \cr
#!   \code{double}     \tab 64 bit float \cr
#!   \code{complex}    \tab 2x64 bit float \cr
#!   \code{raw}        \tab 8 bit unsigned char \cr
#!   \code{character}  \tab character \cr
#!  }
#! }
#! \value{
#!  a vector of the desired vmode initialized with 0
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{as.vmode}}, \code{\link{vector}} }
#! \examples{
#!   vector.vmode("byte",12)
#!   vector.vmode("double",12)
#!   byte(12)
#!   double(12)
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }


#! \name{as.vmode}
#! \alias{as.vmode}
#! \alias{as.vmode.default}
#! \alias{as.vmode.ff}
#! \alias{as.boolean}
#! \alias{as.boolean.default}
#! \alias{as.quad}
#! \alias{as.quad.default}
#! \alias{as.nibble}
#! \alias{as.nibble.default}
#! \alias{as.byte}
#! \alias{as.byte.default}
#! \alias{as.ubyte}
#! \alias{as.ubyte.default}
#! \alias{as.short}
#! \alias{as.short.default}
#! \alias{as.ushort}
#! \alias{as.ushort.default}
#! \title{ Coercing to virtual mode }
#! \description{
#!   \command{as.vmode} is a generic that converts some R ram object to the desired \code{\link{vmode}}.
#! }
#! \usage{
#! as.vmode(x, ...)
#! as.boolean(x, ...)
#! as.quad(x, ...)
#! as.nibble(x, ...)
#! as.byte(x, ...)
#! as.ubyte(x, ...)
#! as.short(x, ...)
#! as.ushort(x, ...)
#! \method{as.vmode}{default}(x, vmode, \dots)
#! \method{as.vmode}{ff}(x, \dots)
#! \method{as.boolean}{default}(x, ...)
#! \method{as.quad}{default}(x, ...)
#! \method{as.nibble}{default}(x, ...)
#! \method{as.byte}{default}(x, ...)
#! \method{as.ubyte}{default}(x, ...)
#! \method{as.short}{default}(x, ...)
#! \method{as.ushort}{default}(x, ...)
#! }
#! \arguments{
#!   \item{x}{ any object }
#!   \item{vmode}{ virtual mode }
#!   \item{\dots}{ The \code{\dots} don't have a function yet, they are only defined to keep the generic flexible. }
#! }
#! \details{
#!  Function \command{as.vmode} actually coerces to one of the usual \code{\link{storage.mode}s} (see \code{\link{.rammode}}) but flags them with an additional attribute 'vmode' if necessary.
#!  The coercion generics can also be called directly:
#!  \tabular{rl}{
#!   \code{as.boolean}    \tab  1 bit logical without NA \cr
#!   \code{as.logical}    \tab  2 bit logical with NA \cr
#!   \code{as.quad}       \tab  2 bit unsigned integer without NA \cr
#!   \code{as.nibble}     \tab  4 bit unsigned integer without NA \cr
#!   \code{as.byte}       \tab  8 bit signed integer with NA      \cr
#!   \code{as.ubyte}      \tab  8 bit unsigned integer without NA \cr
#!   \code{as.short}      \tab 16 bit signed integer with NA      \cr
#!   \code{as.ushort}     \tab 16 bit unsigned integer without NA \cr
#!   \code{as.integer}    \tab 32 bit signed integer with NA      \cr
#!   \code{as.single}     \tab 32 bit float \cr
#!   \code{as.double}     \tab 64 bit float \cr
#!   \code{as.complex}    \tab 2x64 bit float \cr
#!   \code{as.raw}        \tab 8 bit unsigned char \cr
#!   \code{as.character}  \tab character \cr
#!  }
#! }
#! \value{
#!   a vector of the desired vmode containing the input data
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{vmode}}, \code{\link{vector.vmode}} }
#! \examples{
#!   as.vmode(1:3,"double")
#!   as.vmode(1:3,"byte")
#!   as.double(1:3)
#!   as.byte(1:3)
#! }
#! \keyword{ IO }
#! \keyword{ data }
#! \keyword{ attribute }


.vmode <- c(
  boolean   = "boolean"   # {FALSE,TRUE} no NA
, logical   = "logical"   # {FALSE,TRUE,NA}
, quad      = "quad"      # int2bit
, nibble    = "nibble"    # int4bit
, byte      = "byte"      # int8bit
, ubyte     = "ubyte"     # uint8bit
, short     = "short"     # int16bit
, ushort    = "ushort"    # uint16bit
, integer   = "integer"   # int32bit
, single    = "single"    # float32bit
, double    = "double"    # float64bit
, complex   = "complex"   # struct 2x float64bit
, raw       = "raw"       # unsigned char 0:255, no NA
, character = "character" # not used but coded
)

.vset <- list(
  boolean   = function(x)setattr(x, "vmode", "boolean")
, logical   = function(x)setattr(x, "vmode", NULL)
, quad      = function(x)setattr(x, "vmode", "quad")
, nibble    = function(x)setattr(x, "vmode", "nibble")
, byte      = function(x)setattr(x, "vmode", "byte")
, ubyte     = function(x)setattr(x, "vmode", "ubyte")
, short     = function(x)setattr(x, "vmode", "short")
, ushort    = function(x)setattr(x, "vmode", "ushort")
, integer   = function(x)setattr(x, "vmode", NULL)
, single    = function(x){setattr(x, "vmode", NULL); setattr(x, "Csingle", TRUE)}
, double    = function(x)setattr(x, "vmode", NULL)
, complex   = function(x)setattr(x, "vmode", NULL)
, raw       = function(x)setattr(x, "vmode", NULL)
, character = function(x)setattr(x, "vmode", NULL)
)

.vunsigned <- c(
  boolean   = TRUE
, logical   = FALSE
, quad      = TRUE
, nibble    = TRUE
, byte      = FALSE
, ubyte     = TRUE
, short     = FALSE
, ushort    = TRUE
, integer   = FALSE
, single    = FALSE
, double    = FALSE
, complex   = FALSE
, raw       = FALSE
, character = FALSE
)


.vvalues <- c(
  boolean   = 2
, logical   = 3
, quad      = 4
, nibble    = 16
, byte      = 256
, ubyte     = 256
, short     = 65536
, ushort    = 65536
, integer   = 4294967296
, single    = NA
, double    = NA
, complex   = NA
, raw       = 256
, character = NA
)


.rammode <- c(
  boolean   = "logical"
, logical   = "logical"
, quad      = "integer"
, nibble    = "integer"
, byte      = "integer"
, ubyte     = "integer"
, short     = "integer"
, ushort    = "integer"
, integer   = "integer"
, single    = "double"
, double    = "double"
, complex   = "complex"
, raw       = "raw"
, character = "character"
)

# these codes are used in r_api.c to switch the datatype
.ffmode <- c(
  boolean   =  1L         # {FALSE,TRUE} no NA
, logical   =  2L         # {FALSE,TRUE,NA}
, quad      =  3L         # int2bit
, nibble    =  4L         # int4bit
, byte      =  5L         # int8bit
, ubyte     =  6L         # uint8bit
, short     =  7L         # int16bit
, ushort    =  8L         # uint16bit
, integer   =  9L         # int32bit
, single    = 10L         # float32bit
, double    = 11L         # float64bit
, complex   = 12L         # struct 2x float64bit
, raw       = 13L         # unsigned char 0:255, no NA
, character = 14L         # not used but coded
)

.vmin <- c(
  boolean   = 0L
, logical   = 0L
, quad      = 0L
, nibble    = 0L
, byte      = -127L
, ubyte     = 0L
, short     = -32767L
, ushort    = 0L
, integer   = -2147483647L
, single    = NA
, double    = NA
, complex   = NA
, raw       = 0L
, character = NA
)

.vmax <- c(
  boolean   = 1L
, logical   = 1L
, quad      = 3L
, nibble    = 15L
, byte      = 127L
, ubyte     = 255L
, short     = 32767L
, ushort    = 65535L
, integer   = 2147483647L
, single    = NA
, double    = NA
, complex   = NA
, raw       = 255L
, character = NA
)

.vNA <- c(
  boolean   = 0L
, logical   = NA
, quad      = 0L
, nibble    = 0L
, byte      = NA
, ubyte     = 0L
, short     = NA
, ushort    = 0L
, integer   = NA
, single    = NA
, double    = NA
, complex   = NA
, raw       = 0L
, character = NA
)

# R bytes in RAM
.rambytes <- c(
  boolean   = 4L
, logical   = 4L
, quad      = 4L
, nibble    = 4L
, byte      = 4L
, ubyte     = 4L
, short     = 4L
, ushort    = 4L
, integer   = 4L
, single    = 8L
, double    = 8L
, complex   = 16L
, raw       = 1L
, character = 4L          # from R 2.6.0 upwards: a string is a vector of pointers to a global string hash
)

# ff bytes on disk
.ffbytes <- c(
  boolean   = 0.125
, logical   = 0.25
, quad      = 0.25
, nibble    = 0.5
, byte      = 1
, ubyte     = 1
, short     = 2
, ushort    = 2
, integer   = 4
, single    = 4
, double    = 8
, complex   = 16
, raw       = 1
, character = NA
)

# Note: .vimplemented is finally initialized in .onLoad using a .Call into the dll
.vimplemented <- c(
 boolean    = FALSE
, logical   = FALSE
, quad      = FALSE
, nibble    = FALSE
, byte      = FALSE
, ubyte     = FALSE
, short     = FALSE
, ushort    = FALSE
, integer   = FALSE
, single    = FALSE
, double    = FALSE
, complex   = FALSE
, raw       = FALSE
, character = FALSE
)

# Note: .vcoerceable is finally initialized in .onLoad using a .Call into the dll
.vcoerceable <- NULL 

#xx
#rambytes <- function(vmode){
#  ret <- .rambytes[vmode]
#  isFixedWidth <- grepl("char[0123456789]",vmode)
#  if (any(isFixedWidth)){
#    names(ret) <- vmode
#    vmode <- vmode[isFixedWidth]
#    ret[isFixedWidth] <- as.integer(substr(vmode, 5, nchar(vmode)))
#  }
#  ret
#}


vector.vmode <- function(vmode="logical", length=0){
    if (is.null(vmode) || vmode=="NULL")
      stop("you can't create a vector with vmode='NULL'")
    makevmode <- list(
      boolean   = boolean
    , logical   = logical
    , quad      = quad
    , nibble    = nibble
    , byte      = byte
    , ubyte     = ubyte
    , short     = short
    , ushort    = ushort
    , integer   = integer
    , single    = single
    , double    = double
    , complex   = complex
    , raw       = raw
    , character = character
    )[[vmode]]
    makevmode(length)
}


# source("d:/mwp/eanalysis/ff/R/vmode.R")
vmode.default <- function(x, ...){
  rammode <- storage.mode(x)
  v <- attr(x,"vmode")
  if (is.null(v)){
    s <- attr(x, "Csingle")
    if (!is.null(s) && s)
      return("single")
    else
      return(rammode)
  }else{
    if (.rammode[v]!=rammode)
      stop("storage.mode does not match vmode; most likely the ramobject was silently coerced through lazy assignment of value with 'higher' storage.mode")
    return(v)
  }
}

"vmode<-.default" <-
function (x, value)
{
  if (is.null(value) || value=="NULL")
    stop("you can't create a vector with vmode value='NULL'")
  vm <- vmode(x)
  if (vm=="NULL")
    stop("you can't coerce NULL to a different vmode")
  if (vm==value){
    x
  }else{

    atr <- attributes(x)
    atr$Csingle <- NULL
    atr$vmode <- NULL

    asvmode <- list(
      boolean   = as.boolean
    , logical   = as.logical
    , quad      = as.quad
    , nibble    = as.nibble
    , byte      = as.byte
    , ubyte     = as.ubyte
    , short     = as.short
    , ushort    = as.ushort
    , integer   = as.integer
    , single    = as.single
    , double    = as.double
    , complex   = as.complex
    , raw       = as.raw
    , character = as.character
    )[[value]]
    x <- asvmode(x)

    atr2 <- attributes(x)
    atr[names(atr2)] <- atr2
    attributes(x) <- atr

    x
  }
}

as.vmode.default <- function(x, vmode
, ... # dummy to keep R CMD check quiet
){
  if (is.null(vmode) || vmode=="NULL")
    stop("you can't create a vector with vmode='NULL'")
  vm <- vmode(x)
  if (vm=="NULL")
    stop("you can't coerce NULL to a different vmode")
  if (vm==vmode){
    x
  }else{
    asvmode <- list(
      boolean   = as.boolean
    , logical   = as.logical
    , quad      = as.quad
    , nibble    = as.nibble
    , byte      = as.byte
    , ubyte     = as.ubyte
    , short     = as.short
    , ushort    = as.ushort
    , integer   = as.integer
    , single    = as.single
    , double    = as.double
    , complex   = as.complex
    , raw       = as.raw
    , character = as.character
    )[[vmode]]
    asvmode(x)
  }
}


vmode.ff <- function(x
, ... # dummy to keep R CMD check quiet
){
  attr(attr(x, "physical"), "vmode")
}

"vmode<-.ff" <- function(x, value){
  stop("ff objects cannot be coerced to a different vmode")
}

as.vmode.ff <- function(x
, ... # dummy to keep R CMD check quiet
){
  stop("ff objects cannot be coerced to a different vmode")
}




boolean <-
function (length = 0){
  x <- vector("logical", length)
  #attr(x, "vmode") <- "boolean"
  setattr(x, "vmode", "boolean")
  x
}
as.boolean <-
function (x, ...){
  if (vmode(x)=="boolean")
    return(x)
  UseMethod("as.boolean")
}
as.boolean.default <-
function (x, ...)
{
  x <- as.vector(x, "logical")
  x[is.na(x)] <- FALSE
  attr(x, "vmode") <- "boolean"
  x
}


quad <-
function (length = 0){
  x <- vector("integer", length)
  #attr(x, "vmode") <- "quad"
  setattr(x, "vmode", "quad")
  x
}
as.quad <-
function (x, ...){
  if (vmode(x)=="quad")
    return(x)
  UseMethod("as.quad")
}
as.quad.default <-
function (x, ...)
{
  x <- as.vector(x, "integer")
  x[is.na(x)] <- 0L
  attr(x, "vmode") <- "quad"
  x
}

nibble <-
function (length = 0){
  x <- vector("integer", length)
  #attr(x, "vmode") <- "nibble"
  setattr(x, "vmode", "nibble")
  x
}
as.nibble <-
function (x, ...){
  if (vmode(x)=="unibble")
    return(x)
  UseMethod("as.nibble")
}
as.nibble.default <-
function (x, ...)
{
  x <- as.vector(x, "integer")
  x[is.na(x)] <- 0L
  attr(x, "vmode") <- "nibble"
  x
}

byte <-
function (length = 0){
  x <- vector("integer", length)
  #attr(x, "vmode") <- "byte"
  setattr(x, "vmode", "byte")
  x
}
as.byte <-
function (x, ...){
  if (vmode(x)=="ubyte")
    return(x)
  UseMethod("as.byte")
}
as.byte.default <-
function (x, ...)
{
  x <- as.vector(x, "integer")
  attr(x, "vmode") <- "byte"
  x
}

ubyte <-
function (length = 0){
  x <- vector("integer", length)
  #attr(x, "vmode") <- "ubyte"
  setattr(x, "vmode", "ubyte")
  x
}
as.ubyte <-
function (x, ...){
  if (vmode(x)=="ubyte")
    return(x)
  UseMethod("as.ubyte")
}
as.ubyte.default <-
function (x, ...)
{
  x <- as.vector(x, "integer")
  x[is.na(x)] <- 0L
  attr(x, "vmode") <- "ubyte"
  x
}

short <-
function (length = 0){
  x <- vector("integer", length)
  #attr(x, "vmode") <- "short"
  setattr(x, "vmode", "short")
  x
}
as.short <-
function (x, ...){
  if (vmode(x)=="short")
    return(x)
  UseMethod("as.short")
}
as.short.default <-
function (x, ...)
{
  x <- as.vector(x, "integer")
  attr(x, "vmode") <- "short"
  x
}

ushort <-
function (length = 0){
  x <- vector("integer", length)
  #attr(x, "vmode") <- "ushort"
  setattr(x, "vmode", "ushort")
  x
}
as.ushort <-
function (x, ...){
  if (vmode(x)=="ushort")
    return(x)
  UseMethod("as.ushort")
}
as.ushort.default <-
function (x, ...)
{
  x <- as.vector(x, "integer")
  x[is.na(x)] <- 0L
  attr(x, "vmode") <- "ushort"
  x
}




#! \name{ram2ffcode}
#! \alias{ram2ffcode}
#! \alias{ram2ramcode}
#! \title{ Factor codings }
#! \description{
#!   Function \command{ram2ffcode} creates the \emph{internal} factor codes used by ff to store factor levels. Function \code{ram2ramcode} is a compatibility function used instead if \code{RETURN_FF==FALSE}.
#! }
#! \usage{
#! ram2ffcode(value, levels, vmode)
#! ram2ramcode(value, levels)
#! }
#! \arguments{
#!   \item{value}{ factor or character vector of values }
#!   \item{levels}{ character vector of factor levels }
#!   \item{vmode}{ one of the integer vmodes in \code{\link{.rammode}} }
#! }
#! \details{
#!   Factors stored in unsigned vmodes \code{\link{.vunsigned}} have their first level represented as 0L instead of 1L.
#! }
#! \value{
#!   A vector of integer values representing the correspnding factor levels.
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{factor}}, \code{\link{levels.ff}}, \code{\link{vmode}} }
#! \examples{
#!  ram2ffcode(letters, letters, vmode="byte")
#!  ram2ffcode(letters, letters, vmode="ubyte")
#!  ram2ffcode(letters, letters, vmode="nibble")
#!  message('note that ram2ffcode() does NOT warn that vmode="nibble" cannot store 26 levels')
#! }
#! \keyword{ IO }
#! \keyword{ data }

ram2ffcode <- function(value, levels, vmode){
  if (.vunsigned[vmode]){
    if (is.factor(value)){
      if (identical(levels(value), levels)){
        ret <- unclass(value) - 1L
      }else{
        ret <- match(as.character(value), levels) - 1L
      }
    }else{
      ret <- match(value, levels) - 1L
    }
    if (any(is.na(ret)))
      warning("unknown factor values mapped to first level (no NAs in unsigned vmode)")
  }else{
    if (is.factor(value)){
      if (identical(levels(value), levels)){
        ret <- unclass(value)
      }else{
        ret <- match(as.character(value), levels)
      }
    }else{
      ret <- match(value, levels)
    }
    if (any(is.na(ret) & !is.na(value)))
      warning("unknown factor values mapped to NA")
  }
  ret
}
ram2ramcode <- function(value, levels){
  if (is.factor(value)){
    if (identical(levels(value), levels)){
      ret <- unclass(value)
    }else{
      ret <- match(as.character(value), levels)
    }
  }else{
    ret <- match(value, levels)
  }
  if (any(is.na(ret)) & !is.na(value))
    warning("unknown factor values mapped to NA")
  ret
}


#! \name{maxffmode}
#! \alias{maxffmode}
#! \title{ Lossless vmode coercability }
#! \description{
#!   \command{maxffmode} returns the lowest \code{\link{vmode}} that can absorb all input vmodes without data loss
#! }
#! \usage{
#! maxffmode(...)
#! }
#! \arguments{
#!   \item{\dots}{ one or more vectors of vmodes }
#! }
#! \value{
#!   the smallest \code{\link{.ffmode}} which can absorb the input vmodes without data loss
#! }
#! \note{
#!   The output can be larger than any of the inputs (if the highest input vmode is an integer type without NA and any other input requires NA).
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{  \code{\link{.vcoerceable}}, \code{\link{.ffmode}}, \code{\link{ffconform}} }
#! \examples{
#!    maxffmode(c("quad","logical"), "ushort")
#! }
#! \keyword{ IO }
#! \keyword{ data }

# returns lowest ffmode to which all arguments can be coerced without information loss or zero if not coerceable
# NOTE: vmode is in names of return value
maxffmode <- function(...){
  l <- unlist(lapply(list(...), function(x).ffmode[x]))
  n <- length(l)
  if (n==1)
    return(.ffmode[l[[1]]])
  if (n){
    m <- .vcoerceable[[l[[1]]]]
    if (n>1){
      for (i in 2:n){
        m <- intersect(m, .vcoerceable[[l[[i]]]])
      }
      if (length(m))
        m <- min(m)
      else
        return(c("not .vcoerceable"=0))
    }
  }else{
    m <- 2L
  }
  return(.ffmode[m])
}


if (FALSE){
  # xx experimental, not used yet
  charmode <- function(x, n=.Machine$integer.max){
      n <- min(length(x), n)
      if (n<1) return("character")
      nc <- unique(nchar(x[1:n]))
      if (length(nc)==1)
        paste("char",nc,sep="")
      else
        "character"
  }
}

regtest.vmode <- function(){
  testmodes <- .vmode[!is.na(.vmin)]
  testmodes <- testmodes[.vimplemented[testmodes]]
  #testmodes <- setdiff(testmodes, "integer")
  OK <- TRUE
  for (v in testmodes){
    message(v)
    if (v=="raw")
      oldopt <- options(warn=-1)
    if (v=="integer")
      x <- as.vmode(c(NA, c(.vmin[v]:(.vmin[v]+4), (.vmax[v]-4):.vmax[v])), v)
    else
      x <- as.vmode(c(NA, .vmin[v]:.vmax[v]), v)
    if (v=="raw")
      options(oldopt)
    a <- ff(x)
    if (!identical(as.integer(a[1]), as.vector(.vNA[v]))){
      warning('NA mapped to wrong value in vmode "',v,'"')
      OK <- FALSE
    }
    if (!identical(x[-1], a[-1])){
      warning('Non-NA mapped to wrong value in vmode "',v,'"')
      OK <- FALSE
    }
    if (v!="raw")
      x[1] <- .vNA[v]
    a[,add=TRUE] <- 0L
    if (!identical(as.integer(x[1:length(a)]), as.integer(a[1:length(a)]))){
      warning('+=0 does not give the expected result in vmode "',v,'"')
      OK <- FALSE
    }
    a[,add=TRUE] <- 1L
    if (v=="raw"){
      x2 <- c(1L,1:255,0L)
    }else if(v=="boolean"){
      x2 <- c(1L,1L,0L)       # toggle 0<->1
    }else if(v=="logical"){
      x2 <- c(NA,1L,0L)       # toggle 0<->1 (non-NAs)
    }else{
      x2 <- ifelse(x==.Machine$integer.max, NA, x)
      x2 <- x2 + 1L
      if (v %in% c("byte","short"))
        x2[-1][x2[-1]>.vmax[v]] <- NA         # signed integer types overflow to NA
      else
        x2[-1][x2[-1]>.vmax[v]] <- .vmin[v]   # unsigned integer types wrap-around
    }
    if (!identical(as.vector(x2), as.integer(a[1:length(a)]))){
      warning('+=1 does not give the expected result in vmode "',v,'"')
      OK <- FALSE
    }
    a[] <- x
    if (v=="raw")
      oldopt <- options(warn=-1)
    a[,add=TRUE] <- as.integer(NA)
    if (v=="raw")
      options(oldopt)
    if (v %in% c("boolean","quad","nibble","ubyte","ushort","raw")){
      x2 <- as.integer(x[1:length(a)])
    }else{
      x2 <- as.integer(rep(NA, length(a)))
    }
    if (!identical(as.vector(x2), as.integer(a[1:length(a)]))){
      warning('+=NA does not give the expected result in vmode "',v,'"')
      OK <- FALSE
    }
  }
  OK
}
#regtest.vmode()



if (FALSE){
  regtest.vmode()
}
