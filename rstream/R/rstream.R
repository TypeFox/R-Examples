#############################################################################
##                                                                         ##
##   virtual Class: rstream                                                ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   The superclass of the package                                         ##
##                                                                         ##
#############################################################################


## Initialize global variables ----------------------------------------------

## .rstream.envir
##    Environment for rstream global variables and constants.
.rstream.envir <- new.env(parent=.GlobalEnv)

.rstream.init <- function () {

	## .rstream.Count
	##    Counter for Rstream objects
	##    Used for default names of Rstream objects
	if (!exists(".rstream.Count", envir=.rstream.envir))
		assign(".rstream.Count", 0, envir=.rstream.envir)

        ## .rstream.version
        ##    Store rstream version for compatibilty function.
        ##    Use 9999000 to indicate installed version ("default")
        assign(".rstream.version", 9999000, envir=.rstream.envir)
}


## Class --------------------------------------------------------------------

setClass( "rstream", 
         ## slots:
         representation( 
                        type	= "character",	  # type of RNG Rstream
                        info	= "character",	  # short description of Rstream type
                        stream	= "externalptr",  # pointer to Rstream object
                        xstate   = "environment", #
                        is.packed = "logical",    # whether the Rstream object is packed
                        pack	= "list" ),       # contains packed data of Rstream
         ## virtual class:
         contains = "VIRTUAL") 


## Methods ------------------------------------------------------------------

## Access and Replacement methods ...........................................

## rstream.name
##    get and set name of Rstream object
if(!isGeneric("rstream.name"))
	setGeneric("rstream.name", function(stream) standardGeneric("rstream.name"))
setMethod("rstream.name", "rstream", 
          function(stream) { "unknown" } )

if(!isGeneric("rstream.name<-"))
        setGeneric("rstream.name<-", function(stream, value) standardGeneric("rstream.name<-"))
setReplaceMethod("rstream.name", "rstream", 
                 function(stream, value) {
                         stop(class(stream)[1],": No names for this Rstream class")
                         stream
                 } )


## rstream.antithetic
##   get and set flag for antithetic random numbers
if(!isGeneric("rstream.antithetic"))
	setGeneric("rstream.antithetic", function(stream) standardGeneric("rstream.antithetic"))
setMethod("rstream.antithetic", "rstream", 
          function(stream) { FALSE } )

if(!isGeneric("rstream.antithetic<-"))
        setGeneric("rstream.antithetic<-", function(stream, value) standardGeneric("rstream.antithetic<-"))
setReplaceMethod("rstream.antithetic", "rstream", 
                 function(stream, value) {
                         stop(class(stream)[1],": No antithetic variates available for this Rstream class")
                         stream
                 } )


## rstream.incprecision
##    get and set flag for increased precision of random numbers
if(!isGeneric("rstream.incprecision"))
	setGeneric("rstream.incprecision", function(stream) standardGeneric("rstream.incprecision"))
setMethod("rstream.incprecision", "rstream", 
          function(stream) { FALSE } )

if(!isGeneric("rstream.incprecision<-"))
        setGeneric("rstream.incprecision<-", function(stream, value) standardGeneric("rstream.incprecision<-"))
setReplaceMethod("rstream.incprecision", "rstream", 
                 function(stream, value) {
                         stop(class(stream)[1],": No increased precision available for this Rstream class")
                         stream
                 } )


## Sampling methods .........................................................

## rstream.sample
##    make a random sample
if(!isGeneric("rstream.sample"))
        setGeneric("rstream.sample", function(stream,...) standardGeneric("rstream.sample"))

setMethod("rstream.sample", "rstream", 
          function(stream,n=1) { 
                  stop(class(stream)[1],": Cannot sample from this Rstream class")
                  NULL
          } )


## r
##    alias for rstream.sample  (slow!!)
if(!isGeneric("r"))
        setGeneric("r", function(stream,...) standardGeneric("r"))

setMethod("r", "rstream",
          function(stream,n=1) { rstream.sample(stream,n) } )


## Jump methods .............................................................

## rstream.resetsubstream
##   reset current substream of Rstream object
if(!isGeneric("rstream.resetsubstream"))
        setGeneric("rstream.resetsubstream", function(stream) standardGeneric("rstream.resetsubstream"))
setMethod("rstream.resetsubstream", "rstream", 
          function(stream) { 
                  stop(class(stream)[1],": No substreams for this Rstream class") } )


## rstream.nextsubstream
##   skip to next substream of Rstream object
if(!isGeneric("rstream.nextsubstream"))
        setGeneric("rstream.nextsubstream", function(stream) standardGeneric("rstream.nextsubstream"))
setMethod("rstream.nextsubstream", "rstream", 
          function(stream) { 
                  stop(class(stream)[1],": No substreams for this Rstream class") } )


# Reset and copy methods ....................................................

## rstream.reset
##   reset Rstream object
if(!isGeneric("rstream.reset"))
        setGeneric("rstream.reset", function(stream) standardGeneric("rstream.reset"))
setMethod("rstream.reset", "rstream", 
          function(stream) { 
                  stop(class(stream)[1],": No reset for this Rstream class") } )


## rstream.clone
##    clone (copy) Rstream object
if(!isGeneric("rstream.clone"))
	setGeneric("rstream.clone", function(stream) standardGeneric("rstream.clone"))
setMethod("rstream.clone", "rstream", 
          function(stream) { 
                  stop(class(stream)[1],": Cannot make clone for this Rstream class")
                  stream
          } )


## rstream.pack, rstream.unpack
##    pack and unpack Rstream object such that all data are contained in R object
##    (and can be easily copied within R)
if(!isGeneric("rstream.packed"))
	setGeneric("rstream.packed", function(stream) standardGeneric("rstream.packed"))
setMethod("rstream.packed", "rstream", 
          function(stream) { stream@is.packed } )

if(!isGeneric("rstream.packed<-"))
        setGeneric("rstream.packed<-", function(stream, value) standardGeneric("rstream.packed<-"))
setReplaceMethod("rstream.packed", "rstream", 
                 function(stream, value) {
                         stop(class(stream)[1],": Cannot pack/unpack object of this Rstream class")
                         stream
                 } )


## Printing et al. ..........................................................

## print:
##    print all data of a Rstream object
setMethod( "print", "rstream",
          function(x, ...) { .rstream.PrintData(x) } )

setMethod( "show", "rstream",
          function(object) { print(object) } )


## Rstream objects <-> R generators -----------------------------------------

## rstream.RNG
##    get and set Rstream object from/to R global generator.
##    it returns a clone of the current/old global Rstream.
##
##    the current stream is stored in variable '.rstream.current'.

rstream.RNG <- function (stream=NULL) {
	## if there is current stream yet, create one.
	if (!exists(".rstream.current", envir=.rstream.envir) ) {
		s <- new("rstream.runif")
		assign(".rstream.current", s, envir=.rstream.envir) }

	## make a new object of the current Rstream object
	current <- .rstream.getRNG(get(".rstream.current", envir=.rstream.envir))
        
	## replace current Rstream object by given one (if provided)
	if (!is.null(stream)) {
		## input must be Rstream object
		if (!is(stream,"rstream")) stop ("object is not an 'rstream' object")

		## we use a clone of the stream object
		Rs <- rstream.clone(stream)

		## set Rstream object as R generator
		Rs <- .rstream.setRNG(Rs)

		## store current Rstream object
		assign(".rstream.current", Rs, envir=.rstream.envir) }

	## return clone of current or old Rstream object
	if (is.null(stream)) return(current) else return(invisible(current))
}


## .rstream.getRNG
##    get Rstream object for current R generator
##    (internal method; not exported)
if(!isGeneric(".rstream.getRNG"))
	setGeneric(".rstream.getRNG", function(stream) standardGeneric(".rstream.getRNG"))

setMethod(".rstream.getRNG", "rstream", 
          function(stream) {
                  stop(class(stream)[1],": Cannot get object for this Rstream class") } )


## .rstream.setRNG
##    set R generator to given Rstream object
##    (internal method; not exported)
if(!isGeneric(".rstream.setRNG"))
	setGeneric(".rstream.setRNG", function(stream) standardGeneric(".rstream.setRNG"))
setMethod(".rstream.setRNG", "rstream", 
          function(stream) {
                  stop(class(stream)[1],": Cannot set object of this Rstream class as R global generator")
                  stream
          } )


## Auxiliary functions ------------------------------------------------------

## print data 
.rstream.PrintData <- function (stream, ...) {
	name <- rstream.name(stream)
	anti <- rstream.antithetic(stream)
	incp <- rstream.incprecision(stream)
	cat(	"\nObject is random stream\n\n")
	if (stream@is.packed)
		cat("\tRstream is PACKED\n\n")
	cat("\ttype = ", class(stream)[1], "\n",
            "\t       [", stream@type, ": ", stream@info, "]\n", 
            "\tname = \"", name, "\"\n",
            "\tantithetic = ", anti, "\n",
            "\tincreased Precision = ", incp, "\n",
            sep="") 
} 

## Compatibility function to set behavior to earlier version of Rstream -----

## Set version:
##   The version must be given as string.
##   The string is then converted into an number and stored in
##   '.rstream.version'

rstream.version <- function(version) {

  if (missing(version)) {
    num <- get(".rstream.version", envir=.rstream.envir)
    if (num < 9999000) {
      return (paste(num %/% 1000, num %% 1000, sep=".") )
    } else {
      return ("default")  ## default: installed version
    }

  } else {

    if (version=="default") version <- "9999.00"

    num <- as.numeric(strsplit(as.character(version),".", fixed=TRUE)[[1]])

    if ((length(num)!=2 && length(num)!=3)
        || !isTRUE(num[1]>=0) || !isTRUE(num[2]>=0) || !isTRUE(num[2]<=999)) 
      stop("malformed version string")
    
    assign(".rstream.version", 1000*num[1] + num[2], envir=.rstream.envir)
  }
}


## End ----------------------------------------------------------------------
