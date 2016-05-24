#############################################################################
##                                                                         ##
##   Class: rstream.lecuyer                                                ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Interface to Pierre L'Ecuyers RngStreams Library                      ##
##                                                                         ##
#############################################################################


## Initialize global variables ----------------------------------------------

.rstream.lecuyer.init <- function () {

	## Indicate whether a package seed is already set
	if (!exists(".rstream.lecuyer.HasSeed", envir=.rstream.envir))
		assign(".rstream.lecuyer.HasSeed", FALSE, envir=.rstream.envir)

	## Default Package seed
	if (!exists(".rstream.lecuyer.DefaultSeed", envir=.rstream.envir))
		assign(".rstream.lecuyer.DefaultSeed", rep(as.double(Sys.time()),6), envir=.rstream.envir)
}


## Class --------------------------------------------------------------------

setClass( "rstream.lecuyer", representation(), contains = "rstream" )


## Initialize ---------------------------------------------------------------

setMethod( "initialize", "rstream.lecuyer",  
          function(.Object, name=NULL, seed=NULL, force.seed=FALSE, 
                   antithetic=FALSE, incprecision=FALSE) {

                  ## first increment counter for Rstreams
                  count <- get(".rstream.Count", envir=.rstream.envir) + 1
                  assign(".rstream.Count", count, envir=.rstream.envir)

                  ## name of the Rstream object.
                  ## by default we use type + number
                  if (is.null(name)) name <- paste("lecuyer", count, sep="")

                  ## type of Rstream object
                  .Object@type <- "lecuyer"

                  ## add info about Rstream type
                  .Object@info <- "RngStreams - library for multiple independent streams of Random Numbers"

                  ## at creation a Rstream is never packed
                  .Object@is.packed <- FALSE 

                  ## set (package) seed
                  hasSeed <- get(".rstream.lecuyer.HasSeed", envir=.rstream.envir)
                  defaultSeed <- get(".rstream.lecuyer.DefaultSeed", envir=.rstream.envir)
                  if (!is.null(seed) && !force.seed && hasSeed)
                          stop("rstream.lecuyer: Package already seeded! Seed ignored!\n",
                               "rstream.lecuyer: Set force.seed = TRUE if you really want to reseed.")
                  if (is.null(seed)) seed <- defaultSeed
                  seed <- .rstream.lecuyer.CheckSeed(seed)
                  .Call("R_RngStreams_SetPackageSeed", as.double(seed), PACKAGE="rstream")

                  ## Create Rstream object
                  .Object@stream <- .Call("R_RngStreams_Init", .Object, as.character(name), PACKAGE="rstream")
                  .Call("R_RngStreams_SetAntithetic", .Object@stream, as.integer(antithetic), PACKAGE="rstream")
                  .Call("R_RngStreams_SetIncreasedPrecis", .Object@stream, as.integer(incprecision), PACKAGE="rstream")

                  ## save seed as R variavble
                  assign(".rstream.lecuyer.DefaultSeed",
                         as.numeric(.Call("R_RngStreams_GetPackageSeed", PACKAGE="rstream")), envir=.rstream.envir)
                  assign(".rstream.lecuyer.HasSeed", TRUE, envir=.rstream.envir)

                  ## return new Rstream object
                  .Object
          } )


## Validity -----------------------------------------------------------------

## .rstream.lecuyer.CheckSeed
##    make simple check on given seed
.rstream.lecuyer.CheckSeed <- function(seed) {
	ll <- length(seed)
	if (ll < 6)
		stop("rstream.lecuyer: seed too short; 6 numbers required")
	if (ll > 6) {
		warning("rstream.lecuyer: seed too long; truncated") 
		seed<-seed[1:6] }
###	for (i in 1:6)
###		if (seed[i] < 0) 
###			stop("invalid seed")
	seed
}


## Methods ------------------------------------------------------------------

## Access and Replacement methods ...........................................

## rstream.name
##    get and set name of Rstream object
setMethod("rstream.name", "rstream.lecuyer", 
          function(stream) { 
                  if (stream@is.packed) 
                          return (stream@pack$name)
                  else
                          return (.Call("R_RngStreams_GetName", stream@stream, PACKAGE="rstream")) 
          } )

setReplaceMethod("rstream.name", "rstream.lecuyer", 
                 function(stream, value) {
                         if (stream@is.packed) stop("Cannot change name for PACKED Rstream") 
                         .Call("R_RngStreams_SetName", stream@stream, as.character(value), PACKAGE="rstream")
                         stream
                 } )


## rstream.antithetic
##   get and set flag for antithetic random numbers:  
setMethod("rstream.antithetic", "rstream.lecuyer", 
          function(stream) {
                  if (stream@is.packed) 
                          return (stream@pack$anti)
                  else 
                          return (as.logical(.Call("R_RngStreams_GetAntithetic", stream@stream, PACKAGE="rstream")))
          } )

setReplaceMethod("rstream.antithetic", "rstream.lecuyer",
                 function(stream, value) { 
                         if (stream@is.packed) stop("Cannot change antithetic flag for PACKED Rstream") 
                         .Call("R_RngStreams_SetAntithetic", stream@stream, as.integer(value), PACKAGE="rstream")
                         stream
                 } )


## rstream.incprecision
##    get and set flag for increased precision of random numbers
setMethod("rstream.incprecision", "rstream.lecuyer", 
          function(stream) {
                  if (stream@is.packed) 
                          return (stream@pack$incp)
                  else 
                          return (as.logical(.Call("R_RngStreams_GetIncreasedPrecis", stream@stream, PACKAGE="rstream")))
          } )

setReplaceMethod("rstream.incprecision", "rstream.lecuyer",
                 function(stream, value) { 
                         if (stream@is.packed) stop("Cannot change increased precision flag for PACKED Rstream") 
                         .Call("R_RngStreams_SetIncreasedPrecis", stream@stream, as.integer(value), PACKAGE="rstream")
                         stream
                 } )


## Sampling methods .........................................................

## rstream.sample
##    make a random sample
setMethod("rstream.sample", "rstream.lecuyer",
          function(stream,n=1) { 
                  if (stream@is.packed) stop("Cannot sample from PACKED Rstream") 
                  .Call("R_RngStreams_Sample", stream@stream, as.integer(n), PACKAGE="rstream") } )

setMethod("r", "rstream.lecuyer",
          function(stream,n=1) { 
                  if (stream@is.packed) stop("Cannot sample from PACKED Rstream") 
                  .Call("R_RngStreams_Sample", stream@stream, as.integer(n), PACKAGE="rstream") } )


## Jump methods .............................................................

## rstream.resetsubstream
##   reset current substream of Rstream object
if(!isGeneric("rstream.resetsubstream"))
        setGeneric("rstream.resetsubstream", function(stream) standardGeneric("rstream.resetsubstream"))
setMethod("rstream.resetsubstream", "rstream.lecuyer", 
          function(stream) { 
                  if (stream@is.packed) stop("Cannot reset PACKED Rstream") 
                  dummy <- .Call("R_RngStreams_ResetStartSubstream", stream@stream, PACKAGE="rstream")
          } )


## rstream.nextsubstream
##   skip to next substream of Rstream object
if(!isGeneric("rstream.nextsubstream"))
        setGeneric("rstream.nextsubstream", function(stream) standardGeneric("rstream.nextsubstream"))
setMethod("rstream.nextsubstream", "rstream.lecuyer", 
          function(stream) { 
                  if (stream@is.packed) stop("Cannot skip substream of PACKED Rstream") 
                  dummy <- .Call("R_RngStreams_ResetNextSubstream", stream@stream, PACKAGE="rstream")
          } )


## Reset and copy methods ...................................................

## rstream.reset
##   reset Rstream object
setMethod("rstream.reset", "rstream.lecuyer", 
          function(stream) { 
                  if (stream@is.packed) stop("Cannot reset PACKED Rstream") 
                  dummy <- .Call("R_RngStreams_ResetStartStream", stream@stream, PACKAGE="rstream")
          } )


## rstream.clone
##    clone (copy) Rstream object
setMethod("rstream.clone", "rstream.lecuyer", 
          function(stream) { 
                  if (stream@is.packed) stop("Cannot clone PACKED Rstream") 
                  clone <- stream
                  name <- paste(rstream.name(stream),".",sep="")
                  clone@stream <- .Call("R_RngStreams_Clone", clone, stream@stream, name, PACKAGE="rstream")
                  clone
          } )


## rstream.pack, rstream.unpack
##    pack and unpack Rstream object such that all data are contained in R object
##    (and can be easily copied within R)
setReplaceMethod("rstream.packed", "rstream.lecuyer", 
                 function(stream, value) {
                         value <- as.logical(value)
                         if (value && stream@is.packed)   return(stream)
                         if (!value && !stream@is.packed) return(stream)
                         if (value) {	# pack
                                 name <- rstream.name(stream)
                                 anti <- rstream.antithetic(stream)
                                 incp <- rstream.incprecision(stream)
                                 stream@is.packed <- TRUE
                                 stream@pack <- list()
                                 stream@pack$state <- as.double(.Call("R_RngStreams_GetData", stream@stream, PACKAGE="rstream"))
                                 stream@pack$name <- name 
                                 stream@pack$anti <- anti
                                 stream@pack$incp <- incp
                                 .Call("R_RngStreams_Free", stream@stream, PACKAGE="rstream")
                         }
                         else {		# unpack
                                 stream@is.packed <- FALSE
                                 .Call("R_RngStreams_SetData", stream,
                                       stream@stream, stream@pack$state,
                                       stream@pack$name, PACKAGE="rstream")
                         }
                         stream
                 } )


## Printing et al. ..........................................................

## print:
##    print all data of a Rstream object
setMethod( "print", "rstream.lecuyer",
          function(x, ...) { 
                  .rstream.PrintData(x) 
                  if (!x@is.packed) {
                          state <- .Call("R_RngStreams_GetData", x@stream, PACKAGE="rstream")
                          cat("\n\tInternal state:\n",
                              "\t  initial state:\n",
                              "\t\t", state[13:18], "\n",
                              "\t  starting point of current substream:\n",
                              "\t\t", state[7:12], "\n",
                              "\t  current state:\n",
                              "\t\t", state[1:6], "\n\n" ) }
                  else
                          cat("\n")
          } )


## Rstream objects <-> R generators -----------------------------------------

## .rstream.getRNG
##    get Rstream object for current R generator
##    (internal method; not exported)
setMethod(".rstream.getRNG", "rstream.lecuyer", 
          function(stream) { rstream.clone(stream) } )


## .rstream.setRNG
##    set R generator to given Rstream object
##    (internal method; not exported)
setMethod(".rstream.setRNG", "rstream.lecuyer", 
          function(stream) {
            if (isTRUE(get(".rstream.version", envir=.rstream.envir) < 1003)) {
              ## pre 1.3:
              .Call("R_RngStreams_setRNG", stream@stream, PACKAGE="rstream")
              RNGkind(kind="user-supplied")
            } else {
              ## new version:
              RNGkind(kind="user-supplied")
              .Call("R_RngStreams_setRNG", stream@stream, PACKAGE="rstream")
            }
            stream
          } )

## End ----------------------------------------------------------------------
