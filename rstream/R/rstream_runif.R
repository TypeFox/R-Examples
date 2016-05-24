#############################################################################
##                                                                         ##
##   Class: rstream.runif                                                  ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Interface to R uniform random number generators                       ##
##                                                                         ##
#############################################################################



## Initialize global variables ----------------------------------------------

.rstream.runif.init <- function () { }
## nothing to do


# Class ---------------------------------------------------------------------

setClass( "rstream.runif", 
         ## additional slots	
         representation(
                        kind	= "character",	  # kind of RNG
                        seed	= "integer" ),	  # seed of generator
         ## extends rstream virtual class
         contains = "rstream" )


## Initialize ---------------------------------------------------------------

setMethod( "initialize", "rstream.runif",  
          function(.Object, name=NULL, kind="current", seed=NULL, antithetic=FALSE) {

                  ## check input
                  if (as.logical(pmatch(kind,"user-supplied",nomatch=FALSE)))
                          stop("kind=\"user-supplied\" not allowed for Rstream object") 

                  ## first increment counter for Rstreams
                  count <- get(".rstream.Count", envir=.rstream.envir) + 1
                  assign(".rstream.Count", count, envir=.rstream.envir)

                  ## Create an environment for storing state of the generator
                  .Object@xstate <- new.env(parent=.rstream.envir)

                  ## name of the Rstream object.
                  ## by default we use type + number
                  if (is.null(name)) name <- paste("runif", count, sep="")
                  assign("name", as.character(name), envir=.Object@xstate)

                  ## type of Rstream object
                  .Object@type <- "runif"

                  ## add info about Rstream type
                  .Object@info <- "R uniform random number generator"

                  ## at creation a Rstream is never packed
                  .Object@is.packed <- FALSE 

                  ## set antithetic flag
                  assign("anti", as.logical(antithetic), envir=.Object@xstate)

                  ## save current state of runif.
                  ## however, .Random.seed might not exist.
                  ## then we run runif() one time to initialize it.
                  if (!exists(".Random.seed", envir=.GlobalEnv)) runif(1)
                  current.seed <- get(".Random.seed", envir=.GlobalEnv)

                  ## we have to run runif() once to get the correct answer
                  ## from RNGkind()
                  ## (see also the comments in rstream.sample)
                  dummy <- runif(1)
                  current.kind <- RNGkind()[1]

                  ## we cannot use the "current" RNG if it is a user-supplied.
                  ## then we use the "default" one
                  if (kind == "current" &&
                      as.logical(pmatch(current.kind,"user-supplied",nomatch=FALSE))) {
                          warning("kind=\"user-supplied\" not allowed for Rstream object; using \"default\".") 
                          kind <- "default"
                  }

                  ## "current" RNG
                  ## we only have to set the seed (if one is provided by user)
                  if (kind == "current") {
                          if(!is.null(seed)) {
                                  set.seed(seed) 
                                  .Object@seed <- get(".Random.seed", envir=.GlobalEnv) }
                          else
                                  .Object@seed <- current.seed
                  }

                  ## other RNG
                  ## we use RNGkind() to check the input
                  else {
                          RNGkind(kind=kind)
                          if (is.null(seed)) seed <- Sys.time()
                          set.seed(seed) 
                          .Object@seed <- get(".Random.seed", envir=.GlobalEnv) }
                  
                  ## store kind of stream
                  .Object@kind <- RNGkind()[1]
                  
                  ## store state of generator in its own environment
                  assign("state", .Object@seed, envir=.Object@xstate)
                  
                  ## restore global RNGkind
                  RNGkind(kind=current.kind)
                  assign(".Random.seed", current.seed, envir=.GlobalEnv)
                  
                  ## return new Rstream object
                  .Object
          } )


## Validity -----------------------------------------------------------------


## Methods ------------------------------------------------------------------

## Access and Replacement methods ............................................

## rstream.name
##    get and set name of Rstream object
setMethod("rstream.name", "rstream.runif", 
          function(stream) { 
                  if (stream@is.packed) 
                          return (stream@pack$name)
                  else
                          return (get("name", envir=stream@xstate))
          } )

setReplaceMethod("rstream.name", "rstream.runif", 
                 function(stream, value) {
                         if (stream@is.packed) stop("Cannot change name for PACKED Rstream") 
                         assign("name", as.character(value), envir=stream@xstate)
                         stream
                 } )


## rstream.antithetic
##   get and set flag for antithetic random numbers:  
setMethod("rstream.antithetic", "rstream.runif", 
          function(stream) { 
                  if (stream@is.packed) 
                          return (stream@pack$anti)
                  else
                          return (get("anti", envir=stream@xstate))
          } )

setReplaceMethod("rstream.antithetic", "rstream.runif",
                 function(stream, value) { 
                         if (stream@is.packed) stop("Cannot change antithetic flag for PACKED Rstream") 
                         assign("anti", as.logical(value), envir=stream@xstate)
                         stream
                 } )


## Sampling methods .........................................................

## rstream.sample
##    make a random sample
setMethod("rstream.sample", "rstream.runif",
          function(stream,n=1) .rstream.runif.sample(stream,n) )

setMethod("r", "rstream.runif",
          function(stream,n=1) .rstream.runif.sample(stream,n) )

.rstream.runif.sample <- 
	function(stream,n=1) { 
		if (stream@is.packed) stop("Cannot sample from PACKED Rstream") 
		## store current state of R generator
		old.seed <- get(".Random.seed", envir=.GlobalEnv)
                
		## read state from rstream table and
		## change R generator
		## Notice: The RNGkind is stored in .Random.seed[1]
		## Thus there is no need to make an extensive call to 
		## RNGkind(kind=kind). However, RNGkind() will result
		## a wrong kind until the first call to runif().
		state <- get("state", envir=stream@xstate)
		assign(".Random.seed", state, envir=.GlobalEnv)
                
		## sample
		x <- runif(n)
                
		## store R generator state in rstream table and
		## reset R generator to original state
		state <- get(".Random.seed", envir=.GlobalEnv)
		assign("state", state, envir=stream@xstate)
		assign(".Random.seed", old.seed, envir=.GlobalEnv)
                
		## return result (depending on the antithetic variable flag)
		if (get("anti", envir=stream@xstate))
			return (1-x)
		else
			return (x)
	}


## Reset and copy methods ...................................................

## rstream.reset
##   reset Rstream object
setMethod("rstream.reset", "rstream.runif", 
          function(stream) { 
                  if (stream@is.packed) stop("Cannot reset PACKED Rstream") 
                  ## copy seed into state variable
                  assign("state", stream@seed, envir=stream@xstate) } )


## rstream.clone
##    clone (copy) Rstream object
setMethod("rstream.clone", "rstream.runif", 
          function(stream) { 
                  if (stream@is.packed) stop("Cannot clone PACKED Rstream") 
                  ## copy R object
                  clone <- stream
                  ## we need a new environment for this clone
                  clone@xstate <- new.env(parent=.rstream.envir)
                  ## copy data from old to new environment
                  ## add a period '.' to name to mark this as clone
                  name <- paste(get("name", envir=stream@xstate),".",sep="")
                  assign("name", name, envir=clone@xstate)
                  state <- get("state", envir=stream@xstate)
                  assign("state", state, envir=clone@xstate)
                  anti <- get("anti", envir=stream@xstate)
                  assign("anti", anti, envir=clone@xstate)
                  ## return clone of Rstream object
                  clone
          } )


## rstream.pack, rstream.unpack
##    pack and unpack Rstream object such that all data are contained in R object
##    (and can be easily copied within R)
setReplaceMethod("rstream.packed", "rstream.runif", 
                 function(stream, value) {
                         value <- as.logical(value)
                         ## check whether there is something to do
                         if (value && stream@is.packed)   return(stream)
                         if (!value && !stream@is.packed) return(stream)
                         
                         if (value) {
                                 ## pack
                                 stream@is.packed <- TRUE
                                 stream@pack <- as.list(stream@xstate)
                         }
                         else {
                                 ## unpack
                                 stream@is.packed <- FALSE
                                 ## we need a new environment for the unpacked stream
                                 stream@xstate <- new.env(parent=.rstream.envir)
			## copy data into new environment
                                 assign("name", stream@pack$name, envir=stream@xstate)
                                 assign("state", stream@pack$state, envir=stream@xstate)
                                 assign("anti", stream@pack$anti, envir=stream@xstate)
                         }
                         ## return un/packed Rstream object
                         stream
                 } )


## Printing et al. ..........................................................

## print:
##    print all data of a Rstream object
setMethod( "print", "rstream.runif",
          function(x, ...) { 
                  .rstream.PrintData(x) 
                  cat(	"\tRNGkind = ", x@kind, "\n\n", sep="")
          } )


## Rstream objects <-> R generators -----------------------------------------

## .rstream.getRNG
##    get Rstream object for current R generator
##    (internal method; not exported)
setMethod(".rstream.getRNG", "rstream.runif", 
          function(stream) { new("rstream.runif", kind="current") } )


## .rstream.setRNG
##    set R generator to given Rstream object
##    (internal method; not exported)

setMethod(".rstream.setRNG", "rstream.runif", 
          function(stream) {
                  if (rstream.antithetic(stream)) {
                          warning ("antithetic DISABLED")
                          rstream.antithetic(stream) <- FALSE }	
                  RNGkind(kind=stream@kind)
                  state <- get("state", envir=stream@xstate)
                  assign(".Random.seed", state, envir=.GlobalEnv)
                  stream
          } )


## End ----------------------------------------------------------------------
