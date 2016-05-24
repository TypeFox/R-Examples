#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007-2010, Josef Leydold and Wolfgang Hoermann                    ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran                                                         ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran", 
         ## slots:
         representation( 
                        unur       = "externalptr",   # pointer to UNU.RAN object
                        data       = "list",          # list for packed data
                        dom        = "numeric",       # domain of distribution
                        distr      = "unuran.distr",  # pointer to S4 distribution object
                        distr.str  = "character",     # distribution
                        method.str = "character",     # generation method
                        inversion  = "logical"        # inversion method ?
                        ),
         ## defaults for slots
         prototype = list(
           unur       = NULL,
           data       = NULL,
           dom        = NULL,
           distr      = NULL,
           distr.str  = character(),
           method.str = character(),
           inversion  = FALSE
           ),
         ## misc
         sealed = TRUE
         )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran",  
          function(.Object, distr, method="auto") {

                  ## Check entries
                  if (missing(distr)) {
                          stop("no distribution given", call.=FALSE)
                  }
                  if (!is.character(method)) {
                          stop("'method' must be a character string", call.=FALSE)
                  }

                  ## Create an empty object if distr equals NULL
                  if (is.null(distr)) {
                          return(.Object)
                  }
                  
                  ## Store informations 
                  .Object@distr.str <- ifelse(is.character(distr), distr, "[S4 class]")
                  .Object@method.str <- method

                  ## Create UNU.RAN object
                  if (is.character(distr)) {
                          .Object@unur <-.Call("Runuran_init", .Object, distr, method, PACKAGE="Runuran")
                  } else { if (class(distr)=="unuran.discr" ||
                               class(distr)=="unuran.cont"  ||
                               class(distr)=="unuran.cmv") {
                          .Object@unur <-.Call("Runuran_init", .Object, distr@distr, method, PACKAGE="Runuran")
                          .Object@distr <- distr
                  } else {
                          stop("'distr' must be a character string or a Runuran distribution object", call.=FALSE)
                  }}

                  ## Check UNU.RAN object
                  if (is.null(.Object@unur)) {
                          stop("Cannot create UNU.RAN object", call.=FALSE)
                  }
                  
                  ## return new UNU.RAN object
                  .Object
          } )

## Shortcut
unuran.new <- function(distr,method="auto") {
        new("unuran",distr,method)
}

## Validity -----------------------------------------------------------------

## Sampling -----------------------------------------------------------------

## ur
## ( We avoid using a method as this has an expensive overhead. )
ur <- function(unr,n=1) { 
        .Call("Runuran_sample", unr, n, PACKAGE="Runuran")
}

## unuran.sample: deprecated name for ur()
unuran.sample <- function(unr,n=1) { 
        .Call("Runuran_sample", unr, n, PACKAGE="Runuran")
}

## Quantile -----------------------------------------------------------------

## uq
uq <- function(unr,U) { 
        .Call("Runuran_quantile", unr, U, PACKAGE="Runuran")
}

## PDF & PMF ----------------------------------------------------------------

## ud
ud <- function(obj,x,islog=FALSE) {
  if ( ! (is(obj,"unuran.cont") || is(obj,"unuran.discr") ||
          is(obj,"unuran") ) )
    stop("argument 'obj' must be UNU.RAN object")
  .Call("Runuran_PDF", obj, x, islog, PACKAGE="Runuran")
}

## CDF ----------------------------------------------------------------------

## up
up <- function(obj,x) {
  if ( ! (is(obj,"unuran.cont") || is(obj,"unuran.discr") ||
          is(obj,"unuran") ) )
    stop("argument 'obj' must be UNU.RAN object")
  .Call("Runuran_CDF", obj, x, PACKAGE="Runuran")
}

## Packing ------------------------------------------------------------------

## We have the following two situtations for
## slots 'unur', 'data', and 'dom':
##
## 1. Runuran object is NOT PACKED:
##    'unur' ... contains pointer to UNU.RAN object
##    'data' ... is set to NULL
##    'dom'  ... is ignored
##
## 2. Runuran object is PACKED:
##    'unur' ... contains NULL pointer '(nil)'
##    'data' ... contains list of coefficients
##    'dom'  ... contains domain of distribution
##
## Otherwise, if 'unur' points to '(nil)' and 'data' is still set to NULL
## then the Runuran object is assumed to be 'empty'.
## (This happens if an unuran object is created with distr=NULL.)
##
## It should not happen that both 'unur' and 'data' contain non-NULL objects.
##
## ..........................................................................

## unuran.packed
##    Some Runuran object can be packed such that all data are stored as 
##    R object (and thus can be copied and saved within R)
if(!isGeneric("unuran.packed"))
        setGeneric("unuran.packed", function(unr) standardGeneric("unuran.packed"))

setMethod("unuran.packed", "unuran", 
          function(unr) {
            if (is.null(unr@data)) return(FALSE)
            return(TRUE)
          } )

if(!isGeneric("unuran.packed<-"))
        setGeneric("unuran.packed<-", function(unr, value) standardGeneric("unuran.packed<-"))

setReplaceMethod("unuran.packed", "unuran", 
                 function(unr, value) {
                   value <- as.logical(value)
                   is.packed <- !is.null(unr@data)

                   if (value && is.packed) {
                     warning("[UNU.RAN - warning] object already PACKED", call.=FALSE)
                     
                   }
                   if (!value && is.packed) {
                     ## we cannot unpack object data
                     stop("[UNU.RAN - error] Cannot unpack 'unuran' object", call.=FALSE)
                   }
                   if (value && !is.packed) {
                     ## pack data
                     .Call("Runuran_pack", unr, PACKAGE="Runuran")
                   }
                   ## otherwise: nothing to do
                   
                   return (unr)
                 } )


## Second (auxiliary) URNG  -------------------------------------------------

if(!isGeneric("use.aux.urng"))
  setGeneric("use.aux.urng", function(unr) standardGeneric("use.aux.urng"))

setMethod("use.aux.urng", "unuran", 
          function(unr) {
                  .Call("Runuran_use_aux_urng", unr, NULL, PACKAGE="Runuran")
          } )

          
if(!isGeneric("use.aux.urng<-"))
        setGeneric("use.aux.urng<-", function(unr, value) standardGeneric("use.aux.urng<-"))
          
setReplaceMethod("use.aux.urng", "unuran", 
                 function(unr, value) {
                         value <- as.logical(value)
                         .Call("Runuran_use_aux_urng", unr, value, PACKAGE="Runuran")
                         return (unr)
                 } )

set.aux.seed <- function(seed) {
        seed <- as.integer(seed)
        if (seed <= 0) stop("seed must be positive integer");
        invisible(.Call("Runuran_set_aux_seed", seed, PACKAGE="Runuran"))
}


## Printing -----------------------------------------------------------------

## print strings of UNU.RAN object
setMethod( "print", "unuran",
          function(x, ...) {
            cat("\nObject is UNU.RAN object:\n")
            cat("\tmethod:   ",x@method.str,"\n")
            cat("\tdistr:    ",x@distr.str,"\n")
            cat("\tinversion:",x@inversion,"\n\n")
            cat(.Call("Runuran_print", x, FALSE, PACKAGE="Runuran"))
            cat("")
} )

setMethod( "show", "unuran",
          function(object) { print(object) } )

## unuran.details
## (print for information and hints)
unuran.details <- function(unr, show=TRUE, return.list=FALSE, debug=FALSE) {
  if (isTRUE(show)) {
    cat("\nObject is UNU.RAN object:\n")
    cat("\tmethod:   ",unr@method.str,"\n")
    cat("\tdistr:    ",unr@distr.str,"\n")
    cat("\tinversion:",unr@inversion,"\n\n")
    info <- .Call("Runuran_print", unr, TRUE, PACKAGE="Runuran")
    cat(info)
  }
  if (isTRUE(return.list) || isTRUE(debug)) {
    data <- .Call("Runuran_performance", unr, debug, PACKAGE="Runuran")
    invisible(data)
  }
}


## Mixture ------------------------------------------------------------------

## UNU.RAN meta method for sampling from a mixture of distributions

mixt.new <- function (prob, comp, inversion=FALSE) {

  ## Check arguments
  if (length(prob) != length(comp))
    stop ("'prob' and 'comp' must have same length")
  if (! is.numeric(prob))
    stop ("invalid argment 'prob'")
  if (! is.list(comp))
    stop ("invalid argment 'comp'")

  if (! is(comp[[1]],"unuran")) {
    if (is(comp[[1]],"unuran.distr")) {
      stop ("invalid argment: 'comp' must be list of generators not of distributions")
    } else {
      stop ("invalid argment: 'comp' must be list of generators")
    }
  }
  
  ## Create empty "unuran" object.
  obj <- new("unuran",distr=NULL)
             
  ## Store informations 
  obj@distr.str <- "mixture of distributions"
  obj@method.str <- "mixt"

  ## Create UNU.RAN object
  obj@unur <- .Call("Runuran_mixt", obj, prob, comp, inversion, PACKAGE="Runuran")
  if (is.null(obj@unur)) {
    stop("Cannot create UNU.RAN object", call.=FALSE)
  }

  ## Return new UNU.RAN object
  obj
}


## Check hat function -------------------------------------------------------

## verify hat and squeeze of a rejection method.
## it counts the number of violations, i.e.,
## when (hat(x) < density(x) or squeeze(x) > density(x))

unuran.verify.hat <- function (unr, n=1e5, show = TRUE) {

  ## check arguments
  if (! is(unr,"unuran"))
    stop ("invalid argument 'unr'");
  if (! (is.numeric(n) && n>10) )
    stop ("invalid sample size 'n'");
  
  ## run test
  failed <- .Call("Runuran_verify_hat", unr, n, PACKAGE="Runuran")
  ratio <- failed / n
  perc <- round(100*ratio,digits=2)

  ## print diagnostics
  if (isTRUE(show)) {
    cat ("\n  Check inequality squeeze(x) <= density(x) <= hat(x) \n  for automatic rejection method:\n\n")
    cat ("\t ",failed," out of ",n," (= ",perc,"%) points failed.\n\n", sep="")

    if (isTRUE(all.equal(failed,0))) {
      cat ("\t No problems have been detected!\n\n")
    }
    else if (ratio < 1e-4) {
      cat ("\t Some points failed!\n")
      cat ("\t Probably there are a few round-off errors.\n\n")
    }
    else if (ratio < 1e-3) {
      cat ("\t Some points failed!\n")
      cat ("\t These might be round-off errors.\n\n")
    }
    else if (ratio < 1e-2) {
      cat ("\t Points failed!\n")
      cat ("\t These might be round-off errors but the\n")
      cat ("\t sampling distribution could deviate from the requested one.\n\n")
    }
    else {
      cat ("\t Verifying hat and squeeze failed!\n")
      cat ("\t The generator does not sample from the \n")
      cat ("\t requested distribution!\n\n")
    }
  }

  ## return ratio
  invisible(ratio)
}


## Test for inversion method ------------------------------------------------

## Test whether given unuran object implements an (approximate) inversion
## method.

unuran.is.inversion <- function (unr) {

  ## check arguments
  if (! is(unr,"unuran"))
    stop ("invalid argument 'unr'");

  ## return result
  unr@inversion
}

## End ----------------------------------------------------------------------
