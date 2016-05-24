#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran.cmv                                                     ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.cmv", 
         ## add slots for continuous multivariate distributions
         representation = representation(
                 ndim = "integer",    # dimensions of distribution
                 pdf  = "function"    # PDF of distribution
                 ),
         ## defaults for slots
         prototype = list(
                 ndim = as.integer(1),
                 pdf  = NULL
                 ),
         ## superclass
         contains = "unuran.distr",
         ## seal this class
         sealed = TRUE )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran.cmv",
          function(.Object, dim=1, pdf=NULL, ll=NULL, ur=NULL, mode=NULL, center=NULL,
                   name=NA, empty=FALSE) {
            ## dim  ... dimension of distribution
            ## pdf  ... probability density function (PDF)
            ## ll   ... lower left vertex of rectangular domain
            ## ur   ... upper right vertex of rectangular domain
            ## mode ... mode of distribution
            ## name ... name of distribution
            ## empty .. if TRUE only return empty object (for internal use only)

            if (isTRUE(empty)) return (.Object)
            
            ## Check entries
            if (! is.numeric(dim))
              stop("invalid argument 'dim'", call.=FALSE)
            ndim <- as.integer(dim)
            if (!isTRUE(all.equal(ndim,dim)) || ndim < 1 || ndim > 100000)
              stop("invalid argument 'dim'", call.=FALSE)

            if(! (is.function(pdf) || is.null(pdf)) )
              stop("invalid argument 'pdf'", call.=FALSE)

            if(! (is.numeric(ll) || is.null(ll)) )
              stop("invalid argument 'll'", call.=FALSE)
            if( (! is.null(ll)) && length(ll)!=ndim ) 
              stop("argument 'll' must have length 'dim'", call.=FALSE)

            if(! (is.numeric(ur) || is.null(ur)) )
              stop("invalid argument 'ur'", call.=FALSE)
            if( (! is.null(ur)) && length(ur)!=ndim ) 
              stop("argument 'ur' must have length 'dim'", call.=FALSE)

            if(! (is.numeric(mode) || is.null(mode)) )
              stop("invalid argument 'mode'", call.=FALSE)
            if( (! is.null(mode)) && length(mode)!=ndim ) 
              stop("argument 'mode' must have length 'dim'", call.=FALSE)

            if(! (is.numeric(center) || is.null(center)) )
              stop("invalid argument 'center'", call.=FALSE)
            if( (! is.null(center)) && length(center)!=ndim ) 
              stop("argument 'center' must have length 'dim'", call.=FALSE)

            if(! (is.character(name) || is.na(name)) )
              stop("invalid argument 'name'", call.=FALSE)
                  
            ## Store informations (if provided)
            .Object@ndim <- ndim
            if (is.function(pdf))  .Object@pdf <- pdf
            if (!is.na(name))      .Object@name <- name
            
            ## We need an evironment for evaluating R expressions
            .Object@env <- new.env()
            
            ## Create UNUR_DISTR object
            .Object@distr <-.Call("Runuran_cmv_init",
                                  .Object, .Object@env,
                                  .Object@ndim, .Object@pdf, mode, center, ll, ur, name,
                                  PACKAGE="Runuran")
            
            ## Check UNU.RAN object
            if (is.null(.Object@distr)) {
              stop("Cannot create UNU.RAN distribution object", call.=FALSE)
            }
            
            ## return new UNU.RAN object
            .Object
          } )


## Shortcut
unuran.cmv.new <- function(dim=1, pdf=NULL, ll=NULL, ur=NULL, mode=NULL, center=NULL, name=NA) {
        new("unuran.cmv", dim=dim, pdf=pdf, ll=ll, ur=ur, mode=mode, center=center, name=name)
}

## End ----------------------------------------------------------------------
