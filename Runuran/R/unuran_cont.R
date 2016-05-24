#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran.cont                                                    ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.cont", 
         ## add slots for continuous univariate distributions
         representation = representation(
                 cdf  = "function",    # CDF of distribution
                 pdf  = "function",    # PDF of distribution
                 dpdf = "function"     # derivative of PDF of distribution
                 ),
         ## defaults for slots
         prototype = list(
                 cdf  = NULL,
                 pdf  = NULL,
                 dpdf = NULL
                 ),
         ## superclass
         contains = "unuran.distr",
         ## seal this class
         sealed = TRUE )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran.cont",
          function(.Object, cdf=NULL, pdf=NULL, dpdf=NULL, islog=FALSE,
                   lb=NA, ub=NA, mode=NA, center=NA, area=NA, name=NA, empty=FALSE) {
            ## cdf .... cumulative distribution function (CDF)
            ## pdf .... probability density function (PDF)
            ## dpdf ... derivative of PDF
            ## islog .. whether CDF and PDF are given as logarithms
            ##          (the dpdf is then the derative of log(pdf)!)
            ## lb ..... lower bound of domain
            ## ub ..... upper bound of domain
            ## mode ... mode of distribution
            ## center . "center" (typical point) of distribution
            ## area ... area below PDF
            ## name ... name of distribution
            ## empty .. if TRUE only return empty object (for internal use only)
            
            if (isTRUE(empty)) return (.Object)
            
            ## Check entries
            if(! (is.numeric(lb) && is.numeric(ub) && lb < ub) )
              stop("domain ('lb','ub') missing or invalid", call.=FALSE)
            
            if(! (is.function(cdf) || is.null(cdf)) )
              stop("invalid argument 'cdf'", call.=FALSE)
            if(! (is.function(pdf) || is.null(pdf)) )
              stop("invalid argument 'pdf'", call.=FALSE)
            if(! (is.function(dpdf) || is.null(dpdf)) )
              stop("invalid argument 'dpdf'", call.=FALSE)
            
            if(! is.logical(islog))
              stop("argument 'islog' must be boolean", call.=FALSE)
            
            if(! (is.numeric(mode) || is.na(mode)) )
              stop("invalid argument 'mode'", call.=FALSE)
                  
            if(! (is.numeric(center) || is.na(center)) )
              stop("invalid argument 'center'", call.=FALSE)
                  
            if(! is.na(area) && !(is.numeric(area) && area>0) )
              stop("invalid argument 'area'", call.=FALSE)

            if(! (is.character(name) || is.na(name)) )
              stop("invalid argument 'name'", call.=FALSE)
                  
            ## Store informations (if provided)
            if (is.function(cdf))  .Object@cdf  <- cdf
            if (is.function(pdf))  .Object@pdf  <- pdf
            if (is.function(dpdf)) .Object@dpdf <- dpdf
            if (!is.na(name))      .Object@name <- name
            
            ## We need an evironment for evaluating R expressions
            .Object@env <- new.env()
            
            ## Create UNUR_DISTR object
            .Object@distr <-.Call("Runuran_cont_init",
                                  .Object, .Object@env,
                                  .Object@cdf, .Object@pdf, .Object@dpdf, islog,
                                  mode, center, c(lb,ub), area, name,
                                  PACKAGE="Runuran")
            
            ## Check UNU.RAN object
            if (is.null(.Object@distr)) {
              stop("Cannot create UNU.RAN distribution object", call.=FALSE)
            }
            
            ## return new UNU.RAN object
            .Object
          } )

## Shortcut
unuran.cont.new <- function(cdf=NULL, pdf=NULL, dpdf=NULL, islog=FALSE,
                            lb=NA, ub=NA, mode=NA, center=NA, area=NA, name=NA) {
  new("unuran.cont", cdf=cdf, pdf=pdf, dpdf=dpdf, islog=islog,
      lb=lb, ub=ub, mode=mode, center=center, area=area, name=name)
}

## End ----------------------------------------------------------------------
