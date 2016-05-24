#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran.discr                                                   ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.discr", 
         ## add slots for discrete distributions
         representation = representation(
                 cdf  = "function",    # CDF of distribution
                 pmf  = "function"     # PMF of distribution
                 ),
         ## defaults for slots
         prototype = list(
                 cdf  = NULL,
                 pmf  = NULL
                 ),
         ## superclass
         contains = "unuran.distr",
         ## seal this class
         sealed = TRUE )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran.discr",
          function(.Object, cdf=NULL, pv=NULL, pmf=NULL, lb=NA, ub=NA,
                   mode=NA, sum=NA, name=NA, empty=FALSE) {
            ## cdf .... cumulative distribution function (CDF)
            ## pv ..... probability vector (PV)
            ## pmf .... probability mass function (PMF)
            ## lb ..... lower bound of domain
            ## ub ..... upper bound of domain
            ## mode ... mode of distribution
            ## sum .... sum over PV / PMF
            ## name ... name of distribution
            ## empty .. if TRUE only return empty object (for internal use only)
            
            if (isTRUE(empty)) return (.Object)
            
            ## Check entries
            if(! (is.function(cdf) || is.null(cdf)) )
              stop("invalid argument 'cdf'", call.=FALSE)

            if (! (is.numeric(pv) || is.null(pv) ))
              stop("invalid argument 'pv'", call.=FALSE)

            if (! (is.function(pmf) || is.null(pmf)) )
              stop("invalid argument 'pmf'", call.=FALSE)

            if (! is.null(pv) && is.na(ub) )
              ## we can use the default 'ub=Inf' if only the PV is given
              ub <- Inf

            ## we need both 'lb' and 'ub'
            if (! (is.numeric(lb) && is.numeric(ub) && lb < ub) )
              stop("domain ('lb','ub') missing or invalid", call.=FALSE)
            
            if(! (is.numeric(mode) || is.na(mode)) )
              stop("invalid argument 'mode'", call.=FALSE)
                  
            if(! is.na(sum) && !(is.numeric(sum) && sum>0) )
              stop("invalid argument 'sum'", call.=FALSE)

            if(! (is.character(name) || is.na(name)) )
              stop("invalid argument 'name'", call.=FALSE)
                  

            ## Store informations (if provided)
            if (is.function(cdf)) .Object@cdf  <- cdf
            if (is.function(pmf)) .Object@pmf  <- pmf
            if (!is.na(name))     .Object@name <- name
            ## (There is no need to store the PV)
            
            ## We need an evironment for evaluating R expressions
            .Object@env <- new.env()
            
            ## Create UNUR_DISTR object
            .Object@distr <-.Call("Runuran_discr_init",
                                  .Object, .Object@env,
                                  .Object@cdf, pv, .Object@pmf,
                                  mode, c(lb,ub), sum, name,
                                  PACKAGE="Runuran")
            
            ## Check UNU.RAN object
            if (is.null(.Object@distr)) {
              stop("Cannot create UNU.RAN distribution object", call.=FALSE)
            }
            
            ## return new UNU.RAN object
            .Object
          } )

## Shortcut
unuran.discr.new <- function(cdf=NULL, pv=NULL, pmf=NULL, lb=NA, ub=NA,
                             mode=NA, sum=NA, name=NA) {
  new("unuran.discr", cdf=cdf, pv=pv, pmf=pmf, lb=lb, ub=ub,
      mode=mode, sum=sum,name=name)
}

## End ----------------------------------------------------------------------
