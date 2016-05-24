
#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Virtual Class: unuran.dist                                            ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.distr", 
         ## slots:
         representation = representation( 
                 distr = "externalptr",   # pointer to UNU.RAN distribution object
                 env   = "environment",   # environment for evaluating R expressions
                 name  = "character"      # name of distribution
                 ),
         ## defaults for slots
         prototype = list(
                 distr = NULL,
                 env   = NULL,
                 name  = NULL
                 ),
         ## indicate virtual class
         contains = "VIRTUAL"
         )


## Printing -----------------------------------------------------------------

## print strings of UNU.RAN object
setMethod( "print", "unuran.distr",
          function(x, ...) {
                  cat("\nObject is UNU.RAN distribution object\n\n")
                  if (!is.null(x@name)) {
                    cat("name:",x@name,"\n\n")
                  }
          } )

setMethod( "show", "unuran.distr",
          function(object) { print(object) } )

## End ----------------------------------------------------------------------
