#==============================================================================
# File: strumMarker.R
#
# Author: Nathan Morris
#
# Notes: strumMarker class definition & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of strumMarker class
#------------------------------------------------------------------------------
setClass("strumMarker",
         representation(
           markerFacts          = "data.frame",
           haplotypes           = "matrix",
           populationRecombRate = "numeric",
           errorRate            = "numeric",
           mutationRate         = "numeric",
           missingRate          = "numeric",
           coding               = "numeric",
           returnIBD            = "logical",
           intervalIBD          = "numeric"),
         prototype = list(
           markerFacts          = data.frame(),
           haplotypes           = matrix(NA), 
           populationRecombRate = 50,
           errorRate            = numeric(0),
           mutationRate         = numeric(0),
           missingRate          = numeric(0),
           coding               = c(0, 1, 2),
           returnIBD            = logical(0),
           intervalIBD          = 5)
)

#------------------------------------------------------------------------------
# Check validity of strumMarker object
#------------------------------------------------------------------------------
setValidity ("strumMarker",
             function(object)
             {
               retval <- NULL
               if( ( object@errorRate < 0 | object@errorRate > 1) )
                 retval <- c( retval , "errorRate is invalid")

               if( ( object@mutationRate < 0 | object@mutationRate > 1) )
                 retval <- c( retval , "mutationRate is invalid")

               if( ( object@missingRate < 0 | object@missingRate > 1) )
                 retval <- c( retval , "missingRate is invalid")

               if( ( object@returnIBD != FALSE) & (object@returnIBD != TRUE) )
                 retval <- c( retval , "returnIBD invalid")

               if( object@intervalIBD <= 0 )
                 retval <- c( retval , "intervalIBD is invalid")

               if( is.null(retval) ) return(TRUE)
               else return(retval)
             }
)

#------------------------------------------------------------------------------
# 'markerFacts' accessor functions:
#------------------------------------------------------------------------------
setGeneric('markerFacts', function(object) standardGeneric('markerFacts'))
setMethod('markerFacts', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@markerFacts)
          }
)

setGeneric('markerFacts<-', function(object, value) standardGeneric('markerFacts<-'))
setMethod('markerFacts<-', signature(object = 'strumMarker'),
          function(object, value)
          {
            object@markerFacts <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'haplotypes' accessor functions:
#------------------------------------------------------------------------------
setGeneric('haplotypes', function(object) standardGeneric('haplotypes'))
setMethod('haplotypes', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@haplotypes)
          }
)

setGeneric('haplotypes<-', function(object, value) standardGeneric('haplotypes<-'))
setMethod('haplotypes<-', signature(object = 'strumMarker'),
          function(object, value)
          {
            object@haplotypes <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'populationRecombRate' accessor functions:
#------------------------------------------------------------------------------
setGeneric('populationRecombRate', function(object) standardGeneric('populationRecombRate'))
setMethod('populationRecombRate', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@populationRecombRate)
          }
)

setGeneric('populationRecombRate<-', function(object, value) standardGeneric('populationRecombRate<-'))
setMethod('populationRecombRate<-', signature(object = 'strumMarker'),
           function(object, value)
           {
             object@populationRecombRate <- value
             return(object)
           }
)

#------------------------------------------------------------------------------
# 'errorRate' accessor functions:
#------------------------------------------------------------------------------
setGeneric('errorRate', function(object) standardGeneric('errorRate'))
setMethod('errorRate', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@errorRate)
          }
)

setGeneric('errorRate<-', function(object, value) standardGeneric('errorRate<-'))
setMethod('errorRate<-', signature(object = 'strumMarker'),
          function(object, value)
          {
            object@errorRate <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'mutationRate' accessor functions:
#------------------------------------------------------------------------------
setGeneric('mutationRate', function(object) standardGeneric('mutationRate'))
setMethod('mutationRate', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@mutationRate)
          }
)

setGeneric('mutationRate<-', function(object, value) standardGeneric('mutationRate<-'))
setMethod('mutationRate<-', signature(object = 'strumMarker'),
          function(object, value)
          {
            object@mutationRate <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'missingRate' accessor functions:
#------------------------------------------------------------------------------
setGeneric('missingRate', function(object) standardGeneric('missingRate'))
setMethod('missingRate', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@missingRate)
          }
)

setGeneric('missingRate<-', function(object, value) standardGeneric('missingRate<-'))
setMethod('missingRate<-', signature(object = 'strumMarker'),
          function(object, value)
          {
            object@missingRate <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'coding' accessor functions:
#------------------------------------------------------------------------------
setGeneric('coding', function(object) standardGeneric('coding'))
setMethod('coding', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@coding)
          }
)

setGeneric('coding<-', function(object, value) standardGeneric('coding<-'))
setMethod('coding<-', signature(object = 'strumMarker'),
          function(object,value)
          {
            object@coding <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'returnIBD' accessor functions:
#------------------------------------------------------------------------------
setGeneric('returnIBD', function(object) standardGeneric('returnIBD'))
setMethod('returnIBD', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@returnIBD)
          }
)

setGeneric('returnIBD<-', function(object, value) standardGeneric('returnIBD<-'))
setMethod('returnIBD<-', signature(object = 'strumMarker'),
          function(object, value)
          {
            object@returnIBD <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'intervalIBD' accessor functions:
#------------------------------------------------------------------------------
setGeneric('intervalIBD', function(object) standardGeneric('intervalIBD'))
setMethod('intervalIBD', signature(object = 'strumMarker'),
          function(object)
          {
            return(object@intervalIBD)
          }
)

setGeneric('intervalIBD<-', function(object, value) standardGeneric('intervalIBD<-'))
setMethod('intervalIBD<-', signature(object = 'strumMarker'),
          function(object, value)
          {
            object@intervalIBD <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# show generic functions
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumMarker"),
          function(object) 
          {
            cat("\nBasic properties of the marker:\n")
            .printInfoLine("populationRecombRate", object@populationRecombRate, 40)
            .printInfoLine("errorRate", object@errorRate, 40)
            .printInfoLine("missingRate", object@missingRate, 40)
            .printInfoLine("mutationRate", object@mutationRate, 40)
            .printInfoLine("coding", paste(object@coding,collapse=' '), 40)
            .printInfoLine("returnIBD", object@returnIBD, 40)
            .printInfoLine("intervalIBD", object@intervalIBD, 40)

            cat("\nMarker facts contains information for ", nrow(object@markerFacts), "markers.\n")
            if( nrow(object@markerFacts) > 5 )
              cat("First 5 rows of marker facts:\n")
            else
              cat("Marker facts:\n")

            print(head(object@markerFacts, 5))

            cat("\nHaplotype contains information for ", ncol(object@haplotypes)," haplotypes.\n")
            if( ncol(object@haplotypes) > 10 )
            {
              if( nrow(object@markerFacts) > 5 )
                cat("First 10 haplotypes for 5 markers:\n")
              else
                cat("First 10 haplotypes:\n")
            } else
            {
              if( nrow(object@markerFacts) > 5 )
                cat("Haplotypes for 5 markers:\n")
              else
                cat("Haplotypes:\n")
            }

            print(object@haplotypes[1:5,1:10])
          }
)
