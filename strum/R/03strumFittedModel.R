#==============================================================================
# File: fittedModel.R
#
# Author: Nathan Morris
#
# Notes: fittedModel class definition & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of fittedModel class
#------------------------------------------------------------------------------
setClass("strumFittedModel",
         representation(
           myStrumModel              = "ANY",
           modelValidity             = "ANY",
           fittedParameters          = "ANY",
           fittedParametersCovMatrix = "ANY",
           deltaParameters           = "ANY",
           deltaParametersCovMatrix  = "ANY",
           parDiff                   = "ANY",
           parDiffCovMatrix          = "ANY",
           chiTestOut                = "ANY",
           fitIndices                = "ANY"),
         prototype = list()
) 

#------------------------------------------------------------------------------
# 'fittedParameters' accessor functions:
#------------------------------------------------------------------------------
setGeneric('fittedParameters', function(object) standardGeneric('fittedParameters'))
setMethod('fittedParameters', signature(object = 'strumFittedModel'),
          function(object)
          {
            return(object@fittedParameters)
          }
)

setGeneric('fittedParameters<-', function(object,value) standardGeneric('fittedParameters<-'))
setMethod('fittedParameters<-', signature(object = 'strumFittedModel'),
          function(object, value)
          {
            object@fittedParameters <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'fittedParametersCovMatrix' accessor functions:
#------------------------------------------------------------------------------
setGeneric('fittedParametersCovMatrix', function(object) standardGeneric('fittedParametersCovMatrix'))
setMethod('fittedParametersCovMatrix', signature(object = 'strumFittedModel'),
          function(object)
          {
            return(object@fittedParametersCovMatrix)
          }
)

setGeneric('fittedParametersCovMatrix<-', function(object,value) standardGeneric('fittedParametersCovMatrix<-'))
setMethod('fittedParametersCovMatrix<-', signature(object = 'strumFittedModel'),
          function(object, value)
          {
            object@fittedParametersCovMatrix<-value
            return(object)
          }
)

#------------------------------------------------------------------------------
# show generic functions
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumFittedModel"),
          function(object) 
          {
            cat("\n\n=========\n")
            cat("  Model  \n")
            cat("=========\n")

            .showModel(object@myStrumModel, "strumFittedModel")

            if( object@modelValidity == 1 )
            {
              cat("\n\n*** Not identifiable, more parameters in model than in saturated model! ***\n\n")

            } else if( object@modelValidity == 2 )
            {
              cat("\n\n*** Parameters are not locally identifiable! ***\n\n")

            } else
            {
              cat("\n\n==========\n")
              cat("  Result  \n")
              cat("==========\n")

              if( object@modelValidity == 3 )
              {
                cat("\nParameter estimates:\n")
                print(object@fittedParameters)

                cat("\n*** Same number of parameters in model and saturated model.  No fit test! ***\n\n")

              } else if( object@modelValidity == 4 )
              {
                cat("\n*** Negative variance detected! Possible numerical instability! ***\n\n")

              } else
              {
                cat("\nParameter estimates:\n")
                print(object@fittedParameters)

                cat("\nChi-square fit statistics:\n")
                print(object@chiTestOut)

                cat("\nModel fit indices:\n")
                print(object@fitIndices)
              }
            }
          }
)
