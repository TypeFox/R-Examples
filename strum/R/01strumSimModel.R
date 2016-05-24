#==============================================================================
# File: strumSimModel.R
#
# Author: Nathan Morris
#
# Notes: strumSimModel class definition & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of strumSimModel class
#------------------------------------------------------------------------------
setClass("strumSimModel",
         representation(
           markerInfo = "ANY",
           traitMissingRate = "numeric"),
         contains = "strumVirtualModel",
         prototype = list(
           markerInfo = NULL,
           traitMissingRate = numeric())
)

#------------------------------------------------------------------------------
# 'markerInfo' accessor functions:
#------------------------------------------------------------------------------
setGeneric('markerInfo', function(object) standardGeneric('markerInfo'))
setMethod('markerInfo', signature(object = 'strumSimModel'),
          function(object)
          {
            return(object@markerInfo)
          }
)

setGeneric('markerInfo<-', function(object,value) standardGeneric('markerInfo<-'))
setMethod('markerInfo<-', signature(object = 'strumSimModel'),
          function(object, value)
          {
            object@markerInfo <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'missingRate' accessor functions:
#------------------------------------------------------------------------------
setGeneric('traitMissingRate', function(object) standardGeneric('traitMissingRate'))
setMethod('traitMissingRate', signature(object = 'strumSimModel'),
          function(object)
          {
            return(object@traitMissingRate)
          }
)

setGeneric('traitMissingRate<-', function(object,value) standardGeneric('traitMissingRate<-'))
setMethod('traitMissingRate<-', signature(object = 'strumSimModel'),
          function(object, value)
          {
            object@traitMissingRate <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# show generic functions
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumSimModel"),
          function(object) 
          {
            .showModel(object, "strumSimModel")

            cat("\nTrait missing rate:\n")
            varY = object@varList$name[object@varList$inY]

            for( i in 1:length(object@traitMissingRate) )
            {
              yName = varY[i]
              mRate = object@traitMissingRate[i]

              .printInfoLine(yName, mRate, 25, 1)
            }

            cat("\nMarker info:\n")
            show(object@markerInfo)
          }
)

#------------------------------------------------------------------------------
# plot generic functions
#------------------------------------------------------------------------------
setMethod("plot", "strumSimModel",
          function(x, y, name="strumSimModel", toFile=TRUE, fileType="dot", ...) 
          {
            if (missing(y))
              y = "dot"

            .plotModel(x, layoutType=y, name=name, toFile=toFile, fileType=fileType, ...)
          }
)
