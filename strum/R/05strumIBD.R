#==============================================================================
# File: strumIBD.R
#
# Author: Nathan Morris
#
# Notes: strumIBD class definition & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of strumIBD class
#------------------------------------------------------------------------------
setClass("strumIBD",
         representation(
           ibdMarker = "data.frame",
           ibdMatrix = "list"),
         prototype = list(
           ibdMarker = data.frame(),
           ibdMatrix = list())
)

#------------------------------------------------------------------------------
# Check validity of strumIBD object
#------------------------------------------------------------------------------
setValidity ("strumIBD",
             function(object)
             {
               retval <- NULL
               if( is.null(object@ibdMarker) )
               {
                 retval <- c( retval , "ibdMarker is invalid")
               }
               if( is.null(object@ibdMatrix) )
               {
                 retval <- c( retval , "ibdMatrix list invalid")
               }

               if( is.null(retval) ) return(TRUE)
               else return(retval)
             }
)

#------------------------------------------------------------------------------
# 'ibdMarker' accessor functions:
#------------------------------------------------------------------------------
setGeneric('ibdMarker', function(object) standardGeneric('ibdMarker'))
setMethod('ibdMarker', signature(object = 'strumIBD'),
          function(object)
          {
            return(object@ibdMarker)
          }
)

setGeneric('ibdMarker<-', function(object,value) standardGeneric('ibdMarker<-'))
setMethod('ibdMarker<-', signature(object = 'strumIBD'),
          function(object, value)
          {
            object@ibdMarker <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'ibdMatrix' accessor functions:
#------------------------------------------------------------------------------
setGeneric('ibdMatrix', function(object) standardGeneric('ibdMatrix'))
setMethod('ibdMatrix', signature(object = 'strumIBD'),
          function(object)
          {
            return(object@ibdMatrix)
          }
)

setGeneric('ibdMatrix<-', function(object,value) standardGeneric('ibdMatrix<-'))
setMethod('ibdMatrix<-', signature(object = 'strumIBD'),
          function(object, value)
          {
            object@ibdMatrix <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# show generic functions
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumIBD"),
          function(object) 
          {
            if( length(object@ibdMarker) == 0 )
              cat("Empty IBD object.\n")
            else
            {
              cat("IBD object contains ", nrow(object@ibdMarker)," markers:\n")
              if( nrow(object@ibdMarker) > 5 )
                cat("First 5 rows of markers:\n")
              else
                cat("Markers:\n")

              print(head(object@ibdMarker,5))

              cat("First matrix: \n")
              print(head(object@ibdMatrix[[1]],1))
            }
          }
)
