#==============================================================================
# File: strumModel.R
#
# Author: Nathan Morris
#
# Notes: strumModel class definition & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of strumModel class
#------------------------------------------------------------------------------
setClass("strumModel",
         contains = "strumVirtualModel"
)

#------------------------------------------------------------------------------
# show generic functions
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumModel"),
          function(object) 
          {
            .showModel(object, "strumModel")
          }
)

#------------------------------------------------------------------------------
# plot generic functions
#------------------------------------------------------------------------------
setMethod("plot", "strumModel",
          function(x, y, name="strumModel", toFile=TRUE, fileType="dot", ...) 
          {
            if (missing(y))
              y = "dot"

            .plotModel(x, layoutType=y, name=name, toFile=toFile, fileType=fileType, ...)
          }
)
