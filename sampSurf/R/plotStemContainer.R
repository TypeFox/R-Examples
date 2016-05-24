#---------------------------------------------------------------------------
#
#   This file holds the definitions for the basic plotting of objects of
#   class "StemContainer" and subclasses...
#
#   1. StemContainer -- virtual, sets up the plot extents if necessary
#   2. downLogs subclass
#   3. standingTrees subclass
#
#Author...									Date: 26-Oct-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#





#================================================================================
#  1.  method for class StemContainer...
#
setMethod('plot',
          signature(x = 'StemContainer', y='missing'),
function(x,
         axes = FALSE,         #not a par() so can't be passed to callNextMethod, so separate it
         add = FALSE,          #no existing plot assumed
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#
#   just a simple plot of the invisible extents from the virtual class...
#
#------------------------------------------------------------------------------
#
#   set up the plot extents via the bounding box...
#
    object = x
    if(!add) {
      bbox = object@bbox
      plot(SpatialPoints(t(bbox)), col=NA, axes=axes, asp=asp)  #set up limits
    }
    
    return(invisible())
}   #plot for 'StemContainer'
)   #setMethod












#================================================================================
#================================================================================
#  2. method for a "downLogs" collection/population...
#
setMethod('plot',
          signature(x = 'downLogs', y='missing'),
function(x, 
         axes = FALSE,        #not a par() so can't be passed to callNextMethod, so separate it
         add = FALSE,          #no existing plot assumed
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots all the downed logs in the collection...
#------------------------------------------------------------------------------
#
    object = x
    numLogs = length(object@logs)

    callNextMethod(x=x, axes=axes, add=add, asp=asp, ...)

    suppressWarnings({                              #for non-plot arguments in ...    
    for(i in seq_len(numLogs)) 
      plot(object@logs[[i]], axes=axes, add = TRUE, ...)
    })   
        
    return(invisible())
}   #plot for 'downLogs'
)   #setMethod








#================================================================================
#================================================================================
#  3. method for a "standingTrees" collection/population...
#
setMethod('plot',
          signature(x = 'standingTrees', y='missing'),
function(x, 
         axes = FALSE,         #not a par() so can't be passed to callNextMethod, so separate it
         add = FALSE,          #no existing plot assumed
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots all the standing trees in the collection...
#------------------------------------------------------------------------------
#
    object = x
    numTrees = length(object@trees)

    callNextMethod(x=x, axes=axes, add=add, asp=asp, ...)

    suppressWarnings({                              #for non-plot arguments in ...    
      for(i in seq_len(numTrees)) 
        plot(object@trees[[i]], axes=axes, add = TRUE, ...)
    })   
        
    return(invisible())
}   #plot for 'standingTrees'
)   #setMethod
    
