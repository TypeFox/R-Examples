#---------------------------------------------------------------------------
#
#   Methods for generic plot() for class for izContainer subclasses...
#
#   All y="missing" signature, x is...
#
#   1. "izContainer" or any subclass
#
#
#Author...									Date: 23-Aug-2010
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
#================================================================================
#  method for a "izContainer" collection/population...
#
setMethod('plot',
          signature(x = 'izContainer', y='missing'),
function(x, 
         axes = FALSE,        #not a par() so can't be passed to callNextMethod, so separate it
         add = FALSE,          #add each log to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the downed logs...
#------------------------------------------------------------------------------
#
    object = x
    numIZs = length(object@iZones)

#
#   set up the plot extents via the bounding box, then plot each inclusion zone...
#
    suppressWarnings({                                #for object-specific parameters not in par() ...
    bbox = object@bbox
    if(!add)
      plot(SpatialPoints(t(bbox)), col=NA, axes=axes, asp=asp, ...)  #set up limits

    for(i in seq_len(numIZs))  #use one of the individual iz methods for plotting...
      plot(object@iZones[[i]], axes=axes, add=TRUE, ...)
    
    })
        
    return(invisible())
}   #plot for 'izContainer'
)   #setMethod
