#---------------------------------------------------------------------------
#
#   Returns the area of an object. Note that the "raster" package defines
#   the generic and defining one here also will result in a warning from R,
#   so I have just used their generic instead...
#
#   1. "InclusionZone" general, for most IZs
#   2. "standUpIZ" uses circularPlot
#   3. "chainSawIZ" is always a pain
#   4. "Tract" and subclasses
#   5. "InclusionZoneGrid"
#
#      ...standingTree related...
#   1. "circularPlotIZ"
#
#Author...									Date: 10-May-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#   generic definition...
#
#if(!isGeneric("area")) 
#  setGeneric('area', package='sampSurf',  
#             function(object, ...) standardGeneric('area'),
#             signature = c('object')
#            )
#^^^^raster has a generic area function already^^^

       
#================================================================================
#  method for general "InclusionZone" objects...
#
setMethod('area',
          signature(x = 'InclusionZone'),
function(x,
         ...
        )
{
#------------------------------------------------------------------------------
#   most of the newer inclusion zone objects just have an area slot...
#
    return(x@area)
}   #"InclusionZone" constructor
)   #setMethod




#================================================================================
#  method for "standUpIZ" objects...
#
setMethod('area',
          signature(x = 'standUpIZ'),
function(x,
         ...
        )
{
#------------------------------------------------------------------------------
#   gotta get it from the circularPlot object...
#
    
    return(x@circularPlot@area)
}   #"standUpIZ" constructor
)   #setMethod



#================================================================================
#  method for "chainsawIZ" objects...
#
setMethod('area',
          signature(x = 'chainSawIZ'),
function(x,
         ...
        )
{
#------------------------------------------------------------------------------
#   chainsaw is always wierd--the area is a point==0...
#
    
    return(0)
}   #"chainSawIZ" constructor
)   #setMethod



       
#================================================================================
#  method for "Tract" objects...
#
setMethod('area',
          signature(x = 'Tract'),
function(x,
         ...
        )
{
#------------------------------------------------------------------------------
#   just get the area slot...
#
    return(x@area)
}   #"Tract" constructor
)   #setMethod

       
#================================================================================
#  method for "InclusionZoneGrid" objects...
#
setMethod('area',
          signature(x = 'InclusionZoneGrid'),
function(x,
         ...
        )
{
#------------------------------------------------------------------------------
#   just determine the area from the extents...
#
    area = nrow(x@grid)*ncol(x@grid)*xres(x@grid)^2
    return(area)
}   #"InclusionZoneGrid" constructor
)   #setMethod




#================================================================================
#  method for "circularPlotIZ" objects...
#
setMethod('area',
          signature(x = 'circularPlotIZ'),
function(x,
         ...
        )
{
#------------------------------------------------------------------------------
#   gotta get it from the circularPlot object...
#
    
    return(x@circularPlot@area)
}   #"circularPlotIZ" constructor
)   #setMethod
