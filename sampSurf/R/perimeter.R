#---------------------------------------------------------------------------
#
#   These methods are used to get the perimeter spatial object from a given
#   InclusionZone object. It could also be used to get a perimeter from a
#   downLog object, etc. but there's really not much need for the latter
#   as that is simple to retrieve...
#
#   The methods include (in no particular order)...
#     1. standUpIZ
#     2. chainSawIZ: either the circular plot location IZ or the sausage IZ
#                    are possible returns
#     3. sausageIZ
#     4. pointRelascopeIZ
#     5. perpendicularDistanceIZ (omnibusPDSIZ, distanceLimitedPDSIZ, omnibusDLPDSIZ)
#     6. distanceLimitedIZ
#     7. downLog
#     8. standingTree
#     9. StemContainer subclasses
#    10. circularPlot of class ArealSampling
#    11. Tract
#    12. sampSurf
#    13. izContainer
#    14. circularPlotIZ
#
#Author...									Date: 21-Sept-2010
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
if(!isGeneric('perimeter')) {
  setGeneric('perimeter',  
             function(object, ...) standardGeneric('perimeter'),
             signature = c('object')
            )
}





#================================================================================
#  method for standUpIZ object...
#
setMethod('perimeter',
          signature(object = 'standUpIZ'),
function(object, ...)
{
    return(object@circularPlot@perimeter)
}   #standUpIZ
) #setMethod






#================================================================================
#  method for chainSawIZ object...
#
setMethod('perimeter',
          signature(object = 'chainSawIZ'),
function(object,
         whatSense = c('point', 'plot', 'sausage'),
         ...)
{
    whatSense = match.arg(whatSense)
    if(whatSense == 'point')                    #return the plot centerpoint
      return(object@circularPlot@location)
    else if(whatSense == 'plot')                #return the plot
      return(perimeter(object@circularPlot))
    else {                                      #full sausage inclusion zone
      sa = sausageIZ(object@downLog, object@circularPlot@radius,
                     spUnits = object@circularPlot@spUnits, ...)
      return(sa@perimeter)
    }
      
}   #chainSawIZ
) #setMethod






#================================================================================
#  method for sausageIZ object...
#
setMethod('perimeter',
          signature(object = 'sausageIZ'),
function(object, ...)
{
    return(object@perimeter)
}   #sausage
) #setMethod



#================================================================================
#  method for pointRelascopeIZ object...
#
setMethod('perimeter',
          signature(object = 'pointRelascopeIZ'),
function(object, ...)
{
    return(object@perimeter)
}   #pointRelascope
) #setMethod



#================================================================================
#  method for perpendicularDistanceIZ object or subclasses  (omnibusPDSIZ,
#  distanceLimitedPDSIZ, omnibusDLPDSIZ)...
#
setMethod('perimeter',
          signature(object = 'perpendicularDistanceIZ'),
function(object, ...)
{
    return(object@perimeter)
}   #perpendicularDistance
) #setMethod



#================================================================================
#  method for distanceLimitedIZ or subclass (distanceLimitedMCIZ) object...
#
setMethod('perimeter',
          signature(object = 'distanceLimitedIZ'),
function(object, ...)
{
    return(object@perimeter)
}   #distanceLimitedIZ
) #setMethod



#================================================================================
#  method for downLog object...
#
setMethod('perimeter',
          signature(object = 'downLog'),
function(object, ...)
{
    return(object@spLog)
}   #downLog
) #setMethod





#================================================================================
#  method for standingTree object...
#
setMethod('perimeter',
          signature(object = 'standingTree'),
function(object, ...)
{
    return(object@spDBH)
}   #standingTree
) #setMethod



#================================================================================
#  method for StemContainer subclass objects...
#
setMethod('perimeter',
          signature(object = 'StemContainer'),
function(object, ...)
{
    return(bboxToPoly(object))
}   #StemContainer subclasses
) #setMethod




#================================================================================
#  method for circularPlot object...
#
setMethod('perimeter',
          signature(object = 'circularPlot'),
function(object, ...)
{
    return(object@perimeter)
}   #circularPlot
) #setMethod



#================================================================================
#  method for Tract object...
#
setMethod('perimeter',
          signature(object = 'Tract'),
function(object, ...)
{
    return(bboxToPoly(object))
}   #Tract
) #setMethod



#================================================================================
#  method for sampSurf object...
#
setMethod('perimeter',
          signature(object = 'sampSurf'),
function(object, ...)
{
    perimeter = perimeter(object@tract)  #a 'Tract' object
  
    return(perimeter)
}   #sampSurf
) #setMethod



#================================================================================
#  method for izContainer object...
#
setMethod('perimeter',
          signature(object = 'izContainer'),
function(object, ...)
{
    return(bboxToPoly(object))
}   #izContainer
) #setMethod





#================================================================================
#  method for circularPlotIZ or horizontalPointIZ object...
#
setMethod('perimeter',
          signature(object = 'circularPlotIZ'),
function(object, ...)
{
    return(object@circularPlot@perimeter)
}   #circularPlotIZ
) #setMethod
