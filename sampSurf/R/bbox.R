#---------------------------------------------------------------------------
#
#   This extends the bbox methods found in both the sp and raster packages
#   to the sampSurf package class objects...
#
#   The methods for classes include...
#     1. InclusionZone & downLogIZs container
#     2. downLog
#     3. standingTree
#     4. StemContainer
#     5. circularPlot
#     6. Tract
#     7. sampSurf
#     8. lineSegment  (3-Oct-2012)
#
#Author...									Date: 17-Sept-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------

          
#================================================================================
#  method for 'InclusionZone' subclasses...
#
setMethod('bbox',
          signature(obj = 'InclusionZone'),
function(obj)
{
    return(obj@bbox)
}   #bbox for 'InclusionZone' only
)   #setMethod


#================================================================================
#  method for 'izContainer' subclasses...
#
setMethod('bbox',
          signature(obj = 'izContainer'),
function(obj)
{
    return(obj@bbox)
}   #bbox for 'izContainer'
)   #setMethod

          
#================================================================================
#  method for 'InclusionZoneGrid' classes...
#
setMethod('bbox',
          signature(obj = 'InclusionZoneGrid'),
function(obj)
{
    return(obj@bbox)
}   #bbox for 'InclusionZoneGrid' only
)   #setMethod



#================================================================================
#  method for 'downLog' class...
#
setMethod('bbox',
          signature(obj = 'downLog'),
function(obj)
{
    return(bbox(obj@spLog))
}   #bbox for 'downLog' only
)   #setMethod


#================================================================================
#  method for 'standingTree' class...
#
setMethod('bbox',
          signature(obj = 'standingTree'),
function(obj)
{
    return(bbox(obj@spDBH))
}   #bbox for 'standingTree' only
)   #setMethod


#================================================================================
#  method for 'StemContainer' class and subclasses...
#
setMethod('bbox',
          signature(obj = 'StemContainer'),
function(obj)
{
    return(obj@bbox)
}   #bbox for 'StemContainer' & subclasses
)   #setMethod
   

#================================================================================
#  method for 'circularPlot' class...
#
setMethod('bbox',
          signature(obj = 'circularPlot'),
function(obj)
{
    return(bbox(obj@perimeter))
}   #bbox for 'circularPlot' only
)   #setMethod
   

#================================================================================
#  method for 'Tract' class and subclasses...
#
setMethod('bbox',
          signature(obj = 'Tract'),
function(obj)
{
    bb = callNextMethod(obj)
    rownames(bb) = c('x','y')

    return(bb)
}   #bbox for 'Tract' and subclasses
)   #setMethod
    

#================================================================================
#  method for 'sampSurf' class and subclasses...
#
setMethod('bbox',
          signature(obj = 'sampSurf'),
function(obj)
{
    bb = bbox(obj@tract)   #should call 'Tract' method

    return(bb)
}   #bbox for 'Tract' and subclasses
)   #setMethod
   


#================================================================================
#  method for 'lineSegment' class...
#
setMethod('bbox',
          signature(obj = 'lineSegment'),
function(obj)
{
    return(bbox(obj@segment))
}   #bbox for 'lineSegment' only
)   #setMethod

