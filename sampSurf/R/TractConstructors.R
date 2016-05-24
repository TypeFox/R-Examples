#---------------------------------------------------------------------------
#
#   This file holds the S4 definition for the constructor methods of the
#   Tract class & subclasses...
#
#   The Tract methods include...
#
#     1. object=='missing': a constructor using cell dimensions
#     2. object=='numeric': for xy extents in vector form always assuming the
#                minimum/origin is at (0,0) (calls #1)
#     3. object=='matrix': for bbox matrix (calls #1) -- use this method (or #1)
#                for making tracts with something other than (0,0) minimum
#     4. object=='RasterLayer': one for a 'RasterLayer' class
#
#   bufferedTract methods...
#
#     1. bufferWidth and 'Tract' object
#
#Author...									Date: 31-Aug-2010
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
#if (!isGeneric("Tract")) 
  setGeneric('Tract',  
             function(object, ...) standardGeneric('Tract'),
             signature = c('object')
            )

  setGeneric('bufferedTract',  
             function(bufferWidth, tract, ...) standardGeneric('bufferedTract'),
             signature = c('bufferWidth', 'tract')
            )




          
#================================================================================
#  1. constructor for class Tract based on SpatialGridDataFrame constructor in sp...
#
setMethod('Tract',
          signature(object = 'missing'),
function(cellSize,
         cellDims = c(x=100, y=100), 
         cellCenter = NULL,
         units = 'metric',
         data = 0,
         spUnits = CRS(projargs=as.character(NA)),
         description = 'object of class Tract',
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   Arguments...
#     cellSize = the x,y cell size extents, uniform (square) in x,y ==> scalar
#     cellDims = number of cells in x,y
#     cellCenter = the origin cell's center x,y
#     runQuiet = TRUE: no summary printing; FALSE: print some info
#     ... = to pass on or ignore
#
#
    if(length(cellSize) > 1)
      stop('cellSize must be a numeric scalar!')
    cellSize = rep(cellSize, 2)  #uniform size in x,y == square cells only!
    names(cellSize) = c('x','y') #this sets up the names for x and y for the grid
    if(length(cellCenter) != 2)
      stop('cellCenter must be a numeric vector of length 2!')
    if(is.null(cellCenter) || is.na(cellCenter))
      cellCenter = cellSize/2
    if(length(cellDims) != 2)
      stop('cellDims argument must be a length 2 vector!')
    gt = GridTopology(cellCenter, cellSize, cellDims)
    
    
#
#   this will add a surface value to the GridTopology class to allow us
#   to convert easily to RasterLayer...
#
    dfLen = prod(gt@cells.dim)                     #x- and y-extends in # of cells
    if(is.na(data) || length(data) == 1)
      surf = rep(data, dfLen)
    else if(length(data) == dfLen)
      surf = data
    else
      stop('bad length of data argument!')
    df = data.frame( surf=surf )
    sgdf = SpatialGridDataFrame(gt, data=df, proj4string = spUnits)   #contains everything in GT plus the surface

    
    if(!runQuiet) {
      cat('\n\nGrid Topology...\n')
      print(gt)
      cat('\n\n')
    }

    ra = raster(sgdf)
    ra = as(ra, 'Tract')
    ra@description = description
    ra@units = units

    ra@area = nrow(ra)*ncol(ra)*xres(ra)^2

    if(validObject(ra))
      return(invisible(ra))
    else
      return(NULL)

}   #Tract constructor for
)   #setMethod






          
#================================================================================
#  2. constructor for class Tract with tract extents always starting at (0,0)...
#
setMethod('Tract',
          signature(object = 'numeric'),
function(object = c(x=10, y=10),  #grid maximum extents
         cellSize,
         units = 'metric',
         data = 0,
         spUnits = CRS(projargs=as.character(NA)),
         description = 'object of class Tract',
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   Arguments...
#     object = the actual x,y extents of the tract with lower bounds (0,0) always
#     cellSize = the x,y cell size extents, uniform (square) in x,y ==> scalar
#     runQuiet = TRUE: no summary printing; FALSE: print some info
#     ... = to pass on or ignore
#
#
    xyExtents = object
    if(!identical(intersect(names(xyExtents), c('x','y')), c('x','y')) || length(xyExtents)!=2)
      stop('xyExtents must be a length 2 vector with named "x" and "y" values!')
  
    cellDims = xyExtents/cellSize                 #overall cell dimensions          
    names(cellDims) = names(xyExtents)
    cellCenter = rep(cellSize/2, 2)               #minimum grid cell center point

    tract = Tract(cellSize=cellSize, cellDims=cellDims, cellCenter=cellCenter, units=units,
                  data=data, spUnits=spUnits, description=description, runQuiet=runQuiet,
                  ...
                 )
    return(invisible(tract))
}   #Tract constructor for maximal extents
)   #setMethod
    







          
#================================================================================
#  3. constructor for class Tract with tract bbox matrix...
#
setMethod('Tract',
          signature(object = 'matrix'),
function(object, 
         cellSize,
         units = 'metric',
         data = 0,
         spUnits = CRS(projargs=as.character(NA)),
         description = 'object of class Tract',
         runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#
#   Arguments...
#     object = the bbox extents of the tract which can have lower limits that
#              are not necessarily (0,0)
#     cellSize = the x,y cell size extents, uniform (square) in x,y ==> scalar
#     runQuiet = TRUE: no summary printing; FALSE: print some info
#     ... = to pass on or ignore
#
#
    bbox = object
    stopifnot(bboxCheck(bbox))


    cellDims = apply(bbox, 1, function(x) x['max']-x['min'])/cellSize  #overall cell dimensions
    cellCenter = bbox[,'min'] + rep(cellSize/2, 2)                     #minimum grid cell center point

    tract = Tract(cellSize=cellSize, cellDims=cellDims, cellCenter=cellCenter, units=units,
                  data=data, spUnits=spUnits, description=description, runQuiet=runQuiet,
                  ...
                 )
    return(invisible(tract))
}   #Tract constructor for bbox
)   #setMethod
    
    


  






          
#================================================================================
#  4. method for construction of class Tract using a RasterLayer object...
#
setMethod('Tract',
          signature(object = 'RasterLayer'),
function(object,
         units = 'metric',
         description = 'object of class Tract',
         ##runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------
#  
# could just set up a call to above, and invoke inheritance here,
# probably should eventually, but this is more efficient...
#
    ra = object
    if(!identical(xres(ra), yres(ra)))
      stop('Cells must be square, not rectangular!')     
   
    ra = as(ra, 'Tract')
    ra@description = description
    ra@units = units

    ra@area = nrow(ra)*ncol(ra)*xres(ra)^2

    if(validObject(ra))
      return(invisible(ra))
    else
      return(NULL)
}   #Tract constructor for 'RasterLayer'
)   #setMethod

    













          
#================================================================================
#================================================================================
#  method for functions and class bufferedTract...
#
setMethod('bufferedTract',
          signature(bufferWidth = 'numeric', tract = 'Tract'),
function(bufferWidth,
         tract,
         ##runQuiet = TRUE,
         ...
        )
{
#---------------------------------------------------------------------------

    if(bufferWidth <= 0)
      stop('bufferWidth must be positive')

#
#   make sure buffer width is rounded to cell resolution so the buffer will align
#   properly with a cell boundary/edge...
#
 #   cellRes = res(tract)[1]
 #   bufferWidth = round(bufferWidth/cellRes)*cellRes

    bb = bbox(tract)
    lims = apply(bb, 1, diff)
    intReg = lims - 2*bufferWidth
    if(any(lims < 2*bufferWidth))
      stop('bufferWidth is larger than tract extents!!')
       

    bt = as(tract, 'bufferedTract')
    
#   buffer rectangle as a matrix...
    bufferRect = matrix(c(bb[1,'min'] + bufferWidth, bb[1,'max'] - bufferWidth,
                          bb[2,'min'] + bufferWidth, bb[2,'max'] - bufferWidth),
                        nrow=2, byrow = TRUE,
                        dimnames = list(c('x','y'),c('min','max'))
                       )
    lims = apply(bufferRect, 1, diff)
    if(any(lims <= 0))
      stop('internal buffer region is zero or negative in either x or y!')
    
    bt@bufferRect = bufferRect
    
#   buffer region as a SpatialPolygon...
    #***use polygonFromExtent() instead ********** does it label x,y as such??
    sr = bufferRect
    #sr = rbind(sr[,1], sr[1,], sr[,2], rev(sr[2,]), sr[,1])  #closed polygon
    sr = rbind(sr[,'min'],
               c(sr['x','max'], sr['y','min']),
               c(sr['x','max'], sr['y','max']),
               c(sr['x','min'], sr['y','max']),
               sr[,'min'])  #closed polygon
    colnames(sr) = c('x','y')
    rownames(sr) = 1:5
    pgSR = Polygon(sr) 
    pgsSR = Polygons(list(sr=pgSR), ID='bufferRect')
    spSR = SpatialPolygons(list(pgsSR=pgsSR))      #takes a list of Polygons objects
    bt@spBuffer = spSR
    

    if(validObject(bt))
      return(invisible(bt))
    else
      return(NULL)

}   #bufferedTract constructor for
)   #setMethod

       
