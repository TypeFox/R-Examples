gridCellEnhance = function(tract, 
                           gridLines = FALSE,      #can be too busy for high resolution
                           gridCenters = FALSE,    #can be too busy for high resolution
                           gridLineColor = .StemEnv$gridLineColor,
                           gridCenterColor = .StemEnv$gridCenterColor,
                           lwdGrid = 1,                                 #lwd for gridLines
                           ...
                          )
{
#---------------------------------------------------------------------------
#
#   This routine will plot gridLines and grid cell centers if desired. To
#   accomplish this, we just convert our "Tract" object to an object of
#   class "SpatialGridDataFrame" and use the built-in sp package routines.
#
#   Arguments...
#     tract = an object of class "Tract" or subclass--although a
#             'Rasterlayer' object will work too.
#     gridLines = TRUE: plot grid lines; FALSE: don't plot
#     gridCenters = TRUE: plot grid cell centers; FALSE: don't
#     gridLineColor = a legal color for grid lines
#     gridCenterColor = same for cell centers
#     lwGrid = the line width (i.e., lwd) for gridlines
#
#   Returns...
#     invisibly
#
#Author...									Date: 1-Oct-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   just a check...
#
    if(!gridLines && !gridCenters)
      return(invisible())
    if(!is(tract,'RasterLayer')) #since a "Tract" object is a "RasterLayer
      stop('Must have a "Tract" or "RasterLayer" object for grid lines!')
  
#
#   add some gridLines if desired, unless chainSawIZ or other tract with just
#   one grid cell...
#
    sg = as(tract, 'SpatialGridDataFrame')           #convert
    if(gridLines && ncell(tract)!= 1) {
      bb.sg = bbox(sg)
      rownames(bb.sg) = c('x','y')
      plot(gridlines(sg, easts= seq(bb.sg['x','min'], bb.sg['x','max'], len=sg@grid@cells.dim[1]+1) ,
           norths = seq(bb.sg['y','min'], bb.sg['y','max'], len=sg@grid@cells.dim[2]+1), 
           ndiscr=2 ), col=gridLineColor, add=TRUE, lwd=lwdGrid)
#     for some reason, not all bbox lines plot with gridLines, so augment...
      lines(rep(bb.sg[1,1],2),bb.sg[2,],col=gridLineColor, lwd=lwdGrid) #left
      lines(rep(bb.sg[1,2],2),bb.sg[2,],col=gridLineColor, lwd=lwdGrid) #right
      lines(bb.sg[1,],rep(bb.sg[2,1],2),col=gridLineColor, lwd=lwdGrid) #bottom
      lines(bb.sg[1,],rep(bb.sg[2,2],2),col=gridLineColor, lwd=lwdGrid) #top
    }
    if(gridCenters)
      plot(sg, col=gridCenterColor, add=TRUE)

    return(invisible())
}   #gridCellEnhance

