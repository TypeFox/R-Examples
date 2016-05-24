#---------------------------------------------------------------------------
#
#   This method should handle plotting of InclusionZoneGrid class objects
#   without modification if they have been constructed correctly. It simply
#   plots the background grid and then the InclusionZone object.
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
#





#================================================================================
#  method for objects of class InclusionZoneGrid and subclasses...
#
setMethod('plot',
          signature(x = 'InclusionZoneGrid', y='missing'),
function(x, 
         axes = TRUE,           #not a par() so can't be passed to callNextMethod, so separate it
         gridColor = .StemEnv$blue.colors(1000),
         asp = 1,
         izColor = NA,                                #default: no shading to obscure
         gridLines = TRUE,                            #plot gridlines?
         gridLineColor = .StemEnv$gridLineColor,
         gridCenters = FALSE,                         #can be too busy
         gridCenterColor = .StemEnv$gridCenterColor,
         estimate = names(c(.StemEnv$puaEstimates, .StemEnv$ppEstimates)),
         lwdGrid = 1,                                 #lwd for gridLines
         ...
        )
{
#------------------------------------------------------------------------------
#   first, establish the bounds by generating a raster object with nothing
#   for data to establish the plot extents; note that we can not just plot the
#   object bounding box using SpatialPoints, because it then messes up the
#   location of the legend in the final raster plot...
#
    ext = extent(bbox(x))
    cr = xres(x@grid)                              #square cells always    
    ext@xmin = ext@xmin - cr                       #pad each way
    ext@xmax = ext@xmax + cr
    ext@ymin = ext@ymin - cr
    ext@ymax = ext@ymax + cr
    nx = (ext@xmax - ext@xmin)/cr                   #remember x=columns
    ny = (ext@ymax - ext@ymin)/cr                   #and y=rows
    sr = raster(ext, nrows=as.integer(ny), ncols=as.integer(nx), crs=NA)
    sr[] = rep(0, ncell(sr))                        #NA gives warnings in plot, use 0 and col=NA below
    suppressWarnings(                               #for non-plot arguments in ...
    plot(sr, col=NA, legend=FALSE, axes=axes, asp=asp, ...) #set up limits first, legend comes later
                     )
 
#
#   now generate the raster plot for the given estimate variable desired...
#
    estimate = match.arg(estimate)
    grid = setValues(x@grid, x@data[,estimate]) 
    suppressWarnings(                                      #for non-plot arguments in ...
    plot(grid, col=gridColor, axes=axes, asp=asp, add=TRUE, ...) #call the raster method, with legend
                     )

#    
#   add some gridlines if desired, unless chainSawIZ for one grid cell...
#
    gridCellEnhance(grid, gridLines = gridLines, gridLineColor = gridLineColor,
                    gridCenters = gridCenters, gridCenterColor = gridCenterColor,
                    lwdGrid = lwdGrid, 
                    ...)
      

#   and the iz object...
    suppressWarnings(                                  #for non-plot arguments in ...
    plot(x@iz, add=TRUE, izColor=izColor, ...)
                     )
    
    
    return(invisible())

}    #plot for 'InclusionZoneGrid' base
) #setMethod



