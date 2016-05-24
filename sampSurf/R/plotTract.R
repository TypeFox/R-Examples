#---------------------------------------------------------------------------
#
#   Methods for generic plot() for class for Tract and subclasses
#
#   Note in the plotting of Spatial objects, the xlab and ylab arguments
#   do not seem to be passed through to plot, so no labels will be shown
#   when these arguments are specified. To get labels, use title(xlab=,
#   ylab=).
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
#





#================================================================================
#  method for data frames and class Tract...
#
setMethod('plot',
          signature(x = 'Tract', y='missing'),
function(x, 
         axes = TRUE,           #not a par() so can't be passed to callNextMethod, so separate it
         gridColor = .StemEnv$blue.colors(1000),
         asp = 1,
         useImage = TRUE,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple plot of the tract portion of a Tract object (well, see comments)...
#
#   for now, set useImage=TRUE so that the graphic window will resize correctly
#   w/r to both the image and the vector graphics. Otherwise useImage=FALSE will
#   generate the gradient legend using the image routines from the fields package
#   and is fine for plotting a hard copy or if you do not resize the graphics device...
#

#
#   now a simple plot of the object, by adding it to the above dummy graphics window setup...
#
    suppressWarnings({                               #for non-plot arguments in ...
      if(useImage)                     
        image(x, asp=asp, col=gridColor, axes=axes, ...)
      else
        callNextMethod(x, col=gridColor, asp=asp, axes=axes, ...)   #RasterLayer plot method               
    })

#
#   this will plot grid lines and cell centers if desired, and if not too many cells...
#
    gridCellEnhance(x, ...)
                     
    
    return(invisible())

}    #plot for 'Tract'
) #setMethod





#================================================================================
#  method for bufferedTract...
#
setMethod('plot',
          signature(x = 'bufferedTract', y='missing'),
function(x, 
         bufferColor = transparentColorBase('blanchedalmond', .StemEnv$alphaTrans),
         axes = TRUE,           #not a par() so can't be passed to callNextMethod, so separate it
         gridColor = .StemEnv$blue.colors(1000),      
         lwd = 2,
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
  
    callNextMethod(x, gridColor=gridColor, axes=axes, ...)

    
    #don't plot the background within the buffer twice (col=NA)...
    suppressWarnings(                               #for non-plot arguments in ...    
      plot(x@spBuffer, add=TRUE, border=bufferColor, lwd = lwd, col=NA, ...)
    )
    
    return(invisible())

}    #plot for 'bufferedTract'
) #setMethod
  
        


