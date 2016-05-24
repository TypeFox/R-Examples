#---------------------------------------------------------------------------
#
#   Methods for generic plot() for class for ArealSampling and subclasses;
#   this includes the container classes
#
#   1. ArealSampling base class
#   2. circularPlot class
#   3. lineSegment class    (3-Oct-2012)
#
#   Note in the plotting of Spatial objects, the xlab and ylab arguments
#   do not seem to be passed through to plot, so no labels will be shown
#   when these arguments are specified. To get labels, use title(xlab=,
#   ylab=).
#
#Author...									Date: 19-Aug-2010
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
#  method for data frames and class ArealSampling...
#
setMethod('plot',
          signature(x = 'ArealSampling', y='missing'),
function(x, 
         pchIZCenter = 20,      #3 is also good
         izCenterColor = .StemEnv$izCenterColor,
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple plot of the x,y location from the virtual class...
#
#   since there will never be an object of class 'ArealSampling', we can safely assume
#   that some components have been plotted already, and we can just
#   add to it with points()
#
#------------------------------------------------------------------------------
#
    object = x

    suppressWarnings(                         #**note, for 'add' parameter passed in ...
      points(object@location, col=izCenterColor, pch=pchIZCenter,
             asp=asp,  ...)                   #don't use plot() for this!
                    )
        
    return(invisible())

}    #plot for 'ArealSampling'
) #setMethod




#================================================================================
#  method for circularPlot subclass...
#
setMethod('plot',
          signature(x = 'circularPlot', y='missing'),
function(x, 
         axes = FALSE,           #not a par() so can't be passed to callNextMethod, so separate it
         izColor = .StemEnv$izColor,
         pchPlotCenter = 3,          #20 is also good
         showPlotCenter = FALSE,
         showPerimeter = TRUE,
         borderColor = .StemEnv$izBorderColor,     #plot perimeter color
         plotCenterColor = .StemEnv$izCenterColor,
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the circularPlot object...
#------------------------------------------------------------------------------
#
    object = x

    suppressWarnings(                                #for object-specific parameters not in par() ...
    if(showPerimeter)
      plot(object@perimeter, col=izColor, axes=axes, border=borderColor, asp=asp, ...)
    )                 

    if(showPlotCenter)
      callNextMethod(x=object, pchIZCenter=pchPlotCenter, asp=asp, ...)
        
    return(invisible())
}   #plot for 'circularPlot'
)   #setMethod






#================================================================================
#  method for lineSegment subclass...
#
setMethod('plot',
          signature(x = 'lineSegment', y='missing'),
function(x, 
         axes = FALSE,           #not a par() so can't be passed to callNextMethod, so separate it
         pchLineCenter = 20,
         showLineCenter = FALSE,
         showLineSegment = TRUE,
         lineColor = .StemEnv$izBorderColor, 
         lineCenterColor = .StemEnv$izCenterColor,
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the lineSegment object...
#------------------------------------------------------------------------------
#
    object = x

    suppressWarnings(                                #for object-specific parameters not in par() ...
    if(showLineSegment)
      plot(object@segment, col=lineColor, axes=axes, asp=asp, ...)
    )                 

    if(showLineCenter)
      callNextMethod(x=object, pchIZCenter=pchLineCenter, asp=asp, ...)
        
    return(invisible())
}   #plot for 'lineSegment'
)   #setMethod
