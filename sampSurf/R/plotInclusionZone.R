#---------------------------------------------------------------------------
#
#   Methods for generic plot() for class for InclusionZone and subclasses;
#
#   All y="missing" signature, x is...
#
#   1. "InclusionZone"
#
#       ...downLogIZ subclasses...
#   2. "StandUpIZ"
#   3. "chainSawIZ"
#   4. "sausageIZ"
#   5. "pointRelascopeIZ"
#   6. "perpendicularDistanceIZ"
#   7. "distanceLimitedPDSIZ'
#   8. "distanceLimitedIZ'
#
#   Note that the plot routines for methods like standup, sausage, PRS, etc.
#   are all very similar, they could probably be collapsed sometime with a
#   "base" function taking away much of the duplication. 18-Jan-2011
#
#       ...standingTreeIZ subclasses..
#   1.  "circularPlotIZ"
#   2.  "horizontalPointIZ" is a subclass of "circularPlotIZ" so has no
#       subclass plot at this point, it just uses the superclass version
#
#   Note in the plotting of Spatial objects, the xlab and ylab arguments
#   do not seem to be passed through to plot, so no labels will be shown
#   when these arguments are specified. To get labels, use title(xlab=,
#   ylab=).
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
#  1. method for data frames and class InclusionZone...
#
setMethod('plot',
          signature(x = 'InclusionZone', y='missing'),
function(x, 
         axes = FALSE,        
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#   this just sets up the limits of the plot w/r to the bbox if required
#------------------------------------------------------------------------------
#
    object = x

    bbox = object@bbox
    suppressWarnings(                                     #for object-specific parameters not in par() ...
    plot(SpatialPoints(t(bbox)), col=NA, axes=axes, asp=asp, ...)  #set up limits
                    )

    return(invisible())

}    #plot for 'InclusionZone'
) #setMethod
     






#================================================================================
#  2. method for standUpIZ subclass...
#
setMethod('plot',
          signature(x = 'standUpIZ', y='missing'),
function(x, 
         axes = FALSE,           #not a par() so can't be passed to callNextMethod, so separate it
         showLog = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the standUpIZ object...
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@circularPlot, izColor=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showLog)
        plot(object@downLog, add=TRUE, ...)
    })
                     
    return(invisible())
}   #plot for 'standUpIZ'
)   #setMethod







#================================================================================
#  3. method for chainSawIZ subclass...
#
setMethod('plot',
          signature(x = 'chainSawIZ', y='missing'),
function(x, 
         axes = FALSE,           #not a par() so can't be passed to callNextMethod, so separate it
         showLog = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         showSliver = TRUE,
         ltySliver = 'dashed',
         sliverBorder = 'black',
         sliverColor = transparentColorBase('coral', .StemEnv$alphaTrans),
         showBolt = TRUE,
         ltyBolt = 'dotted',
         boltBorder =  transparentColorBase('grey20', .StemEnv$alphaTrans),
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the chainSawIZ object...
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@circularPlot, izColor=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showLog) 
        plot(object@downLog, add=TRUE, ...)
      if(showSliver)
        plot(object@sliver, lty=ltySliver, border=sliverBorder, col=sliverColor, add=TRUE)
      if(showBolt)
        plot(object@bolt$spBolt, lty=ltyBolt, border=boltBorder, add=TRUE)
      
    })
                     
    return(invisible())
}   #plot for 'chainSawIZ'
)   #setMethod





#================================================================================
#  4. method for sausageIZ subclass...
#
setMethod('plot',
          signature(x = 'sausageIZ', y='missing'),
function(x, 
         axes = FALSE,                     #not a par() so can't be passed to callNextMethod, so separate it
         showLog = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#   plots the sausageIZ object...
#
#   note: the plot center and log center coincide with sausage, if you want to
#         show one of them, show the log center
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@perimeter, col=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showLog)
        plot(object@downLog, add=TRUE, ...)
    })
                     
    return(invisible())
}   #plot for 'sausageIZ'
)   #setMethod






#================================================================================
#  5. method for pointRelascopeIZ subclass (18-Jan-2011)...
#
setMethod('plot',
          signature(x = 'pointRelascopeIZ', y='missing'),
function(x, 
         axes = FALSE,                     #not a par() so can't be passed to callNextMethod, so separate it
         showLog = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         showDualCenters = FALSE,               #show centers of dual circles
         dcColor = .StemEnv$izBorderColor,      #color for dual circle centers
         ...
        )
{
#------------------------------------------------------------------------------
#   plots the pointRelascopeIZ object...
#
#   note: the IZ center and log center coincide with PRS, if you want to
#         show one of them, show the log center
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@perimeter, col=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showLog)
        plot(object@downLog, add=TRUE, ...)
    })

    if(showDualCenters)
      plot(SpatialPoints(object@dualCenters), add=TRUE, col=dcColor, ...)
                     
    return(invisible())
}   #plot for 'pointRelascopeIZ'
)   #setMethod





#================================================================================
#  6. method for perpendicularDistnaceIZ subclass (18-Jan-2011)...
#
setMethod('plot',
          signature(x = 'perpendicularDistanceIZ', y='missing'),
function(x, 
         axes = FALSE,                     #not a par() so can't be passed to callNextMethod, so separate it
         showLog = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#   plots the perpendicularDistanceIZ object...
#
#   note: the IZ center and log center coincide with PDS, if you want to
#         show one of them, show the log center
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@perimeter, col=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showLog)
        plot(object@downLog, add=TRUE, ...)
    })
                     
    return(invisible())
}   #plot for 'perpendicularDistanceIZ'
)   #setMethod




#================================================================================
#  7. method for distanceLimitedPDSIZ or subclasses (10-Mar-2011)...
#
setMethod('plot',
          signature(x = 'distanceLimitedPDSIZ', y='missing'),
function(x, 
         axes = FALSE,                     #not a par() so can't be passed to callNextMethod, so separate it
         showLog = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         showFullPDSIZ = FALSE,
         showDLSPart = FALSE,
         showPDSPart = FALSE,
         ...
        )
{
#------------------------------------------------------------------------------
#   plots the distanceLimitedPDSIZ object...
#
#   note: the IZ center and log center coincide with PDS, if you want to
#         show one of them, show the log center
#------------------------------------------------------------------------------
#
    object = x

#
#   set up extents and plot the everything PDS...
#
    if(!add) 
      if(showFullPDSIZ && !is.null(object@pdsFull))    #should never be NULL, but just in case
        object@bbox = bbox(object@pdsFull)             #need more room for full inclusion zone to fit
    
    callNextMethod(object, axes=axes, asp=asp, showLog=showLog,  #will show the log if desired
                   izColor=izColor, izBorder=izBorder,
                   add=add, ...)                                #if add=T: don't re-adjust extents

#
#   the following will just plot the perimeters for clarity...
#
    suppressWarnings({                                 #for object-specific parameters not in par() ...
      if(showPDSPart && !is.null(object@pdsPart))
        plot(perimeter(object@pdsPart), col=NULL, border=izBorder, axes=axes, add=TRUE, ...)
      if(showDLSPart && !is.null(object@dlsPart))
        plot(perimeter(object@dlsPart), col=NULL, border=izBorder, axes=axes, add=TRUE, ...)
      if(showFullPDSIZ && !is.null(object@pdsFull))
        plot(perimeter(object@pdsFull), col=NULL, border=izBorder, axes=axes, add=TRUE, ...)
     # if(showLog && !add)                              #only if not already shown above
      #  plot(object@downLog, add=TRUE, ...)
    })
                     
    return(invisible())
}   #plot for 'distanceLimitedPDSIZ'
)   #setMethod





#================================================================================
#  8. method for distanceLimitedIZ or distanceLimitedMCIZ subclass (22-Mar-2011)...
#
setMethod('plot',
          signature(x = 'distanceLimitedIZ', y='missing'),
function(x, 
         axes = FALSE,                     #not a par() so can't be passed to callNextMethod, so separate it
         showLog = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#   plots the distanceLimitedIZ object...
#
#   note: the IZ center and log center coincide with VPMC, if you want to
#         show one of them, show the log center
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@perimeter, col=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showLog)
        plot(object@downLog, add=TRUE, ...)
    })
                     
    return(invisible())
}   #plot for 'distanceLimitedIZ'
)   #setMethod


















#---------------------------------------------------------------------------
#
#   plot methods for "standingTreeIZ" subclasses...
#
#   1. circularPlotIZ
#      Note as above that "horizontalPointIZ" is a subclass of "circularPlotIZ"
#      and therefore simply uses its plot method.
#
#Author...									Date: 1-Dec-2011
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
#  1. method for circularPlotIZ subclass...
#
setMethod('plot',
          signature(x = 'circularPlotIZ', y='missing'),
function(x, 
         axes = FALSE,           #not a par() so can't be passed to callNextMethod, so separate it
         showTree = TRUE,
         izColor = .StemEnv$izColor,
         izBorder = .StemEnv$izBorderColor,
         add = FALSE,                           #add each IZ to overall plot if TRUE
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the circularPlotIZ object...
#------------------------------------------------------------------------------
#
    object = x

    if(!add)
      callNextMethod(object, axes=axes, asp=asp, ...)            #setup extents

    suppressWarnings({                                #for object-specific parameters not in par() ...
      plot(object@circularPlot, izColor=izColor, border=izBorder, axes=axes, add=TRUE, ...)
      if(showTree)
        plot(object@standingTree, add=TRUE, ...)
    })
                     
    return(invisible())
}   #plot for 'circularPlotIZ'
) #setMethod
