#---------------------------------------------------------------------------
#
#   Methods for generic plot() for class...
#     (1) Stem and subclasses...
#     (2) downLogs subclass
#     (3) standingTree subclass
#
#   Note in the plotting of Spatial objects, the xlab and ylab arguments
#   do not seem to be passed through to plot, so no labels will be shown
#   when these arguments are specified. To get labels, use title(xlab=,
#   ylab=).
#
#Author...									Date: 9-Aug-2010
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
#  1.  method for class Stem...
#
setMethod('plot',
          signature(x = 'Stem', y='missing'),
function(x,
         pchStemLocation = 20,      #3 is also good
         stemLocationColor = .StemEnv$logAttributeColor,
         ...
        )
{
#------------------------------------------------------------------------------
#   just a simple plot of the x,y location from the virtual class...
#
#   since there will never be an object of class 'Stem', we can safely assume
#   that some components have already been plotted, and we can just
#   add to it with points()
#
#------------------------------------------------------------------------------
#
    object = x

    suppressWarnings(                                       #**note, for 'add' parameter passed in ...
      points(object@location, col=stemLocationColor,
             pch=pchStemLocation, ...)  #don't use plot() for this!
                    )
        
    return(invisible())

}    #plot for 'Stem'
) #setMethod






#================================================================================
#  2. method for downLog subclass...
#
setMethod('plot',
          signature(x = 'downLog', y='missing'),
function(x,
         axes = FALSE,  #not a par() so can't be passed to callNextMethod, so separate it
         logColor = .StemEnv$logColor,
         showLogCenter = FALSE,
         pchLogCenter = 3,       #20 is also good
         logCenterColor = .StemEnv$logAttributeColor,
         showNeedle = FALSE,
         logNeedleColor = .StemEnv$logAttributeColor,
         logBorderColor = .StemEnv$logBorderColor,   #log perimeter color
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the downed log...
#------------------------------------------------------------------------------
#
    object = x
    suppressWarnings({                              #for non-plot arguments in ...    
    plot(object@spLog, col=logColor, axes=axes, border=logBorderColor, asp=asp, ...)

    if(showNeedle)
      plot(object@slNeedleAxis, col=logNeedleColor,  #also a Spatial object
           add=TRUE)                 
    })   

# 
#   call next method subsequent, adding to the existing plot()...
#
    if(showLogCenter)
      callNextMethod(x=object, pchStemLocation = pchLogCenter,
                     stemLocationColor = logCenterColor, ...)  

        
    return(invisible())
}   #plot for 'downLog'
)   #setMethod


    



#================================================================================
#  3. method for standingTree subclass...
#
setMethod('plot',
          signature(x = 'standingTree', y='missing'),
function(x,
         axes = FALSE,  #not a par() so can't be passed to callNextMethod, so separate it
         treeColor = .StemEnv$treeColor,
         showLocation = TRUE,
         pchLocation = 3,       #20 is also good
         locationColor = .StemEnv$treeAttributeColor,
         treeBorderColor = .StemEnv$treeBorderColor,   #tree perimeter color
         asp = 1,
         ...
        )
{
#------------------------------------------------------------------------------
#  plots the trees...
#------------------------------------------------------------------------------
#
    object = x
    suppressWarnings({                              #for non-plot arguments in ...    
      plot(object@spDBH, col=treeColor, axes=axes, border=treeBorderColor, asp=asp, ...)
    })   


# 
#   call next method subsequent, adding to the existing plot()...
#
    if(showLocation)
      callNextMethod(x=object, pchStemLocation = pchLocation,
                     stemLocationColor = locationColor, ...)  

        
    return(invisible())
}   #plot for 'standingTree'
)   #setMethod






    


#*******************
#showMethods('plot')
#*******************
