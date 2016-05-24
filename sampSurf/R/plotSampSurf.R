#---------------------------------------------------------------------------
#
#   Methods for generic plot() for class for sampSurf class.
#
#   This uses the 'Tract" class plotting routine, so please see the info
#   there on extensions.
#
#Author...									Date: 5-Oct-2010
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
#  method for class 'sampSurf'...
#
setMethod('plot',
          signature(x = 'sampSurf', y='missing'),
function(x,
         showIZs = TRUE,
         izColor = NA,         #don't fill the interiors by default
         ...
        )
{
#------------------------------------------------------------------------------
#
#   plot both...
#
    suppressWarnings({                              #for non-plot arguments in ...    
      plot(x@tract, ...)
      if(showIZs)
        plot(x@izContainer, add=TRUE, izColor = izColor, ...)
    })   

    return(invisible())

}    #plot for 'sampSurf'
) #setMethod
    
