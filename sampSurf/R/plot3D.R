#---------------------------------------------------------------------------
#
#   Methods for generic plot3D() rgl driver for class for "sampSurf" &
#   "Tract" classes.
#
#   "InclusionZoneGrid" class added 18-Apr-2011
#
#   Note: the only reason to have a method below for "Tract" is so that
#         the blue.colors() are used by default rather than terrain.colors(),
#         which it would default to using the "RasterLayer" method, if an
#         object of class "Tract" of subclass is passed directly to plot3D.
#
#Author...									Date: 18-Oct-2010
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
setMethod('plot3D',
          signature(x = 'sampSurf'),
function(x,
         col = .StemEnv$blue.colors, #pass the actual function here, not a call!
         ...
        )
{
#------------------------------------------------------------------------------
#
    suppressWarnings({                              #for non-plot arguments in ...    
      plot3D(x@tract, col=col, ...)
    })   

    return(invisible())

}    #plot3D for 'sampSurf'
) #setMethod



#================================================================================
#  method for class 'Tract'...
#
setMethod('plot3D',
          signature(x = 'Tract'),
function(x,
         col = .StemEnv$blue.colors, #pass the actual function here, not a call!
         ...
        )
{
#------------------------------------------------------------------------------
#
    suppressWarnings({                              #for non-plot arguments in ...    
      callNextMethod(x, col=col, ...)
    })   

    return(invisible())

}    #plot3D for 'Tract'
) #setMethod



#================================================================================
#  method for class 'InclusionZoneGrid'...
#
setMethod('plot3D',
          signature(x = 'InclusionZoneGrid'),
function(x,
         estimate = 'volume', 
         col = .StemEnv$blue.colors, #pass the actual function here, not a call!
         ...
        )
{
#------------------------------------------------------------------------------
#
#   note that for the legal names below, we are not using, e.g.,
#   estimate = names(c(.StemEnv$puaEstimates, .StemEnv$ppEstimates)),
#   but are allowing whatever the creator of the particular sampling object put in
#   the data frame; thus, we are not using match.arg() here either...
#
    legalNames =  colnames(x@data)          
    edx = pmatch(estimate,legalNames)
    if(is.na(edx))
      stop('estimate must be one of:',paste(legalNames,collapse=','))
    else
      estimate = legalNames[edx]
    #estimate = match.arg(estimate)

#
#   check to see if this attribute has any estimates...
#    
    values = x@data[,estimate]
    values = values[values>0]
    if(all(is.na(values)))
      stop('All values for ',estimate,' are NA')

#
#   go ahead and plot it...
#
    x@grid = setValues(x@grid, x@data[,estimate]) 
    suppressWarnings({                              #for non-plot arguments in ...    
      plot3D(x@grid, col=col, ...)                  #no next method
    })   

    return(invisible())

}    #plot3D for 'InclusionZoneGrid'
) #setMethod

