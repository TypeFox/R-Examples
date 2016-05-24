# Function name: nSlidesANDnTuple 
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


nSlidesANDnTuple <- function(nLevels,nSlides,nTuple,bTwoColorArray)
{
   #This is used by designGG main function
   nConditions <- 1
    for(i in 1:length(nLevels)) nConditions <- nConditions*nLevels[i]
    if (is.null(nSlides) & is.null(nTuple))
    {
        stop( "Either nSlides ( number of slides ) or nTuple ( number of RILs
                per condition is required for this algorithm.")
    }
    if ( !is.null(nSlides) & is.null(nTuple) ) {
        if ( bTwoColorArray )
        {
            nTuple <- 2*nSlides/nConditions
        }
        else
        {
            nTuple <- nSlides/nConditions
        }
        if ( nTuple < 1 )
        {stop("The number slides is too less to perform the experiment.")}
    }
    if ( is.null(nSlides) & !is.null(nTuple) ) {
        if ( nTuple < 1 )
        { stop ( "The number of slides is too small to perform the experiment." )}
        if ( bTwoColorArray )
        {   if ( round((nConditions*nTuple/2)%%(floor(nTuple/2)),digits=2) == 0 )
                {
                    nSlides <- nConditions*nTuple/2
                } else
                {
                    nSlides <- ceiling(nConditions * nTuple/2)
                }
        }else
        {
            if ( round(floor(nConditions*nTuple)-(nConditions*nTuple),digits=3) == 0 )
                {
                    nSlides <- nConditions*nTuple
                } else
                {
                    nSlides <- ceiling(nConditions * nTuple)
                }
        }
    }
      return(list(nSlides=nSlides, nTuple=nTuple))
  }