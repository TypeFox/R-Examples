#' function to produce bar plots on a map
#' 
#' The function will produce a map with bars centred on country centroids (or
#' other chosen points). The length of the bars is determined by the sum of the
#' attribute columns and each section is coloured.
#' 
#' Horizontal or vertical bars can be achieved by using the barOrient argument
#' 'horiz' or 'vert'.
#' 
#' @param dF data frame or SpatialPolygonsDataFrame
#' @param nameX name of column containing the X variable (longitude), not
#' needed if dF is a SpatialPolygonsDataFrame
#' @param nameY name of column containing the Y variable (lattitude), not
#' needed if dF is a SpatialPolygonsDataFrame
#' @param nameZs name of columns containing numeric variables to determine bar
#' sections
#' @param zColours colours to apply to the bar section for each attribute
#' column
#' @param barWidth multiple for the width of bar symbols, relative to barOrient
#' see below
#' @param barOrient orientation of bars, options 'horiz' and 'vert'
#' @param barRelative default is TRUE, each variable (column) is scaled to it's
#' maximum value
#' @param ratio the ratio of Y to N in the output map, set to 1 as default
#' @param addCatLegend whether to add a legend for categories
#' @param addSizeLegend whether to add a legend for symbol size
#' @param symbolSize multiplier of default symbol size
#' @param maxZVal the attribute value corresponding to the maximum symbol size,
#' this can be used to set the scaling the same between multiple plots
#' @param xlim map extents c(west,east), can be overidden by mapRegion
#' @param ylim map extents c(south,north), can be overidden by mapRegion
#' @param mapRegion a country name from getMap()[['NAME']] or
#' 'world','africa','oceania','eurasia','uk' sets map extents, overrides
#' xlim,ylim
#' @param borderCol the colour for country borders
#' @param oceanCol a colour for the ocean
#' @param landCol a colour to fill countries
#' @param add whether to add the symbols to an existing map, TRUE/FALSE
#' @param main title for the map
#' @param lwd line width for country borders
#' @param lwdSymbols line width for symbols
#' @param \dots any extra arguments to points()
#' @return currently doesn't return anything
#' @author andy south
#' @keywords aplot
#' @examples
#' 
#' 
#' #getting example data
#' dF <- getMap()@@data 
#'    
#' mapBars( dF,nameX="LON", nameY="LAT",nameZs=c('POP_EST','GDP_MD_EST') )
#' mapBars( dF,nameX="LON", nameY="LAT",nameZs=c('POP_EST','GDP_MD_EST'), mapRegion='africa' )
#' mapBars( dF,nameX="LON", nameY="LAT",nameZs=c('POP_EST','GDP_MD_EST'), 
#'  mapRegion='africa', symbolSize=20 )
#' mapBars( dF,nameX="LON", nameY="LAT",nameZs=c('POP_EST','GDP_MD_EST'), mapRegion='africa', 
#'  symbolSize=20, barOrient = 'horiz' )
#' 
#' 
#' # this does work too
#' #mapBars( dF,nameX="LON", nameY="LAT" 
#' #       , nameZs=c('POP_EST','GDP_MD_EST')
#' #       , mapRegion='africa'
#' #       , symbolSize=4 )       
#' 
#'   
#' 
#' 
#' @export mapBars
`mapBars` <- function( dF = ""
                      , nameX="longitude", nameY="latitude" 
                      , nameZs=c(names(dF)[3],names(dF)[4])
                      , zColours=c(1:length(nameZs))
                      , barWidth = 1
                      , barOrient = 'vert' ## orientation of bars 'vert' as default or 'horiz'
                      , barRelative = TRUE

                      , ratio = 1
                        #,we=0, ea=0, so=0, no=0
                      , addCatLegend = TRUE
                      , addSizeLegend = TRUE
                        
                      , symbolSize = 1 #multiplier relative to the default
                      , maxZVal=NA

                      , xlim=NA
                      , ylim=NA   
                       
                      , mapRegion = "world"   #sets map extents, overrides we,ea etc.                                                    
                      , borderCol = "grey"
                      , oceanCol=NA
                      , landCol=NA
                      , add=FALSE                        

                      , main=''
                      , lwd=0.5
                      , lwdSymbols=1
                      , ... )
   { 
  
    functionName <- as.character(sys.call()[[1]])

    #for example data need to put in here before example dF loaded
    if (length(dF)==1 && dF == "")
    {
      nameZs <- c('POP_EST','GDP_MD_EST')
    }
    
    
    #2013 refactoring
    #this returns either a dF or sPDF
    dF <- rwmCheckAndLoadInput( dF, inputNeeded = "sPDF or dF", callingFunction=functionName ) 
    
    #if sPDF
    #  sPDF <- dF
    #  dF[nameX & nameY] <- coordinates(SPDF)
    #  dF <- dF@data
    
    #else if dF
    #  xlimylim <- max dF[nameX & nameY]
    #  sPDF <- getMap()
    
    # *shared*
    # plot map using sPDF
    # do bars using dF
    
    
    #if rwmCheckAndLoadInput returns a sPDF get the dF bit add columns for centroid coords & set nameX & nameY
    if ( class(dF)=="SpatialPolygonsDataFrame" ) #################################
    {
      #copying map to sPDF to use later
      sPDF <- dF
      
      nameX <- "rwmX"
      nameY <- "rwmY"
      coords <- coordinates(dF)
      #fill columns in dF with centroid coords
      dF[[nameX]] <- coords[,1]
      dF[[nameY]] <- coords[,2]
      #dF bit to be used for bars
      dF <- dF@data
      
    } else if( class(dF)=="data.frame"  ) #######################################
    {   
      #to be used for background map if !add
      sPDF <- getMap()
      
    } else
    {
      stop(functionName," requires a dataFrame or spatialPolygonsDataFrame for the first argument or dF=\n")
      return(FALSE)       
    }
    
    
    #debugging
    #browser()
    

    #background map
    #if user wants finer control they can call rwmNewMapPlot, and then this with add=TRUE 
    if (!add) 
    {      
      lims <- rwmNewMapPlot(sPDF, oceanCol=oceanCol, mapRegion=mapRegion, xlim=xlim, ylim=ylim)
      
      xlim <- lims$xlim #!!! these lims are used later to set symbol sizes
      ylim <- lims$ylim
      plot( sPDF, add=TRUE, border=borderCol, col=landCol, lwd=lwd )
    } #end of if (!add)    

    
    #**BEWARE what happens with symbolMaxSize if add=TRUE ???
    
    #Warning message:
    #  In max(xlim[2] - xlim[1], (ylim[2] - ylim[1]) * ratio) :
    #  no non-missing arguments to max; returning -Inf
    
    #browser()
    
    #1/7/13 adding a relative option so that all bars can be scaled 0-1
    #partly to make it easier to produce an example plot
    if (barRelative)
    {
      for( numZ in 1:length(nameZs))
      {
        #TEMPORARY FIX TO REPLACE -99 with NA for pop & gdp
        #if ( length(which(dF[nameZs][numZ]=="-99") ))
        #  dF[nameZs][numZ][ which(dF[nameZs][numZ]=="-99"),1 ] <- NA  
        
        dF[nameZs][numZ] <- dF[nameZs][numZ] / max(dF[nameZs][numZ],na.rm=TRUE)
        
      }
    }
    
    #browser()
    
    maxSumValues <- 0
    #go through each circle to plot to find maximum value for scaling
    for (locationNum in 1:length(dF[,nameZs[1]]))
      {  
       sumValues <- sum( dF[ locationNum, nameZs ], na.rm=TRUE )
       if ( sumValues > maxSumValues ) maxSumValues <- sumValues
      }
      
    
    #browser()    
    
    #set symbolMaxSize to 2% of max extent 
    symbolMaxSize <- 0.02*max( xlim[2]-xlim[1], (ylim[2]-ylim[1])*ratio, na.rm=TRUE )    
        
    #symbol size
    #maxZVal & symbolSize can be set by user
    #if ( is.na(maxZVal) ) maxZVal <- max( dF[,nameZSize], na.rm=TRUE )
    #4 in here is just a good sizing default found by trial & error
    #fMult = symbolSize * 4 / sqrt(maxZVal)
    #cex= fMult*sqrt(dF[,nameZSize])
    
    #so want maxSumValues to equate to maxSize
    symbolScale <- symbolMaxSize / maxSumValues 
    
    cat("symbolMaxSize=",symbolMaxSize," maxSumValues=",maxSumValues," symbolScale=",symbolScale,"\n")
    
    #for each location (row, got from num rows for first z value)
    for (locationNum in 1:length(dF[,nameZs[1]]))
      {    
       #to get an array of the values for each slice
       sliceValues <- as.numeric( dF[ locationNum, nameZs ] )
       
       #if the total of all values is 0 then skip this circle
       if (sum(sliceValues, na.rm=TRUE)==0) next
       
       #x is a cumulative list of proportions starting at 0 (i.e. 1 greater than num slices)
       cumulatProps <- c(0,cumsum(sliceValues)/sum(sliceValues, na.rm=TRUE))
       #cat("cumulative proportions", cumulatProps,"\n")
    
       #radius <- sqrt(sum(sliceValues, na.rm=TRUE))*symbolScale
       #1/7/2013 removing sqrt
       radius <- sum(sliceValues, na.rm=TRUE)*symbolScale       
       
       radius <- radius*symbolSize
       
       #for each slice
       for ( sliceNum in 1:length(sliceValues) ) {
       
            #rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,col = NA, border = NULL, lty = par("lty"), lwd = par("lwd")
            
            if ( barOrient == 'horiz' )
               {
                #cat('horiz')
                xleft <- dF[ locationNum, nameX ] + ( radius * cumulatProps[sliceNum] )
                ybottom <- dF[ locationNum, nameY ]  
                xright <- dF[ locationNum, nameX ] + ( radius * cumulatProps[sliceNum+1] ) 
                ytop <- dF[ locationNum, nameY ] + barWidth  
               } else
               {
                #cat('vert')
                xleft <- dF[ locationNum, nameX ] 
                ybottom <- dF[ locationNum, nameY ] + ( radius * cumulatProps[sliceNum] ) 
                xright <- dF[ locationNum, nameX ] + barWidth 
                ytop <- dF[ locationNum, nameY ] + ( radius * cumulatProps[sliceNum+1] )  
               }                         
            
            rect( xleft, ybottom, xright, ytop, col=zColours[sliceNum],lwd=lwdSymbols )
            #number of points on the circumference, minimum of 2
            #difference between next cumulative prop & this
 
            #cat("slice coords", P,"\n")
    
            #plot each slice
            #polygon(c(P$x,dF[ locationNum, nameX ]),c(P$y,dF[ locationNum, nameY ]),col=zColours[sliceNum]) #,col=colours()[tc[i]])
           } #end of each slice in a circle
        } #end of each circle
    
    #legend("bottomleft", select, fill=colours()[tc], cex=0.7, bg="white")
    if (addCatLegend)
        legend("bottomleft", legend=nameZs, fill=zColours, cex=0.7, bg="white")#fill=c(1:length(nameZs))
    
    #do I also want to add option for a legend showing the scaling of the symbols
    #legend(x='bottomright', legend=legendVals, pt.cex = legendSymbolSizes, pch=1, col="black", bg="white")
    
    
    } # end of mapBars

    


  
