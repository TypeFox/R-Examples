#' function to produce pie charts on a map
#' 
#' The function will produce a map with pie charts centred on country centroids
#' (or other chosen points). The size of the circles is determined by the sum
#' of the attribute columns and each section is coloured.
#' 
#' Beware of creating plots that are difficult for the reader to interpret.
#' More than 3 or 4 categories may be too many.
#' 
#' @param dF data frame or SpatialPolygonsDataFrame
#' @param nameX name of column containing the X variable (longitude), not
#' needed if dF is a SpatialPolygonsDataFrame
#' @param nameY name of column containing the Y variable (latitude), not needed
#' if dF is a SpatialPolygonsDataFrame
#' @param nameZs name of columns containing numeric variables to determine pie
#' sections
#' @param zColours colours to apply to the pie section for each attribute
#' column
#' @param ratio the ratio of Y to N in the output map, set to 1 as default
#' @param addCatLegend whether to add a legend for categories
# removed hasn't worked for a while
# @param addSizeLegend whether to add a legend for symbol size
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
#' ## these examples repeat the same column in 'nameZs' 
#' ## to show that equal sized pies are created  
#' 
#' #mapPies( dF,nameX="LON", nameY="LAT",nameZs=c('AREA','AREA') )
#' 
#' #mapPies( dF,nameX="LON", nameY="LAT",nameZs=c('AREA','AREA')
#' #       , mapRegion='africa' )
#' 
#' mapPies( dF,nameX="LON", nameY="LAT"
#'        , nameZs=c('POP_EST','POP_EST','POP_EST','POP_EST'),mapRegion='africa' )
#'   
#' 
#' 
#' @export mapPies
`mapPies` <- function( dF
                        ,nameX="LON", nameY="LAT" 
                        ,nameZs=c(names(dF)[3],names(dF)[4])
                        ,zColours=c(1:length(nameZs))

                        ,ratio = 1
                        #,we=0, ea=0, so=0, no=0
                        ,addCatLegend = TRUE
                        #,addSizeLegend = TRUE
                        #,plotRectangles = FALSE
                        #,addTickLabels = FALSE
                        
                        ,symbolSize = 1 #multiplier relative to the default
                        ,maxZVal=NA
                        
                         , xlim=NA
                         , ylim=NA                         
                        
                         , mapRegion = "world"   #sets map extents, overrides we,ea etc.                                                    
                         , borderCol = "grey"
                         , oceanCol=NA
                         , landCol=NA
                         ,add=FALSE                        
                        
                         ,main=''    
                         ,lwd = 0.5
                        ,... )
   {                        
    functionName <- as.character(sys.call()[[1]])

    #perhaps need to replace any na's with zeroes
    #as they will be plotted the same in pies & have a problem in seq with them ?
    #replace all nas in a df
    #31/5/12 removing
    #201403 re-enabling but moving later
    #dF[is.na(dF)] <- 0
 
    #20/7/2010 changed option for region to be set from data
    #201403 should be able to remove this same as mapBars
    if ( mapRegion == 'data' ) #( (we==0 && so==0) ) # || (we==NA && so==NA)) #caused error with some data not other
    {
      xlim <- c( min(dF[,nameX], na.rm=TRUE),max(dF[,nameX], na.rm=TRUE) )
      ylim <- c( min(dF[,nameY], na.rm=TRUE),max(dF[,nameY], na.rm=TRUE) )
    }    
    
    
    #15/10/12 copying code in here from mapBubbles() to cope with sPDF
    
    #allows just a sPDF to be passed and it will get the label points, so doesn't need nameX & nameY to be specified
    if (class(dF)=="SpatialPolygonsDataFrame")
    {
      #10/10/12 moved plotting to after getting & setting centroids
      
      #22/4/10 changing this to get centroid coords from the spdf
      centroidCoords <- coordinates(dF)    
      
      
      #adding extra attribute columns to contain centroids (even though such columns may already be there)
      dF[['nameX']] <- centroidCoords[,1]
      dF[['nameY']] <- centroidCoords[,2]    
      nameX <- 'nameX'
      nameY <- 'nameY'

      if (!add) 
      {
        #use passed sPDF as the background map
        lims <- rwmNewMapPlot(mapToPlot=dF,oceanCol=oceanCol,mapRegion=mapRegion, xlim=xlim, ylim=ylim)
        #26/3/13
        xlim <- lims$xlim
        ylim <- lims$ylim
        #22/10/12 added main=main but doesn't work
        plot( dF, add=TRUE, border=borderCol, col=landCol, main=main, lwd=lwd )
      }
      
      #within this function just need the dF bit of the sPDF
      dF <- dF@data
      
    } else if (!add) 
    {
      #background map
      #these set the most common params, if user wanted finer control over map
      #they can call rwmNewMapPlot, and then call mapBubbles with add=TRUE     
      lims <- rwmNewMapPlot(mapToPlot=getMap(),oceanCol=oceanCol,mapRegion=mapRegion, xlim=xlim, ylim=ylim)
      #26/3/13
      xlim <- lims$xlim
      ylim <- lims$ylim      
      #22/10/12 added main=main but doesn't work
      plot( getMap(), add=TRUE, border=borderCol, col=landCol, main=main, lwd=lwd )
    }        
    
 
    #201403 re-enabling but moving later to avoid entanglement with sPDF
    dF[is.na(dF)] <- 0   
       
    maxSumValues <- 0
    #go through each circle to plot to find maximum value for scaling
    for (locationNum in 1:length(dF[,nameZs[1]]))
      {  
       sumValues <- sum( dF[ locationNum, nameZs ], na.rm=TRUE )
       if ( sumValues > maxSumValues ) maxSumValues <- sumValues
      }
      
    #should set radius to 5% of max extent (diam will be 10%)  
    #symbolMaxSize <- 0.05*max( ea-we, (no-so)*ratio )
    #but seemed to big ? set to 2% instead  
    #symbolMaxSize <- 0.02*max( ea-we, (no-so)*ratio )
    symbolMaxSize <- 0.02*max( xlim[2]-xlim[1], (ylim[2]-ylim[1])*ratio ) 
        
    #symbol size
    #maxZVal & symbolSize can be set by user
    #if ( is.na(maxZVal) ) maxZVal <- max( dF[,nameZSize], na.rm=TRUE )
    #4 in here is just a good sizing default found by trial & error
    #fMult = symbolSize * 4 / sqrt(maxZVal)
    #cex= fMult*sqrt(dF[,nameZSize])
    
    #so want maxSumValues to equate to maxSize (and remember they scale by square root)
    symbolScale <- symbolMaxSize / sqrt( maxSumValues )
    
    cat("symbolMaxSize=",symbolMaxSize," maxSumValues=",maxSumValues," symbolScale=",symbolScale,"\n")
    
    #for each circle to plot (row, got from num rows for first z value)
    for (locationNum in 1:length(dF[,nameZs[1]]))
      {    
       #to get an array of the values for each slice
       sliceValues <- as.numeric( dF[ locationNum, nameZs ] )
       
       #if the total of all values is 0 then skip this circle
       if (sum(sliceValues, na.rm=TRUE)==0) next
       
       #x is a cumulative list of proportions starting at 0 (i.e. 1 greater than num slices)
       cumulatProps <- c(0,cumsum(sliceValues)/sum(sliceValues, na.rm=TRUE))
       #cat("cumulative proportions", cumulatProps,"\n")
    
       #setting base radius of circles, is modified later to get aspect ratio right
       pointsInCircle = 360
       #radius <- sqrt(sum(sliceValues))*0.006
       radius <- sqrt(sum(sliceValues, na.rm=TRUE))*symbolScale
       radius <- radius*symbolSize
       
       #for each slice
       for ( sliceNum in 1:length(sliceValues) ) {
       
            #number of points on the circumference, minimum of 2
            #difference between next cumulative prop & this
            n <- max(2, floor((pointsInCircle * (cumulatProps[sliceNum+1]-cumulatProps[sliceNum]))))
            
            #to check on colours
            #cat(s,"   ", i,"   ",colours()[tc[i]],"\n")
    
            #x=radius/cos(median(LatLong$ShotLat)) determines aspect ratio of circles
            #which can change over the map ??
            #is same as the setting in plot.map above
            #previously was set to radius*1.5, but they were a bit squashed
            #P contains coordinates for the circumference bit of the slice in $x & $y
            #15/10/12 removed na.rm from seq() to correct warning
            P <- list( x= ratio * radius * cos(2*pi*seq(cumulatProps[sliceNum],cumulatProps[sliceNum+1],length=n))+ dF[ locationNum, nameX ],
                       y=         radius * sin(2*pi*seq(cumulatProps[sliceNum],cumulatProps[sliceNum+1],length=n))+ dF[ locationNum, nameY ] )
            #I wonder if this could be done with vectors rather than geometry ?
            #i.e. specifying an angle and a distance, the distance will alwasys be the same
            #and the angle will just be the proportion of the total angle
            #sure that I have methods in Java to do that
            
            #cat("slice coords", P,"\n")
    
            #plot each slice
            polygon(c(P$x,dF[ locationNum, nameX ]),c(P$y,dF[ locationNum, nameY ]),col=zColours[sliceNum]) #,col=colours()[tc[i]])
           } #end of each slice in a circle
        } #end of each circle
    
    #legend("bottomleft", select, fill=colours()[tc], cex=0.7, bg="white")
    if (addCatLegend)
        legend("bottomleft", legend=nameZs, fill=zColours, cex=0.7, bg="white")#fill=c(1:length(nameZs))
    
    #do I also want to add option for a legend showing the scaling of the symbols
    #legend(x='bottomright', legend=legendVals, pt.cex = legendSymbolSizes, pch=1, col="black", bg="white")
    #trying just a single symbol for a start
    #legendVals = maxSumValues
    #getting the size of the max symbol in terms of pt.cex, could be very tricky
    #might need to just draw a circle & text on the map, but getting that right will also be tricky
    #legendSymbolSizes = 
    #legend(x='bottomright', legend=legendVals, pt.cex = legendSymbolSizes, pch=1, col="black", bg="white")    
    #P <- list( x= ratio * radius * cos(2*pi*seq(cumulatProps[sliceNum],cumulatProps[sliceNum+1],length=n))+ dF[ locationNum, nameX ],
    #           y=         radius * sin(2*pi*seq(cumulatProps[sliceNum],cumulatProps[sliceNum+1],length=n))+ dF[ locationNum, nameY ] )
    
    #ratio <- 2
    #radius <- 1
    #to create a whole circle equivalent to the max symbol size
    radius <- symbolMaxSize*symbolSize
    
    #par('usr'), returns extents of plot, so can use to put components in specific places relative to it
    #[1] -7.161063 -1.689103 48.424150 50.899516
    
    plotExtents <- par('usr')
    plotS <- plotExtents[3]  
    plotW <- plotExtents[1] 
    plotN <- plotExtents[4]  
    plotE <- plotExtents[2] 
    #top left
    centreE <- plotW + radius*2*ratio
    centreN <- plotN - radius*2
    #bottom right
    #centreE <- plotE - radius*2*ratio
    #centreN <- plotS + radius*2
    t <- seq(0,2*pi,length=100)
    P <- list( x= ratio * radius * cos(t)+ centreE,
               y=         radius * sin(t)+ centreN )
    
    str(P)
    
    # LEGEND FOR THE SIZE OF CIRCLES
    
    #could also try to create a box to include the circle & text
    #if (addSizeLegend)           
    #    polygon(P$x,P$y,border="red");#,col="red")  without the col it wouldn't be filled
                                            #for some reason didn't fill properly anyway ?
    #to put text in centre of circle ...
    #text(centreE,centreN,labels=maxSumValues) 
    # could try to put below circle                                                  
    #text(centreE,centreN-radius*1.5,labels=maxSumValues) 
    #text(centreE,centreN-radius*1.5,labels=signif(maxSumValues,digits=4)) 
    
    } # end of mapPies



#######################
#testing the function
    
#dF <- getMap()@data   
#mapPies( dF,nameX="LON", nameY="LAT",nameZs=c('POP_EST','AREA') )
#mapPies( dF,nameX="LON", nameY="LAT",nameZs=c('AREA','AREA') )
#mapPies( dF,nameX="LON", nameY="LAT",nameZs=c('AREA','AREA','AREA'),mapRegion='africa' )  
