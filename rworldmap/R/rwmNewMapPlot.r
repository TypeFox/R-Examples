#' Internal function to set up an existing device for plotting maps
#' 
#' Sets the region, aspect and ocean colour for a new map plot
#' 
#' Called by mapCountryData() and mapGriddedData()
#' 
#' @param mapToPlot the worldmap to be plotted
#' @param oceanCol a colour for the ocean
#' @param mapRegion a string specifying the map region, see setMapExtents()
#' @param xlim map extents c(west,east), can be overidden by mapRegion
#' @param ylim map extents c(south,north), can be overidden by mapRegion
#' @param aspect aspect for the map, defaults to 1, if set to 'variable' uses
#' same default as plot.Spatial in sp
#' @return a dataframe containing xlim and ylim
#' @author andy south
#' @keywords misc
#' @export rwmNewMapPlot
rwmNewMapPlot <- function(mapToPlot=getMap(),
         oceanCol=NA,
         mapRegion="world",
         xlim=NA,                  
         ylim=NA,
         aspect=1){

  #browser()

  ## setting map extents if a mapRegion has been specified
  if (mapRegion != "world"){
    dFwesn <- setMapExtents(mapRegion)
    xlim <- c(dFwesn$we, dFwesn$ea)
    ylim <- c(dFwesn$so, dFwesn$no)
  }

  #getting xlim & ylim from bbox of map if they haven't been specified
  if (length(xlim)<2) xlim <- bbox(mapToPlot)['x',]  
  if (length(ylim)<2) ylim <- bbox(mapToPlot)['y',]
  
  
  plot.new()

  #replicate behaviour of plot.Spatial in the sp package regarding aspect
  #if the map is unprojected the aspect is set based upon the mean y coord
  #only if region not 'world'
  if (aspect == 'variable' & mapRegion != "world")
        aspect <- ifelse(is.na(proj4string(mapToPlot)) || is.projected(mapToPlot),
            1, 1/cos((mean(ylim) * pi)/180))

  plot.window(xlim=xlim,ylim=ylim,asp=aspect)#,xaxs='i',yaxs='i')#,bg=oceanCol,xpd=NA)
  
  #rect(xlim[1],ylim[1],xlim[2],ylim[2],col=oceanCol,border=oceanCol)
  #making the rectangle as big as the whole map should ensure it fills
  rect(mapToPlot@bbox[1],mapToPlot@bbox[2],mapToPlot@bbox[3],mapToPlot@bbox[4],col=oceanCol,border=oceanCol)

  #26/3/13 returning xlim & ylim so it can be used elsewhere
  lims <- data.frame(xlim,ylim)
  invisible(lims)
  
  
}
