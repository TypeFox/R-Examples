PlotTrajFreq <-
function( spGridDf, background = T, 
            overlay = NA, overlay.color = "white", 
            pdf = F, file.name = "output", ... )
{
  # This function is designed to display trajectory frequency map output by
  # the function RasterizeTraj.
  # Since the function RasterizeTraj outputs a RasterLayer object, this Object
  # must be converted to SpatialGridDataDataFrame Object using the 
  # as( rasterObject, "SpatialGridDataFrame" ) for example.
  #
  # Args:
  #   spGridDf: SpatialGridDataFrame Object obtened by the convertion of the 
  #             raster Object output by the RasterizeTraj function.
  #   background: Boolean: Indicates whether or not the Canada background map 
  #               should be displayed.
  #   overlay: SpatialPolygonsDataFrame
  #   overlay.color: String. sets the Poligons' color defined by the overlay
  #                  argument e.g. "blue"
  #   file.name: String: If the argument pdf is True, this argument defined
  #              the name of the output file.
  #   pdf: Boolean. Defined whether or not the output map should be saved 
  #                 in a pdf file
  # Results:
  #   Chart   
  
  if (pdf == T){
    pdf(file.name, paper="USr", height=0, width=0)
  }
  
  #get the old par configuration
  oldpar <- par(no.readonly=TRUE)
  
  # it reduces the margin's size
  par(mar=c(0,0,0,0) + 2.0)
  
  plot.add <- F
  
  # get all extra arguments
  extra.args <- list(...)
  
  # if the argument main was not defined
  if (!"main" %in% names(extra.args)) {
    extra.args$main <- NULL
  }
  
  if(background == T){
    
    bb <- bbox(spGridDf)
    
    PlotBgMap(spGridDf, xlim=bb[1,], ylim=bb[2,], axes=TRUE)
    
    # print the grid
    grid(col="white")
    
    plot.add <- T
  }
  
  #spplot(r1, add=T)
  
  grays <- colorRampPalette(c( "light green", "green", "greenyellow", "yellow", "orange", "orangered", "red"))(10)
  image( spGridDf, col=grays, breaks=(c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)), add = plot.add)
  
  legend("topleft", legend=c("0.0 - 0.1", 
                             "0.1 - 0.2",
                             "0.2 - 0.3",
                             "0.3 - 0.4",
                             "0.4 - 0.5",
                             "0.5 - 0.6",
                             "0.6 - 0.7",
                             "0.7 - 0.8",
                             "0.8 - 0.9",
                             "0.9 - 1.0")
         , fill = grays)
  
  do.call(title, extra.args)
  
  if(!missing(overlay)){
    plot(overlay, add = T, col="black", border="black")
  }
  
  # restore the par configuration
  par(oldpar)
  
  if(pdf == T){
    dev.off()
  }
}
