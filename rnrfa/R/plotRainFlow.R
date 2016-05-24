#' Plot rainfall and flow for a given station
#'
#' @description This function retrieves rainfall and flow time series for a given catchment, divides the flow by the catchment area and converts it to mm/day to that it can be comparable with the rainfall (mm/month). Finally it generates a plots combining rainfall and flow information.
#'
#' @param id Station identification number
#'
#' @return Plot rainfall and flow for a given station
#'
#' @examples
#' #' # plotRainFlow(id = 54022)

plotRainFlow <- function(id){

  # Retrieve area (in Km2) from the catalogue
  meta <- catalogue(columnName="id", columnValue = id)
  area <- as.numeric(as.character(meta$catchmentArea))

  # Retrieve rainfall data for station 54022
  rain <- getTS(id, type = "cmr")

  # Retrieve flow data for station 54022
  flow <- getTS(id, type = "gdf")
  convertedFlow <- ((flow*1000)/(area*1000000))*86400 # mm/day

  proportion <- ceiling((max(convertedFlow, na.rm = T) -
                           min(convertedFlow, na.rm = T))/3)

  # opar <- par()
  par(mar=c(2,4,3,4))

  plot(convertedFlow, ylim =c(-proportion/2,max(convertedFlow)+proportion),
       main=meta$name, xlab="",ylab="Flow [mm/d]", lwd=0.5)

  # Add precipitation to the top
  par(bty="n", new=T)
  plot(rain, type="h",
       ylim=rev(range(rain)*5), # downward bars
       yaxt="n", xaxt="n", ann=F, # do not plot x and y axis
       lwd=0.5, col="deepskyblue3" ) # suggested cosmetics

  # add right axis (4) to describe P
  axis(4, pretty(range(rain))[c(2,4,6,8)],
       col.axis="deepskyblue3", col="deepskyblue3", las=1, cex.axis=0.8 )
  mtext("Rain [mm/month]", side = 4, line = 3, cex = par("cex.lab"),
        col="deepskyblue3", adj = 1)
  # reset border and overlay
  par(bty="o", new=F)

  legend("bottom",
         c("GDF", "CMR"),
         horiz = TRUE,
         y.intersp=1,
         bty = "n",
         lwd =c(3, 2),
         col = c("black","deepskyblue3"),
         cex = 1)

  # par() <- opar

}
