#' to print country labels on a world map
#' 
#' Given no arguments it will print country names stored in the 'NAME' column
#' of \code{\link{getMap}} onto an existing map at the centroids of each
#' country polygon, stored in the 'LAT' and 'LON' columns.  Alternatively the
#' user can specifiy a data frame or SpatialPolygonsDataFrame in which case
#' they need to define the column containing the country names
#' (nameCountryColumn) and optionally a 2nd attribute column to print
#' (nameColumnToPlot).  First you need to create a map plot, for example using
#' \code{\link{mapCountryData}} or \code{\link{mapBubbles}}.
#' 
#' 
#' @param dF dataframe or SpatialPolygonsDataFrame
#' @param nameCountryColumn name of column containing country names to be
#' printed on the map (could also be set to any other column in the dataframe)
#' @param nameX name of column containing the X variable (longitude), not
#' needed if dF is a SpatialPolygonsDataFrame
#' @param nameY name of column containing the Y variable (lattitude), not
#' needed if dF is a SpatialPolygonsDataFrame
#' @param nameColumnToPlot name of an attribute column in the data frame the
#' value of which will be appended to the country names
#' @param col colour for labels, default 'grey', can be e.g.
#' rgb(1,1,0,alpha=0.5)
#' @param cex sizing of labels, default = 0.8
#' @param \dots other parameters that can be passed to text(), e.g. pos=4 to
#' right, (1=below, 2=left, 3=above)
#' @return nothing
#' @author andy south
#' @seealso \code{\link{identifyCountries}}
#' @keywords dplot
#' @examples
#' 
#' mapCountryData()
#' labelCountries()
#' 
#' labelCountries(nameColumnToPlot = "POP_EST")
#' 
#' 
#' @export labelCountries
labelCountries <- function( dF=""
                          , nameCountryColumn="NAME" #?maybe change this to nameCountry
                          , nameX="LON"
                          , nameY="LAT"
                          , nameColumnToPlot=""
                          , col = 'grey'
                          , cex = 0.8
                          , ...){
  
#27/10/12 initially make the first part identical to identifyCountries
#can then put common bits into it's own function  
  
  #^^!! start of bit initially copied from identifyCountries()
  
  #also possible to get centroids from a sPDF
  if (class(dF)=="SpatialPolygonsDataFrame")
  {
    #!9/10/12 BUG correction, don't get coords from internal map if an sPDF is passed
    #centroidCoords <- coordinates(getMap())
    centroidCoords <- coordinates(dF)    
    #within this function just need the dF bit of the sPDF
    dF2 <- dF@data
    #adding extra attribute columns to contain centroids (even though such columns may already be there)
    dF2[['nameX']] <- centroidCoords[,1]
    dF2[['nameY']] <- centroidCoords[,2]    
    nameX <- 'nameX'
    nameY <- 'nameY'    
  } else
    #this assumes that the dF already has columns for lat & lon
    if (class(dF)=="data.frame")
    {
      #this assumes that nameX, nameY & nameCountryColumn columns have been passed correctly
      dF2 <- dF
    } else
    {
      #if no object is passed use the internal map
      dF2 <- getMap()@data
      nameX <- 'LON'
      nameY <- 'LAT'    
    }
  
  labels <- dF2[[nameCountryColumn]]
  
  #if an attribute column name is passed paste it's value onto end of country label
  if ( nameColumnToPlot != "" ) labels <- paste(labels,dF2[[nameColumnToPlot]])  
  
  #^^!! end of bit initially copied from identifyCountries()
  

  #plotting the labels
  text( dF2[[nameX]], dF2[[nameY]], labels=labels, col=col, cex=cex, ... )
  
  
}
