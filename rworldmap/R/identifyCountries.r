#' a function that will print country name and attribute values when a user
#' clicks on the map
#' 
#' An interactive function that will print on a map the nearest country name to
#' a user mouse click.  The user can specify nothing and the function will use
#' a map from the package. Alternatively the user can specifiy a data frame or
#' SpatialPolygonsDataFrame in which case they need to define the column
#' containing the country names (nameCountryColumn) and optionally a 2nd
#' attribute column to print (nameColumnToPlot).
#' 
#' Uses the identify() function, which waits for the user to click on the map,
#' and stops when the user right clicks and selects 'stop'.
#' 
#' It uses country centroids, and will give a warning if one is too far away
#' (default value of 0.25 inches).
#' 
#' @param dF data frame or SpatialPolygonsDataFrame
#' @param nameCountryColumn name of column containing country names to be
#' printed on the map (could also be set to any other attribute the user wants
#' to query)
#' @param nameX name of column containing the X variable (longitude), not
#' needed if dF is a SpatialPolygonsDataFrame
#' @param nameY name of column containing the Y variable (lattitude), not
#' needed if dF is a SpatialPolygonsDataFrame
#' @param nameColumnToPlot name of an attribute column in the data frame the
#' value of which will be appended to the country name when it is printed
#' @param plotSelected if set to TRUE a blue outline will be printed around the
#' countries selected when the selection process is finished
#' @param \dots other parameters that can be passed to identify()
#' @return a vector of the indices of the countries selected
#' @author andy south
#' @seealso identify() \code{\link{labelCountries}}
#' @keywords dplot
#' @examples
#' 
#' #mapCountryData()
#' #identifyCountries()
#' 
#' #identifyCountries(nameColumnToPlot = "POP_EST", plotSelected = TRUE)
#' 
#' @export identifyCountries
identifyCountries <- function(dF=""
                             ,nameCountryColumn="NAME"
                             ,nameX="LON"
                             ,nameY="LAT"
                             ,nameColumnToPlot=""
                             ,plotSelected=FALSE
                             ,...){


#also possible to get centroids from a sPDF
if (class(dF)=="SpatialPolygonsDataFrame")
   {
    #get coords from sPDF if it is passed
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

selectedCountryIndices <- identify(x=dF2[[nameX]], y=dF2[[nameY]], labels=labels,...)

#allowing plotting of the boundaries of the selected countries
#this is really just an initial test of something I may develop further later
#!this will only work if the internal or a passed map is used
#!the plots only appear after the locator is stopped
if (plotSelected & length(selectedCountryIndices)>0 ) 
   {
    if (class(dF)=="SpatialPolygonsDataFrame")
       {
        plot(dF[selectedCountryIndices,],border='blue',add=TRUE)
       } else
        #!warning if a dF in a different order to getMap() is passed this will plot wrong countries
        plot(getMap()[selectedCountryIndices,],border='blue',add=TRUE)    
   }



#return the indices, may be useful later
invisible(selectedCountryIndices)
                           
} 
