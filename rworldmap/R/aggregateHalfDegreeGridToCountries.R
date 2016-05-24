#' Aggregates global half degree gridded data to countries
#' 
#' Aggregates global half degree gridded data to countries (options for sum,
#' mean, min, max ). Uses a very simple grid map defining a single country
#' identity for each half degree cell.  (other more sophisticated approaches
#' dividing cells between multiple countries will be investigated in future).
#' The country identity at each cell is specified in
#' data(gridCountriesDegreesHalf).
#' 
#' 
#' @param inFile either a gridascii filename or an sp SpatialGridDataFrame
#' object specifying a global half degree grid dataset
#' @param aggregateOption how to aggregate the data ('sum','mean','min','max')
#' @return a dataframe with 2 columns : numeric country codes and the
#' aggregated value for each country
#' @author andy south
#' @importFrom maptools readAsciiGrid
#' @seealso \code{\link{mapHalfDegreeGridToCountries}}
#' @keywords dplot
#' @examples
#' 
#' 
#' data(gridExData,envir=environment(),package="rworldmap")
#' gridExData <- get("gridExData")
#' #aggregating the gridded data to countries
#' dF <- aggregateHalfDegreeGridToCountries(gridExData)
#' #joining the aggregated data to a country map
#' sPDF <- joinCountryData2Map(dF, nameJoinColumn='UN', joinCode='UN')
#' #plotting the map
#' mapCountryData(sPDF,nameColumnToPlot='sum_pa2000.asc')
#' 
#' 
#' @export aggregateHalfDegreeGridToCountries
`aggregateHalfDegreeGridToCountries` <-
function( inFile=""
                         ,aggregateOption="sum"  #"mean","max","min"
                         )
   {
    
    #a function to aggregate half degree grid data to country level
    #based upon a grid country file obtained from IIASA
    #returns a dataframe with numeric country code & aggregated values
    #can do sum, mean, max or min        

    #added an option to work on an existing loaded spatialGridDataFrame

    if ( is.character(inFile) )
    {
       if ( !file.exists(inFile) )
          {
           warning("the file: ",inFile," seems not to exist, exiting aggregateHalfDegreeGridToCountries()\n")
           return(FALSE)
          }
       #reading file into a SpatialGridDataFrame   
       sGDF <- readAsciiGrid(fname=inFile)
            
    } else if ( class(inFile)=="SpatialGridDataFrame" ) 
    {
       #if its already a SpatialGridDataFrame just copy it
       #! 6/3/09 this allows for the potential for multiple attribute columns 
       #! which is not coped with below
       sGDF <- inFile
    } else
    {
       warning(inFile," seems not to be a valid file name or a SpatialGridDataFrame, exiting aggregateHalfDegreeGridToCountries()\n") 
    }
         
    #further checking grid resolution
    if ( gridparameters(sGDF)$cellsize[1]!=0.5 )
        warning(inFile," seems not to be a half degree grid, in aggregateHalfDegreeGridToCountries()\n")
     
    
    #25/03/2013 replacing grid file with updated countries
    data(gridCountriesDegreesHalf,envir=environment(),package="rworldmap")
    sGDFcountries <- get("gridCountriesDegreesHalf")
    
    #getting the names of the columns containing the data
    attrNameGrid <- names(sGDF)[1]
    attrNameGridCountries <- names(sGDFcountries)[1]
    
    dF <- data.frame(attribute = sGDF[[attrNameGrid]]
                    ,UN = sGDFcountries[[attrNameGridCountries]])
    
    #4 aggregate cell values by numeric country code
    #! later offer option for user to specify which
    dFbyCountry <- aggregate(dF$attribute
                         , by=list(UN = dF$UN)
                         , FUN = aggregateOption
                         , na.rm=TRUE )
                         #, mean, na.rm=TRUE )
                         
    #renaming the aggregated column name from x, with filename if the input was a file (if a sGDF it causes error use names instead)
    if ( is.character(inFile) )
    {
       names(dFbyCountry)[2] <- paste(aggregateOption,"_",basename(inFile), sep='')            
    } else #if inFile is a sGDF 
    {
       #!6/3/09 this causes an error if the SGDF contained multiple attribute columns
       names(dFbyCountry)[2] <- paste(aggregateOption,"_",names(inFile), sep='') 
    }
    
    
    return(dFbyCountry)
      
    } #end of aggregateHalfDegreeGridToCountries()

