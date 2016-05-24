#' Maps user half degree gridded data at country level by first aggregating.
#' 
#' Maps user half degree gridded data at country level by first aggregating.
#' 
#' Aggregates half degree gridded data to countries using the option specified
#' in 'aggregateOption' then maps at a country level.
#' 
#' @param inFile either a gridascii filename or an sp SpatialGridDataFrame
#' object specifying a global half degree grid dataset, if none specified an
#' internal example data is used
#' @param aggregateOption how to aggregate the data ('sum','mean','min','max')
#' @param nameCountryColumn optional name of column containing country names
#' (used in reporting of success/failure)
#' @param suggestForFailedCodes T/F whether you want system to suggest for
#' failed codes NOT YET WORKING
#' @param projection deprecated june 2012
#' @param mapResolution options low, medium, only for projection='none'
#' initially
#' @param numCats number of categories, may be overided e.g. if catMethod
#' ='pretty'
#' @param xlim map extents c(west,east), can be overidden by mapRegion
#' @param ylim map extents c(south,north), can be overidden by mapRegion
#' @param mapRegion 'world','africa','oceania','eurasia','uk' sets map extents,
#' overrides we,ea etc.
#' @param catMethod method for categorisation of data "pretty", any vector
#' defining breaks, "fixedWidth","quantiles"
#' @param colourPalette "heat","white2Black","palette":for current palette
#' @param addLegend whether to add a legend or not T/F
#' @param lwd line width for country borders
#' @return invisibly returns a list containing the data and main options used
#' for the map, the list can be passed to \code{\link{addMapLegend}} along with
#' additional options to allow greater flexibility in legend creation.
#' @author andy south
#' @importFrom maptools readAsciiGrid
#' @seealso \code{\link{aggregateHalfDegreeGridToCountries}}
#' @keywords aplot
#' @examples
#' 
#' 
#' data(gridExData,envir=environment(),package="rworldmap")
#' gridExData <- get("gridExData")
#' mapHalfDegreeGridToCountries(gridExData)             
#' 
#' #different aggregate option
#' mapHalfDegreeGridToCountries( gridExData, aggregateOption="mean" )
#' 
#' 
#' @export mapHalfDegreeGridToCountries
`mapHalfDegreeGridToCountries` <- function( inFile=""
                          , aggregateOption="sum"  #"mean","max","min"
                          , nameCountryColumn = ""
                          , suggestForFailedCodes = FALSE 
                          , projection=NA 
                          , mapResolution="low" #options low, medium, only for projection='none' initially                                                 
                          , numCats = 7 # *may be ignored by catMethod
                           #we=-160,
                           #ea=160,
                           #so=-80,
                           #no=90,
                          , xlim=c(-160,160)
                          , ylim=c(-80,90)
                          , mapRegion = "world"   #sets map extents, overrides we,ea etc.
                          , catMethod="quantiles"   #any vector defining breaks, "fixedWidth","quantiles"
                          , colourPalette= "heat" #"heat","white2Black","palette" for current palette
                          , addLegend=TRUE
                          , lwd=0.5
                         )
   { 
                       
    #map global half degree gridded data 
    #will work on a gridascii file or a SpatialGridDataFrame or use the example data if none specified
    if ( class(inFile)=="SpatialGridDataFrame" ) 
    {
       sGDF <- inFile           
    }
    else if ( inFile == "" )
    {
       message("using example data because no file specified in mapHalfDegreeGridToCountries()\n")
       data("gridExData",envir=environment(),package="rworldmap")
       sGDF <- get("gridExData") # copying from the example data
    } else if ( is.character(inFile)) 
    {           
       if ( !file.exists(inFile) )
          {
           warning("the file: ",inFile," seems not to exist, exiting mapHalfDegreeGridToCountries()\n")
           return(FALSE)
          }
       #reading file into a SpatialGridDataFrame   
       sGDF <- readAsciiGrid(fname=inFile)    

    } else
    {
       warning(inFile," seems not to be a valid file name or a SpatialGridDataFrame, exiting mapHalfDegreeGridToCountries()\n") 
    }
         
    #further checking grid resolution
    if ( gridparameters(sGDF)$cellsize[1]!=0.5 )
        warning(inFile," seems not to be a half degree grid, in mapHalfDegreeGridToCountries()\n")
                         
     #aggregate the data to countries (pass the sGDF because already read in above)
     dF <- aggregateHalfDegreeGridToCountries(inFile=sGDF, aggregateOption=aggregateOption)
     
     mapToPlot <- joinCountryData2Map(dF, 
                          , joinCode = "UN" #options "ISO2","ISO3","FIPS","NAME", "UN" = numeric codes
                          , nameJoinColumn = "UN"
                          , nameCountryColumn = "" #there is no country column from gridascii data
                          , suggestForFailedCodes = suggestForFailedCodes 
                          , projection=projection  #"none", "EqualArea" 
                          , mapResolution=mapResolution #options low, medium, only for projection='none' initially                                                
                          )     

    ## if join has failed, then exit this function too, message from join should be enough
    #if ( class(mapToPlot)!="SpatialPolygonsDataFrame" ) return(FALSE)
     
     #call map plotting function
     mapParams <- mapCountryData( mapToPlot
                          , nameColumnToPlot = names(dF)[2]                                             
                          , numCats = numCats # *may be ignored by catMethod
                          , xlim=xlim
                          , ylim=ylim
                          , catMethod=catMethod   #any vector defining breaks, "fixedWidth","quantiles"
                          , colourPalette = colourPalette
                          , addLegend=addLegend 
                          , lwd=lwd        
                          )      

    #returning mapParams so they can be used by addMapLegend()
    invisible(mapParams)
                         
                         
    } #end of mapHalfDegreeGridToCountries()

