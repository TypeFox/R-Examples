#' Example code for plot creation
#' 
#' Example code to demonstrate creation of a series of plots
#' 
#' 
#' @author andy south
#' @keywords aplot
#' @export rworldmapExamples
rworldmapExamples <- function()
   {  
    #displays examples - i should probably put somewhere else
    
    #get the data
    data("countryExData",envir=environment(),package="rworldmap")
    countryExData<-get("countryExData")
   
    #should check operating system for windows()
    #or just remove
    #windows()
   
    #this prompts user for keypress between plots 
    par(ask = TRUE)
      

    #1
    #joining the example data to a map
    sPDF <- joinCountryData2Map(countryExData
              , joinCode = "ISO3"
              , nameJoinColumn = "ISO3V10"
              )
    mapCountryData( sPDF
                  , nameColumnToPlot="EPI" 
                  )                 


    #2
    #joining the example data to a map
    sPDF <- joinCountryData2Map(countryExData
              , joinCode = "ISO3"
              , nameJoinColumn = "ISO3V10"
              )
    mapCountryData( sPDF
                  , nameColumnToPlot="ENVHEALTH"
                  , catMethod="quantiles"
                  , numCats = 10 
                  )
    #mtext('mapCountryData( nameColumnToPlot="ENVHEALTH", joinCode = "ISO3"
    #, nameJoinColumn = "ISO3V10", catMethod="quantiles", numCats = 10)',line=-1)                  
                  
    #windows()
    #3
    sPDF <- joinCountryData2Map(countryExData
          , joinCode = "ISO3"
          , nameJoinColumn = "ISO3V10"
          )
    mapCountryData( sPDF
                  , nameColumnToPlot="BIODIVERSITY"
                  )            
                 
                  
    #windows()
    #aggregating gridded data to country level
    #with no file specified it uses internal example data
    #5
    mapHalfDegreeGridToCountries( )              
    mtext('mapHalfDegreeGridToCountries( )',outer=TRUE,line=-1)
    
    #windows()
    #different agggregate option
    #6
    mapHalfDegreeGridToCountries( aggregateOption="mean" )
    mtext('mapHalfDegreeGridToCountries( aggregateOption="mean" )',outer=TRUE,line=-1)
        
    #windows()
    #no parameters : default dataset
    #7
    mapGriddedData()
    mtext('mapGriddedData()')
        
    #windows()
    #adding histogram
    #8
    #mapGridAscii(catMethod="logFixedWidth",addHist=TRUE)    
    #mtext('mapGridAscii(catMethod="logFixedWidth",addHist=TRUE)',outer=TRUE,line=-1)
        
    #windows() 
    #9
    mapGriddedData(mapRegion="africa") 
    mtext('mapGriddedData(mapRegion="africa")',outer=TRUE,line=-1)
            
    #switching off prompting user for keypress between plots 
    par(ask = FALSE)     
                  
    } #end of rWorldMapExamples

