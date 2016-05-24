#rwmCheckAndLoadInput
#andy south 28/6/2013

#to check and load the input data for any rwm function
#if "" is passed it should load example data (which later could be specific to the calling function)

#I could set up two main alternatives
#if requireSPDF
#   =="" : get example SPDF
#   is SPDF : use it
#   not SPDF : give warning about using joinCountry2Map
#if !requireSPDF
#   =="" : get example dF
#   is SPDF : use dF bit
#   is DF : use it
#   is something else : error

#then maybe later put nameColumnToPlot checking in here
#BUT PROB better to put in another function



#' internal function to check and load input data to mapping functions
#' 
#' Internal function checking and loading dFs or sPDFs to
#' \code{\link{mapCountryData}}, \code{\link{mapPolys}}, \code{\link{mapPies}},
#' \code{\link{mapBubbles}}, \code{\link{mapBars}}.
#' 
#' a rworldmap internal function, unlikely to be of use to users
#' 
#' @param inputData a dF, sPDF or "", for latter an internal example data is
#' used
#' @param inputNeeded "sPDF", "sPDF or dF", "dF"
#' @param callingFunction optional : name of the calling function
#' @return invisibly returns a dF or sPDF
#' @author andy south
#' @keywords aplot
#' @export rwmCheckAndLoadInput
rwmCheckAndLoadInput <- function(
    inputData =         ""
  #, nameColumnToPlot =  ""
  #, requireSPDF =       TRUE
  , inputNeeded =       "sPDF" #options "sPDF", "dF", "sPDF or dF"    
  , callingFunction = "" #currently optional may be useful later  
){
  
  functionName <- as.character(sys.call()[[1]])
  
  #message(paste("In", functionName, "called by", callingFunction))  
  
  #browser()  #n to enter the step through debugger, Q to exit
  
  #if (requireSPDF)
  if (inputNeeded == "sPDF")
  {
    if ( class(inputData)=="SpatialPolygonsDataFrame" ) 
    {
      ## checking if there is any data in the dataFrame
      if ( length(inputData@data[,1]) < 1 ){
        stop("seems to be no data in your chosen input in ",functionName, "from", callingFunction) 
        return(FALSE)
      } 
    } else if ( length(inputData)==1 && inputData == "" ) 
    {
      message(paste("using example data because no file specified in",functionName))
      
      inputData <- getMap(resolution="coarse")
      
      ## also setting a default nameColumnToPlot if it isn't set
      #can't have this here because not passed
      #if ( nameColumnToPlot == "" ) nameColumnToPlot <- "POP_EST" #
    } else 
    {
      stop(callingFunction," requires a SpatialPolygonsDataFrame object created by the joinCountryData2Map() or joinData2Map() functions \n")
      return(FALSE) 
    }
     
  } else if (inputNeeded == "sPDF or dF") 
  {
    #   =="" : get example sPDF
    #   is SPDF : return it
    #   is DF : use it
    #   is something else : error    
    if ( length(inputData)==1 && inputData == "" ) 
    {
      message(paste("using example data because no file specified in",functionName))
      
      inputData <- getMap(resolution="coarse")
      
      ## also setting a default nameColumnToPlot if it isn't set
      #can't have this here because not passed
      #if ( nameColumnToPlot == "" ) nameColumnToPlot <- "POP_EST" #
    } else if ( class(inputData)=="SpatialPolygonsDataFrame" ) 
    {
      ## checking if there is any data in the dataFrame
      if ( length(inputData@data[,1]) < 1 ){
        stop("seems to be no data in your chosen input in ",functionName, "from", callingFunction) 
        return(FALSE)
      } 
    } else if ( class(inputData)=="data.frame" ) 
    {
      ## checking if there is any data in the dataFrame
      if ( length(inputData[,1]) < 1 ){
        stop("seems to be no data in your chosen input in ",functionName, "from", callingFunction) 
        return(FALSE)
      }       
    } else  
    {
      stop(callingFunction," requires a dataFrame or spatialPolygonsDataFrame for the first argument or dF=\n")
      return(FALSE) 
    }   
     
  } else
  {
    stop("internal rworldmap error inputNeeded should be sPDF, sPDF or dF, dF, here it is:",inputNeeded,"in",functionName, "from", callingFunction)     
  }
  
  
  #later may add this here or put in another function
  ## check that the column name exists in the data frame
  #if ( is.na(match(nameColumnToPlot, names(mapToPlot@data)) )){
  #  stop("your chosen nameColumnToPlot :'",nameColumnToPlot,"' seems not to exist in your data, columns = ",paste(names(mapToPlot@data),""))
  #  return(FALSE)
  #} 
  
  
  invisible(inputData)
  
} # end of rwmCheckAndLoadInput 
