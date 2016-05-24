#' Map country-level data.
#' 
#' Draw a map of country-level data, allowing countries to be coloured, from an
#' object created in \code{\link{joinCountryData2Map}}.
#' 
#' Certain catMethod and colourPalette options go well together. e.g.
#' "diverging" and "diverging", "categorical" and "rainbow"
#' 
#' There are two styles of legend available.  If catMethod='categorical' or the
#' packages fields and spam are not installed a simple legend with coloured
#' boxes is created. Otherwise a colour bar legend is created. Finer control
#' can be achieved by \code{\link{addMapLegendBoxes}} or
#' \code{\link{addMapLegend}} repectively.
#' 
#' @param mapToPlot a spatial polygons dataframe from joinCountryData2Map()
#' containing country polygons and data, if none specified an internal example
#' data is used
#' @param nameColumnToPlot name of column containing the data you want to plot
#' @param numCats number of categories to put the data in, may be modified if
#' this number is incompatible with the catMethod chosen
#' @param xlim map extents c(west,east), can be overidden by mapRegion
#' @param ylim map extents c(south,north), can be overidden by mapRegion
#' @param mapRegion a country name from getMap()[['NAME']] or
#' 'world','africa','oceania','eurasia','uk' sets map extents, overrides
#' xlim,ylim
#' @param catMethod method for categorisation of data : \enumerate{
#' \item "categorical" - each unique value is treated as a separate category
#' \item for numeric data : "pretty", "fixedWidth", "diverging",
#' "logFixedWidth", "quantiles" \item a numeric vector defining breaks e.g.
#' c(0:5), note that a value of 2 goes into 1-2 not 2-3, uses
#' cut(include.lowest=TRUE) }
#' @param colourPalette string describing the colour palette to use, choice of:
#' \enumerate{ \item "palette" for the current palette \item a vector of valid
#' colours, e.g. =c('red','white','blue') or output from RColourBrewer \item one
#' of "heat", "diverging", "white2Black", "black2White", "topo", "rainbow",
#' "terrain", "negpos8", "negpos9" }
#' @param addLegend whether to add a legend or not
#' @param borderCol the colour for country borders
#' @param mapTitle title to add to the map, any string or 'columnName' to set
#' it to the name of the data column
#' @param oceanCol a colour for the ocean
#' @param aspect aspect for the map, defaults to 1, if set to 'variable' uses
#' same method as plot.Spatial in sp
#' @param missingCountryCol a colour for missing countries
#' @param add whether to add this map on top of an existing map, TRUE/FALSE
#' @param nameColumnToHatch allows hatching of country fills (e.g. to represent
#' uncertainty) , specify a column containing numeric data , highest values
#' will be solid and lower values will have a decreasing density of hatching ,
#' new feature more documentation will be added soon
#' @param lwd line width for country borders
#' @return invisibly returns a list containing the data and main options used
#' for the map, the list can be passed to \code{\link{addMapLegend}} or
#' \code{\link{addMapLegendBoxes}} along with additional options to allow
#' greater flexibility in legend creation.
#' @section Warning: will generate unhelpful errors in data categorisation if
#' inappropriate options are chosen, e.g. with catMethod:Quantiles if numCats
#' too high so that unique breaks cannot be defined.
#' @author andy south
#' @seealso classInt, RColorBrewer
#' @keywords aplot
#' @examples
#' 
#' mapCountryData()
#' data("countryExData",envir=environment(),package="rworldmap")
#' sPDF <- joinCountryData2Map(countryExData
#'               , joinCode = "ISO3"
#'               , nameJoinColumn = "ISO3V10"
#'               )
#' mapCountryData( sPDF
#'               , nameColumnToPlot="BIODIVERSITY" 
#'               )
#'               
#' #user defined map colour scheme for categorical data              
#' mapParams <- mapCountryData(nameColumnToPlot='GEO3major'
#'               , catMethod='categorical'
#'               , addLegend='FALSE'
#'               , colourPalette=c('white','green','red','yellow','blue','black') 
#'               )
#' #changing legendText
#' mapParams$legendText <- c('antarctic','africa','oceania'
#'                          ,'americas','s.asia','eurasia')              
#' do.call( addMapLegendBoxes, c(mapParams,x='bottom',title="Region",horiz=TRUE))
#' 
#' ##showing how rworldmap can be used with the classInt and RColorBrewer packages
#' library(classInt)
#' library(RColorBrewer)
#' #getting example data and joining to a map
#' data("countryExData",envir=environment(),package="rworldmap")
#' sPDF <- joinCountryData2Map(countryExData,joinCode = "ISO3"
#'                            ,nameJoinColumn = "ISO3V10")
#' #getting class intervals using a 'jenks' classification in classInt package
#' classInt <- classIntervals( sPDF$EPI, n=5, style="jenks")
#' catMethod = classInt$brks
#' #getting a colour scheme from the RColorBrewer package
#' colourPalette <- brewer.pal(5,'RdPu')
#' #calling mapCountryData with the parameters from classInt and RColorBrewer
#' mapParams <- mapCountryData( sPDF, nameColumnToPlot="EPI", addLegend=FALSE
#'                            , catMethod = catMethod, colourPalette=colourPalette )
#' do.call(addMapLegend, c(mapParams
#'                        ,legendLabels="all"
#'                        ,legendWidth=0.5
#'                        ,legendIntervals="data"))
#' 
#'               
#' 
#' 
#' @export mapCountryData
mapCountryData <- function(
                mapToPlot =         ""
              , nameColumnToPlot =  ""
              , numCats =           7
              , xlim =              NA
              , ylim =              NA
              , mapRegion =         "world"
              , catMethod =         "quantiles"
              , colourPalette =     "heat"
              , addLegend =         TRUE
              , borderCol =         'grey'
              , mapTitle =          'columnName'
              , oceanCol =          NA
              , aspect =            1
              , missingCountryCol = NA
              , add =               FALSE
              , nameColumnToHatch = ""    
              , lwd =               0.5  
              ){
                           
  functionName <- as.character(sys.call()[[1]])
                           

  #28/6/2013 refactoring
  new <- TRUE
  if (new)
  {
    #mapToPlot <- rwmCheckAndLoadInput( mapToPlot, requireSPDF = TRUE, callingFunction=functionName )    
    mapToPlot <- rwmCheckAndLoadInput( mapToPlot, inputNeeded="sPDF", callingFunction=functionName )  
  } else
  {
    if ( class(mapToPlot)=="SpatialPolygonsDataFrame" ) {
      ## checking if there is any data in the dataFrame
      if ( length(mapToPlot@data[,1]) < 1 ){
        stop("seems to be no data in your chosen file or dataframe in ",functionName) 
        return(FALSE)
      } 
    } else if ( mapToPlot == "" ) {
      message(paste("using example data because no file specified in",functionName))
      mapToPlot <- getMap(resolution="coarse")
      
      ## also setting a default nameColumnToPlot if it isn't set
      if ( nameColumnToPlot == "" ) nameColumnToPlot <- "POP_EST" #
    } else {
      stop(functionName," requires a SpatialPolygonsDataFrame object created by the joinCountryData2Map() function \n")
      return(FALSE) 
    } 
       
  } #end of replaced bit 28/6/2013
  
  ## setting a default nameColumnToPlot if it isn't set
  if ( nameColumnToPlot == "" ) nameColumnToPlot <- "POP_EST" #  
  
  ## check that the column name exists in the data frame
  if ( is.na(match(nameColumnToPlot, names(mapToPlot@data)) )){
    stop("your chosen nameColumnToPlot :'",nameColumnToPlot,"' seems not to exist in your data, columns = ",paste(names(mapToPlot@data),""))
    return(FALSE)
  } 
  
  ##classify data into categories   
  dataCategorised <- mapToPlot@data[[nameColumnToPlot]]

  #if data are not numeric then set catMethod to categorical
  if ( ! is.numeric(dataCategorised) && catMethod != "categorical" )
     {
     catMethod = "categorical"
     message(paste("using catMethod='categorical' for non numeric data in",functionName))
     }
  
  #checking whether method is categorical, length(catMethod)==1 needed to avoid warning if a vector of breaks is passed  
  if( length(catMethod)==1 && catMethod=="categorical" ) #if categorical, just copy the data, add an as.factor() to convert any data that aren't yet as a factor   
    { 
      dataCategorised <- as.factor( dataCategorised )
      cutVector <- levels(dataCategorised) #doesn't do cutting but is passed for use in legend
      numColours <- length(levels(dataCategorised))
    
    }else if( is.numeric(catMethod) ) 	 
    {	
      #if catMethod is numeric it is already a vector of breaks
      cutVector <- catMethod
      #set numColours from the passed breaks
      numColours <- -1 + length(catMethod)
      #Categorising the data, using a vector of breaks.	
      dataCategorised <- cut( dataCategorised, cutVector, include.lowest=TRUE)          
      
  	} else if( is.character(catMethod) )
  	{
  	  cutVector <- rwmGetClassBreaks( dataCategorised, catMethod=catMethod, numCats=numCats, verbose=TRUE )
  	  #Categorising the data, using a vector of breaks.	
  	  dataCategorised <- cut( dataCategorised, cutVector, include.lowest=TRUE)   
  	  #set numColours from the classified data
      #numColours <- length(dataCategorised) #!*! BUG 7/11/12
  	  numColours <- length(levels(dataCategorised))
  	}
      
  
  ## add extra column to map attribute data
  colNameRaw <- nameColumnToPlot
  colNameCat <- paste(colNameRaw,"categorised",sep='')    
  mapToPlot@data[[colNameCat]] <- dataCategorised     
  
  ## how many colours : numCats may be overriden (e.g. for 'pretty') 	
  ## get vector of the colours to be used in map (length=num categories)    
  colourVector <- rwmGetColours(colourPalette,numColours)
  
  ## get numeric index of which category each datapoint is in (length = num points)  
  dataCatNums <- as.numeric(dataCategorised)
  
  #adding missing country colour
  if(!is.na(missingCountryCol)){
    #adding missing country colour as the last element
    colourVector<- c(colourVector,missingCountryCol)
    #setting all missing values to the last element
    dataCatNums[is.na(dataCatNums)]<-length(colourVector)
  }

  #Scale hatching variable (and invert). Then threshold above a certain value to secure solid status
  hatchVar = NULL
  if (nameColumnToHatch=='')
     {
      #setting up the map plot
      if (!add) rwmNewMapPlot(mapToPlot,mapRegion=mapRegion,xlim=xlim,ylim=ylim,oceanCol=oceanCol,aspect=aspect)
      #plotting the map
      plot(mapToPlot, col=colourVector[dataCatNums], border=borderCol, add=TRUE, usePolypath=FALSE, lwd=lwd)#29/9/2012
      } else  
    {
    #*HATCHING OPTION*
     
      hatchVar = mapToPlot@data[[nameColumnToHatch]]
      
      hatchVar = (hatchVar - min(hatchVar, na.rm=TRUE))/max(hatchVar, na.rm=TRUE)
      hatchVar = 1-hatchVar
      hatchVar = (hatchVar*50) + 30
      hatchVar[hatchVar > 79] = -1
      #hatchVar = (hatchVar*70) + 40
      #hatchVar = (hatchVar*70) + (hatchVar^2)/1000

      #setting up the map plot
      if(!add)  rwmNewMapPlot(mapToPlot,mapRegion=mapRegion,xlim=xlim,ylim=ylim,oceanCol=oceanCol,aspect=aspect)
      #plotting the map
      plot(mapToPlot,col=colourVector[dataCatNums],border=borderCol, density=hatchVar, angle=135, lty=1,add=TRUE,usePolypath=FALSE, lwd=lwd)#29/9/2012
      plot(mapToPlot,col=colourVector[dataCatNums],border=borderCol, density=hatchVar, angle=45, lty=1,add=TRUE,usePolypath=FALSE, lwd=lwd)#29/9/2012
                 
     }  #end of hatching option



  if (addLegend){
      

			if((length(catMethod)==1 && catMethod=="categorical") ){
      
        # simpler legend for categorical data OR if you don't have packages spam or fields.
        addMapLegendBoxes(colourVector=colourVector,cutVector=cutVector,catMethod=catMethod) #,plottedData=dataCategorised)          
      
      }else{
        #colour bar legend based on fields package
        addMapLegend(cutVector=cutVector,colourVector=colourVector) #,catMethod=catMethod) # ,plottedData=mapToPlot@data[[nameColumnToPlot]],catMethod=catMethod,colourPalette=colourPalette)
      }  
  
  } #end of addLegend
  
  ## add title
  if ( mapTitle == 'columnName' ){
    title(nameColumnToPlot)
  } else {
    title( mapTitle )
  }
   
  ##returning parameter list that can be used by do.call(addMapLegend,*)  
  invisible(list(colourVector=colourVector
                ,cutVector=cutVector
                ,plottedData=mapToPlot[[nameColumnToPlot]]
                ,catMethod=catMethod
                ,colourPalette=colourPalette
                )
           ) 
  
  #failed attempt at creating something that could be directly used in addMapLegend()         
  #invisible(list(plottedData=paste("'",sys.call()[[2]],"'",sep='')
  #              ,nameColumnToPlot=paste("'",nameColumnToPlot,"'",sep='')
  #              ,catMethod=paste("'",catMethod,"'",sep='')
  #              ,colourPalette=paste("'",colourPalette,"'",sep='')
  #              ,numCats=numCats
  #              )
  #         )            
           
            
} #end of mapCountryData()

