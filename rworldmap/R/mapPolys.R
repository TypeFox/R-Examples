#' Map polygon data.
#' 
#' Plot a map of polygons, from a spatialPolygonsDataFrame, coloured according
#' to one a specified attribute column.
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
#' @param mapToPlot a spatial polygons dataframe (e.g. from joinData2Map())
#' containing polygons and associated data, if none specified an internal
#' example data is used
#' @param nameColumnToPlot name of column containing the data you want to plot
#' @param numCats number of categories to put the data in, may be modified if
#' this number is incompatible with the catMethod chosen
#' @param xlim map extents c(west,east), can be overidden by mapRegion
#' @param ylim map extents c(south,north), can be overidden by mapRegion
#' @param mapRegion a country name from getMap()[['NAME']] or
#' 'world','africa','oceania','eurasia','uk' sets map extents, overrides
#' xlim,ylim
#' @param catMethod for categorisation of data "pretty", "fixedWidth",
#' "diverging", "logFixedWidth", "quantiles", "categorical", or a numeric
#' vector defining breaks
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
#' @param lwd line width for country borders
#' @return invisibly returns a list containing the data and main options used
#' for the map, the list can be passed to \code{\link{addMapLegend}} or
#' \code{\link{addMapLegendBoxes}} along with additional options to allow
#' greater flexibility in legend creation.
#' @author andy south
#' @seealso joinData2Map, classInt, RColorBrewer
#' @keywords aplot
#' @examples
#' 
#' 
#' ## this example uses downloaded files
#' ## to run it download the files
#' ## and remove the comment symbols '#' from all the lines starting with a single '#'
#' 
#' ## US states map downloaded from :
#' ## http://www2.census.gov/cgi-bin/shapefiles2009/national-files
#' 
#' #inFile <- 'tl_2009_us_stateec.shp'
#' #sPDF <- readShapePoly(inFile)
#' #str(sPDF@@data)
#' 
#' ##################
#' ## use mapPolys to map the sPDF
#' #mapPolys(sPDF,nameColumnToPlot = "ALANDEC")
#' #mapPolys(sPDF,nameColumnToPlot = "AWATEREC",mapRegion='North America')
#' 
#' ##################
#' ## join some other data to it
#' ## education data downloaded from here as xls then saved as csv
#' ## http://nces.ed.gov/ccd/drpcompstatelvl.asp
#' 
#' #dataFile <- 'SDR071A_xls.csv'
#' #dF <- read.csv(dataFile,as.is=TRUE)
#' #str(dF)
#' ## STATENAME
#' ## DRP912 Dropout Rate, Grades 9 through 12
#' 
#' ## joining the data to the map
#' ## based upon state names (column NAMEEC in map, and STATENAME in the data)
#' #sPDF2 <- joinData2Map(dF
#' #        , nameMap = sPDF
#' #        , nameJoinIDMap = "NAMEEC"
#' #        , nameJoinColumnData = "STATENAME")
#' 
#' #################
#' ## plot one of the attribute variables
#' #mapDevice()# to set nice shape map window
#' #mapPolys(sPDF2,nameColumnToPlot = "DRP912",mapRegion='North America')
#' 
#' 
#' #################
#' ###to map US counties data (Tiger) downloaded from :
#' ##http://www2.census.gov/cgi-bin/shapefiles2009/national-files
#' 
#' #inFile <- 'tl_2009_us_county.shp'
#' #sPDF <- readShapePoly(inFile)
#' #str(sPDF@@data)
#' #mapPolys(sPDF,nameColumnToPlot='AWATER',xlim=c(-140,-65), ylim=c(25,45))
#' 
#' 
#' 
#' @export mapPolys
mapPolys <- function(
                           mapToPlot = "",
                           nameColumnToPlot = "",
                           
                           numCats = 7, # *may be overridden by catMethod
                           xlim=NA, 
                           ylim=NA, 
                           mapRegion = "world",   #sets map extents, overrides xlim, ylim
                           catMethod="quantiles",   #any vector defining breaks, "fixedWidth","quantiles","logFixedWidth"
                           colourPalette= "heat", #"heat","white2Black","topo","palette" for current palette
                           addLegend=TRUE,
                           borderCol = 'grey',
                           mapTitle = 'columnName', #this sets it to the name of the column, any other string can be passed too
                           oceanCol=NA,
                           aspect=1,
                           missingCountryCol=NA,
                           add=FALSE, 
                           lwd=0.5
                           ){
                           
  functionName <- as.character(sys.call()[[1]])
                           
  #browser()  
  
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
        ## also setting a defsult nameColumnToPlot if it isn't set
        #if ( nameColumnToPlot == "" ) nameColumnToPlot <- "POP_EST" #
      } else {
        stop(functionName," requires a SpatialPolygonsDataFrame object created by the joinCountryData2Map() function \n")
        return(FALSE) 
      }

  } #end of replaced bit 28/6/2013
  
  ## setting a default nameColumnToPlot if it isn't set
  # moved out of above loop replaced by rwmCheckAndLoadInput
  if ( nameColumnToPlot == "" ) nameColumnToPlot <- "POP_EST" #  
 
  ## check that the column name exists in the data frame
  if ( is.na(match(nameColumnToPlot, names(mapToPlot@data)) )){
    stop("your chosen nameColumnToPlot :'",nameColumnToPlot,"' seems not to exist in your data, columns = ",paste(names(mapToPlot@data),""))
    return(FALSE)
  } 
    

  
  dataCategorised <- mapToPlot@data[[nameColumnToPlot]]

  #1/10/12 if the data are not numerical then set catMethod to categorical
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
    }else
    { 
      if(is.character(catMethod)==TRUE)
    	{	
    		cutVector <- rwmGetClassBreaks( dataCategorised, catMethod=catMethod, numCats=numCats, verbose=TRUE )
    	} else if(is.numeric(catMethod)==TRUE)
    	#if catMethod is numeric it is already a vector of breaks	
    	{
    		cutVector <- catMethod
    	}
  	#Categorising the data, using a vector of breaks.	
  	dataCategorised <- cut( dataCategorised, cutVector, include.lowest=TRUE)    	
	  } #end of if data are not categorical
 
  
  ## add extra column to map attribute data
  colNameRaw <- nameColumnToPlot
  colNameCat <- paste(colNameRaw,"categorised",sep='')    
  mapToPlot@data[[colNameCat]] <- dataCategorised     
  
  ## how many colours : numCats may be overriden (e.g. for 'pretty') 	
  numColours <- length(levels(dataCategorised))
  
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
  
  #setting up the map plot
  if (!add) rwmNewMapPlot(mapToPlot,mapRegion=mapRegion,xlim=xlim,ylim=ylim,oceanCol=oceanCol,aspect=aspect)
  #plotting the map
  plot(mapToPlot,col=colourVector[dataCatNums],border=borderCol,add=TRUE,lwd=lwd) #,density=c(20:200))#angle=c(1:360),)
  
  #trying out shading with density & angle
  #plot(mapToPlot,col=colourVector[dataCatNums],border=borderCol,add=TRUE,density=c(20:200))#angle=c(1:360),)
  #plot(mapToPlot,col='white',border=borderCol,add=TRUE,density=c(20:200),bg=colourVector[dataCatNums])#angle=c(1:360),)  
  #the above might need : xaxs="i",yaxs="i") #xaxs="i" ensures maps fill plot area
   
  if (addLegend){
      
      #if((length(catMethod)==1 && catMethod=="categorical") || !require("spam") || !require("fields")){
      #20/8/13 removed require bits
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
   
  ##29/10/09 returning parameter list that can be used by do.call(addMapLegend,*)  
  #sys.call()[[2]] gets the name of the first argument
        
  #invisible(list(plottedData=eval( parse(text=paste(sys.call()[[2]],"[['",nameColumnToPlot,"']]",sep='')))
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
           
            
} #end of mapPolys()

