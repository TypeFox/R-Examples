#' Produce maps of global gridded data at half degree resolution
#' 
#' Produce maps of global gridded data at half degree resolution
#' 
#' Plots a map of global half degree gridded data, allowing classification,
#' colours and regions to be set.
#' 
#' Certain catMethod and colourPalette options go well together. e.g.
#' "diverging" and "diverging", "categorical" and "rainbow"
#' 
#' @param dataset gridded data either as a : \enumerate{ \item
#' SpatialGridDataFrame (R object defined in package sp) \item file name of a
#' GridAscii file - this is an Esri format \item 2D R matrix or array (rows by
#' columns) }
#' @param nameColumnToPlot name of column containing the data to plot
#' @param numCats number of categories to put the data in, may be overidden if
#' catMethod ='pretty'
#' @param catMethod method for categorisation of data "pretty",
#' "fixedWidth","diverging", "logFixedWidth","quantiles","categorical", or a
#' numeric vector defining breaks
#' @param colourPalette a string describing the colour palette to use, choice
#' of : \enumerate{ \item "palette" for the current palette \item a vector of
#' valid colours, e.g. =c('red','white','blue') or output from RColourBrewer
#' \item one of "heat", "diverging", "white2Black", "black2White", "topo",
#' "rainbow", "terrain", "negpos8", "negpos9" }
#' @param xlim map extents c(west,east), can be overidden by mapRegion
#' @param ylim map extents c(south,north), can be overidden by mapRegion
#' @param mapRegion a country name from getMap()[['NAME']] or
#' 'world','africa','oceania','eurasia','uk' sets map extents, overrides
#' xlim,ylim
#' @param addLegend whether to add a legend or not
#' @param addBorders options for country borders, 'low','coarse' = low or
#' coarse resolution, 'coasts' = coasts only, 'none' or NA for none
#' @param borderCol the colour for country borders
#' @param oceanCol a colour for the ocean if the grid values are NA
#' @param landCol a colour to fill countries if the grid values are NA over
#' land
#' @param plotData whether to plotData, if FALSE a legend can be added on its
#' own
#' @param aspect aspect for the map, defaults to 1, if set to 'variable' uses
#' same method as plot.Spatial in sp
#' @param lwd line width for country borders
#' @return invisibly returns a list containing the data and main options used
#' for the map, the list can be passed to \code{\link{addMapLegend}} along with
#' additional options to allow greater flexibility in legend creation.
#' @author andy south
#' @importFrom maptools readAsciiGrid
#' @seealso classInt, RColorBrewer
#' @keywords hplot
#' @examples
#' 
#' ## mapping continuous data
#' data(gridExData,envir=environment(),package="rworldmap")
#' gridExData <- get("gridExData")
#' mapGriddedData(gridExData)
#' 
#' ## reclassing continuous data to categorical & mapping
#' data(gridExData,envir=environment(),package="rworldmap")
#' #find quartile breaks
#' cutVector <- quantile(gridExData@@data[,1],na.rm=TRUE)
#' #classify the data to a factor
#' gridExData@@data$categories <- cut( gridExData@@data[,1]
#'                                       , cutVector, include.lowest=TRUE)
#' #rename the categories
#' levels(gridExData@@data$categories) <- c('low', 'med', 'high', 'vhigh')
#' #mapping
#' mapGriddedData( gridExData, nameColumnToPlot= 'categories'
#'               , catMethod='categorical')
#' 
#' 
#' 
#' @export mapGriddedData
mapGriddedData <- function(
                           dataset =          ""
                         , nameColumnToPlot = "" 
                         , numCats =          5  
                         , catMethod =        "quantiles"   
                         , colourPalette =    "heat"  
                         , xlim =             c(-180,180)
                         , ylim =             c(-80,90) 
                         , mapRegion =        "world"   
                         , addLegend =        TRUE
                         , addBorders =       'coarse' 
                         , borderCol =        'grey'
                         , oceanCol =         NA
                         , landCol =          NA
                         , plotData =         TRUE
                         , aspect =           1
                         , lwd =              1  
                         )
   {

    #browser()
    #print("test")

    functionName <- as.character(sys.call()[[1]])

    
    ## filename or nothing ##
    if (class(dataset)=='character')
       {
        if (dataset=="") #if no dataset passed
           {
            #open example one
            data(gridExData,envir=environment(),package="rworldmap")
            sGDF <- get("gridExData")
           } else
            sGDF <- readAsciiGrid(dataset)
    ## matrix or array ##        
       } else if (class(dataset)=='matrix' || class(dataset)=='array')   
       {
        if ( length(dim(dataset)) == 2 )
           {
            #assume that it's just 2 dimensions
            gridVals <- data.frame(att=as.vector(dataset))
           
            #this assumes that the matrix covers the whole globe -180 to 180, -90 to 90
            gt <- GridTopology( cellcentre.offset = c(-180,-90)
                      , cellsize = c( 360/dim(dataset)[1], 180/dim(dataset)[2] )
                      , cells.dim = dim(dataset)[1:2] )
    
            sGDF <- SpatialGridDataFrame(gt, data = gridVals)
           } else
           {
            stop("the first argument to ",functionName," if a matrix or array should have 2 dimensions, yours has, ", length(dim(dataset))) 
            return(FALSE)
           }           
    ## SGDF passed ##           
       } else if (class(dataset)=='SpatialGridDataFrame')
       {   
        sGDF <- dataset
       } else
       {
    ## !! I could add option here for dataFrame with nameX & nameY columns   
        stop("the first argument to ",functionName," should be a file name, 2D array or matrix, or SpatialGridDataFrame, yours is, ", class(dataset)) 
        return(FALSE)
       } 


     #if the sGDF contains multiple attribute columns, decide which to plot
     if ( length(sGDF@data) == 1 ) attrName <- names(sGDF)[1] #to be able to get at data using original filename
     else if ( length(sGDF@data) > 1 && nameColumnToPlot != "" )
     {
        attrName <- nameColumnToPlot
        #!then want to check wkether this column is an attribute column in the sGDF
        
     } else if ( length(sGDF@data) > 1 && nameColumnToPlot == "" )
     {
         attrName <- names(sGDF)[1]
         message("plotting the first data column because nameColumnToPlot not specified in mapGridAscii()\n")
     }

    #CLASSIFYING THE DATA
    
    #first get the raw data    

    dataCategorised <- sGDF[[attrName]]
    
    #print(dataCategorised[1:10])  
    
    #browser()
    
    ####################### categorical 
    #length(catMethod)==1 needed to avoid warning if a vector of breaks is passed  
    if( length(catMethod)==1 && catMethod=="categorical" )    
      {
       #if categorical, just copy the data, add an as.factor() to convert any data that aren't yet as a factor 
       dataCategorised <- as.factor( dataCategorised )
       cutVector <- levels(dataCategorised) #doesn't do cutting but is passed for use in legend
  		 #6/5/10 bug fix
       #numColours <- -1 + length(levels(dataCategorised))
       #12/10/10 the above seemed to cause one less colour than cat
       numColours <- length(levels(dataCategorised))
       
       #1/11/10 STILL A PROBLEM WITH THIS 
       #SEEMS THAT NOW WITH ELISABETHS EXAMPLE IT GIVES SAME NUM BREAKS AS COLOURS WHICH
       #GENERATES AN ERROR ???
       #******temp fix********
       #numColours <- -1 + length(levels(dataCategorised))
       
       #12/10/10 added in here to avoid problem when doing the same as non-categorical
       #this deliberately gets the numeric index of each factor
       sGDF$indexToPlot <-  as.numeric(dataCategorised) 
       
       #13/11/10 fixing breaks difference between categorical and non-categorical
       breaks <- c(0:(length(cutVector)) )
       
      }else
      ####################### all other catMethods except categorical 
      {
       #if catMethod is not a vector of numbers 
       if(is.character(catMethod)==TRUE)
      	{	
      		cutVector <- rwmGetClassBreaks( dataCategorised, catMethod=catMethod, numCats=numCats, verbose=TRUE )
      		#6/5/10 bug fix
      		#28/7/10 removed to avoid bug with default params
          #numColours <- -1 + length(levels(dataCategorised))
          numColours <- -1 + length(cutVector)
        } else if(is.numeric(catMethod)==TRUE)
      	#if catMethod is numeric it is already a vector of breaks	
      	{
      		cutVector <- catMethod
      		#6/5/10 bug fix
      		numColours <- -1 + length(cutVector)
      	}
    	#Categorising the data, using a vector of breaks.	     
      dataCategorised <- cut( dataCategorised, cutVector, include.lowest=TRUE,labels=FALSE)
      
      #12/10/10 moved into non-categorical loop from below
      sGDF$indexToPlot <- as.numeric( as.character( dataCategorised )) 
      #cut : labels for the levels of the resulting category.  By default,
      #    labels are constructed using ?"(a,b]"? interval notation.
      #    If ?labels = FALSE?, simple integer codes are returned instead of a factor.
      
      #13/11/10 fixing breaks difference between categorical and non-categorical
      breaks <- c(0:(length(cutVector)-1) )
        	
  	  } #end of if data are not categorical
  

    
    colourVector <- rwmGetColours(colourPalette,numColours)

    #setting up the map plot 
    #fills in the ocean but will get overpainted by the grid data if it doesn't have NAs in the ocean
    rwmNewMapPlot(mapToPlot=sGDF,oceanCol=oceanCol,mapRegion=mapRegion,aspect=aspect,xlim=xlim,ylim=ylim)


    #to fill in any countries with NA values in the grid
    if(!is.na(landCol))
       {
        plot( getMap(), add=TRUE, border=borderCol, col=landCol, lwd=lwd )
       }

    #only plot ascii data if plotData=T (allows legend to be plotted on its own by setting plotData=F)
    if (plotData)
       {
        #image(sGDF,add=TRUE,attr='indexToPlot',col=colourVector, xaxs='i', yaxs='i' ) #xaxs=i ensures maps fill plot area
        #7/5/10 !BUG without a breaks param then image fits colours to just those cats present
        #image(sGDF,add=TRUE,attr='indexToPlot',col=colourVector, xaxs='i', yaxs='i',breaks=c(0:(length(cutVector)-1) )) #xaxs=i ensures maps fill plot area
        #31/10/10 setting add=FALSE otherwise breaks doesn't work
        #browser()
        #image(sGDF,add=FALSE,attr='indexToPlot',col=colourVector, xaxs='i', yaxs='i',breaks=c(0:(length(cutVector)-1) ))
        #breaks still didn't work
        #the below does by seprating out the sp bit from the graphics bit : perhaps is a problem in sp ?
        xyz <- as.image.SpatialGridDataFrame(sGDF,attr='indexToPlot')
        #image(xyz,add=TRUE,col=colourVector, xaxs='i', yaxs='i',breaks=c(0:(length(cutVector)-1) ))
        #13/11/10 fixing breaks difference between categorical and non-categorical
        image(xyz,add=TRUE,col=colourVector, xaxs='i', yaxs='i',breaks=breaks)
       }
       
       
    borderOptions = c('low','coarse','coasts',NA,'','none')
    if (addBorders=='low'){
       plot( getMap(resolution='low'), add=TRUE, border=borderCol, lwd=lwd )
       } else
    if (addBorders=='coarse'){
       plot( getMap(resolution='coarse'), add=TRUE, border=borderCol, lwd=lwd )
       } else
    if (addBorders=='coasts'){
       data(coastsCoarse,envir=environment(),package="rworldmap")
       coastsCoarse <- get("coastsCoarse")
       plot(coastsCoarse, add=TRUE, col=borderCol, lwd=lwd ) 
       } else 
    if ( ! addBorders %in% borderOptions){
       warning("unrecognised addBorders = ",addBorders, "none plotted, choose one of",paste(borderOptions,""))
       }             
    
    ## adding a default legend, can be modified by calling addMapLegend() independently  
    if (addLegend){
    
      #if((length(catMethod)==1 && catMethod=="categorical") || !require("spam") || !require("fields")){
      #20/8/13 removed require bits
			if((length(catMethod)==1 && catMethod=="categorical") ){
        
        #legend(x='bottomleft', legend=c(rev(levels(dataCategorised)),"no data"), pch = 22, pt.cex=2, col=borderCol,pt.bg=c(coloursForMap[numColours:1],"white"), title="category",bg="white" )
        addMapLegendBoxes(colourVector=colourVector,cutVector=cutVector,plottedData=dataCategorised,catMethod=catMethod)          
         
      }else{
        #colour bar legend based on fields package
        addMapLegend(colourVector=colourVector,cutVector=cutVector) #,plottedData=sGDF[[attrName]],catMethod=catMethod,colourPalette=colourPalette)   
        }
      }

    #could add title
    #!but need to set it to the filename for gridascii files,
    #!and the column name for multi-attribute sGDFs
    #if ( mapTitle == 'columnName' ) title(nameColumnToPlot)
    #else title( mapTitle )

    #returning data to be used by addMapLegend
    invisible(list(plottedData=sGDF[[attrName]]
                  ,catMethod=catMethod
                  ,colourVector=colourVector
                  ,cutVector=cutVector
                  ,colourPalette=colourPalette
                  )
             )              

    } #end of mapGriddedData()





















