#' Add a legend of coloured boxes to a map
#' 
#' Creates a colour box legend, showing the range of colours and the values the
#' colours correspond to. This works well for categorical data with relatively
#' few categories.
#' 
#' This creates a legend with separate boxes of colour rather than
#' addMapLegend() which creates a colour bar. This method is used as the
#' default for categorical data.
#' 
#' See the examples for how to use a parameter list returned from mapping
#' functions.
#' 
#' @param cutVector the categories or breaks used in the map
#' @param colourVector colours used in the map
#' @param x positioning of legend e.g. 'bottomleft', 'topright'
#' @param horiz if TRUE horizontal legend
#' @param title title for Legend
#' @param cex controls the font size, default is 1
#' @param pt.cex controls size of colour boxes relative to cex, default is 2
#' @param col colour for boundary of colour boxes, default is "gray"
#' @param bg colour for legend background, default is "white", NA makes the
#' legend background transparent
#' @param legendText the text to put against each legend box, if left blank
#' cutVector is used, needs to be a vector the same length as length cutVector
#' @param catMethod the categorisation method used influences what text added
#' to legend elements, for 'categorical' just the category names are used for
#' other options limits are used
#' @param plottedData not used yet but maybe in future
#' @param colourPalette not used yet but maybe in future
#' @param sigFigs not used yet but maybe in future
#' @param missingCountryCol not used yet but maybe in future
#' @param \dots to allow other params to be set in legend
#' @return Adds a legend to a plot.
#' @author Andy South
#' @seealso addMapLegend, mapCountryData, mapGriddedData
#' @keywords aplot
#' @examples
#' 
#' #Set up the plot so the world map uses the full width.
#' mapDevice() 
#' #map example categorical data with no legend              
#' mapParams <- mapCountryData(nameColumnToPlot='GEO3major'
#'               , catMethod='categorical'
#'               , addLegend='FALSE' 
#'               )
#'               
#' #add default legend using the same parameters as mapCountryData               
#' do.call( addMapLegendBoxes, c( mapParams))
#' 
#' #adding a modified legend by specifying extra parameters               
#' do.call( addMapLegendBoxes, c(mapParams,x='bottom',horiz=TRUE,title="Region"))
#' 
#' #user defined map colour sceme              
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
#' #or this way
#' #do.call( addMapLegendBoxes
#' #        , c(mapParams
#' #           ,list(legendText=c('antarctic','africa','oceania'
#' #                              ,'americas','s.asia','eurasia')
#' #                ,x='bottom',title="Region",horiz=TRUE)))
#' 
#' 
#' 
#' @export addMapLegendBoxes
`addMapLegendBoxes`<- function(
                    cutVector=""    # the categories or breaks used in the map
                    ,colourVector = "" #colours used in the map
                    ,x='bottomleft'
                    ,horiz=FALSE
                    ,title="category"
                    ,cex=1 #cex controls font size
                    ,pt.cex=2 #pt.cex controls size of colour boxes
                    ,col="gray" #boundary of boxes
                    ,bg="white" #legend background
                    ,legendText="" #if this is left as empty then the cut vector is used
                    
                    ,catMethod="categorical"           #
                    
                    ,plottedData=""               #not used yet but maybe in future
                    ,colourPalette="heat"         #not used yet but maybe in future
                    ,sigFigs=2                    #not used yet but maybe in future
                    ,missingCountryCol="white"    #not used yet but maybe in future

                    ,... #to allow other params to be set in legend
                    ){
                    
#function for categorical legend or if user doesn't have fields package

#!? deal with what happens if non categorical data get through
if ( catMethod!="categorical" )
   {
 	  #to set cutVector to 1-2,2-3,3-4 etc. for legend 
    #this is an ugly way of doing but it does work    
    func <- function(x,y) c(paste(x,"-",y[1+which(y==x)],sep=""))
    tmp <- sapply(cutVector,cutVector,FUN=func)
    cutVector <- tmp[1:length(tmp)-1] #removing last element from vector
   }

if (length(legendText)==1 && legendText=="") legendText=cutVector

legend( x=x
      , horiz=horiz
      , legend=legendText #cutVector #c(levels(plottedData),"no data")
      , pch = 22
      , cex=cex
      , pt.cex=pt.cex
      , pt.bg=c(colourVector)#,missingCountryCol) #or is missingCountryCol already added on ?
      , title=title
      , col=col
      , bg=bg
      , ...)


} #end of addMapLegendBoxes
