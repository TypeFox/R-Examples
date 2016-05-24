#' Add a legend to a map
#' 
#' Creates a colour bar legend, showing the range of colours and the values the
#' colours correspond to. Relies heavily on image.plot() from the package
#' fields. For simple use, simply use addLegend=TRUE in a rworldmap map
#' function. Or users can call addMapLegend seperately to fine tune the legend.
#' The user should insure that data, catMethod,numCats and colourPalette match
#' the values used in the plot.  The legend is designed to be useful for the
#' variety of classification methods that exist.
#' 
#' The default legend is a horizontal colour bar, with labels only at the
#' extremes.
#' 
#' Can use a parameter list returned from mapping functions, e.g.
#' mapCountryData(). mapCountryData(addLegend=TRUE) produces same results as:
#' mapParams <- mapCountryData(addLegend=FALSE) do.call(addMapLegend,
#' mapParams)
#' 
#' Using the following allows the modification of the legend : mapParams <-
#' mapCountryData(addLegend=FALSE) do.call(addMapLegend, c(mapParams,
#' legendLabels="all", legendWidth=0.5))
#' 
#' @param cutVector the categories or breaks used in the map
#' @param colourVector colours used in the map
#' @param legendLabels Controls the style of the labels on the legend. Choose
#' "none" for no labels, "limits" for the two end values, and "all" to show all
#' the break values if they fit.
#' @param labelFontSize Controls font size of the labels. A multiplier, so use
#' 2 to double the size, 0.5 to halve it, etc.
#' @param legendWidth Controls the width of the colour bar.
#' @param legendShrink Controls the length of the colour bar. 1 means full
#' width of the plot.
#' @param legendMar Moves the legend away from the side of the plot. Measured
#' in character widths.
#' @param horizontal If TRUE the legend is horizontal, if FALSE, vertical.
#' @param legendArgs For producing titles and labels. A list of arguments to be
#' passed to mtext.
#' @param tcl Controls the length of the tick marks.Useful when labelFontSize
#' is changed.
#' @param mgp Numeric vector length 3. The second element controls the distance
#' between labels and the axis. Useful when labelFontSize is changed.
#' @param sigFigs The number of significant figures for legend labels.
#' @param digits An argument to the formatting of the labels
#' @param legendIntervals "page" or "data". Controls the division of the colour
#' bar, "page" sets the intervals equal on the page, "data" sets them to be
#' equal in the units of the data.
#' @param plottedData unused but are passed with mapParams
#' @param catMethod unused but are passed with mapParams
#' @param colourPalette unused but are passed with mapParams
#' 
#' @return Adds a legend to a plot.
#' @note Can have the unintentional effect of modifying graphical parameters,
#' e.g. mfcol reverts to mfrow.
#' @author Andy South
#' @importFrom fields image.plot
#' @seealso mapCountryData, mapGriddedData, image.plot
#' @keywords aplot
#' @examples
#' 
#' #Set up the plot so the world map uses the full width.
#' mapDevice() 
#' #join eaxmple data to a map  
#' data("countryExData",envir=environment())
#' sPDF <- joinCountryData2Map(countryExData
#'               , joinCode = "ISO3"
#'               , nameJoinColumn = "ISO3V10"
#'               )
#' #map the data with no legend              
#' mapParams <- mapCountryData( sPDF
#'               , nameColumnToPlot="BIODIVERSITY"
#'               , addLegend='FALSE' 
#'               )
#'               
#' #add a modified legend using the same initial parameters as mapCountryData               
#' do.call( addMapLegend, c( mapParams
#'                         , legendLabels="all"
#'                         , legendWidth=0.5
#'                         ))
#' 
#' 
#' 
#' 
#' @export addMapLegend
addMapLegend <- function(            
                     colourVector=""
                    ,cutVector=""
                    
                    ,legendLabels="limits"              #style of bar labels 'all','none','limits'
                    ,labelFontSize=1                    #cex.axis by another name
                    ,legendWidth=1.2                    #Width in characters of the legend
                    ,legendShrink=0.9                   #Shrinks the legend, lengthwise
                    ,legendMar=3                       #shifts the legend upwards when horiz, measured in characters.
                    ,horizontal=TRUE                    #orientation
                    ,legendArgs=NULL                   #allows a title above legend
                    ,tcl=-.5                            #as per par(tcl=) tick length par default = -.5
                    ,mgp=c(3,1,0)                       #as per par(mgp=) margin position par default= c(3,1,0)
                    ,sigFigs=4                    #controls how numbers get rounded
                    ,digits=3                           #controls how numbers get formatted into neater numbers.
                    ,legendIntervals='data'       #page" or "data"."page"=intervals equal on page, "data"= equal in data units
                    
                    ,plottedData=""               #not used yet but maybe in future
                    ,catMethod="pretty"           #not used yet but maybe in future
                    ,colourPalette="heat"         #not used yet but maybe in future
                    #,missingCountryCol="white"    #not used yet but maybe in future                    
                                        
                    ){

#BEWARE image.plot from fields package at end modifies the par settings
#seemingly not possible to stop this, e.g. can't query whether mfrow or mfcol
#oldPar <- par(no.readonly = TRUE)


#this checks that length of colour vector is one less than length of cutVector
#if it isn't could be because a missingCountryCol has been added by mapCountryData
#for now remove the last colour, later may want to try to deal with missingCountryCol
if ( length(colourVector)  == length(cutVector) ) colourVector <- colourVector[-length(colourVector)]


#Simplify the plotBreaks. By rounding the numbers, it becomes easier to read.
colourBarBreaks = as.numeric(cutVector)
tidyPlotBreaks <- signif(colourBarBreaks,sigFigs) #

#The image.plot zlim argument only requires the min and max
zlim <- range(colourBarBreaks,na.rm=TRUE)

# adding in equal scale intervals options
if ( legendIntervals == 'page' )
   {
    colourBarBreaks <- colourBarBreaks[1] + (colourBarBreaks[length(colourBarBreaks)]-colourBarBreaks[1])* seq(from=0,to=1,length.out=length(colourBarBreaks))
   }

#axis.args at = positioning, labels = text displayed.
if(legendLabels=="limits"){
  limitsIndex=c(1,length(colourBarBreaks))
  axis.args=list(at=colourBarBreaks[limitsIndex],cex.axis=labelFontSize,mgp=mgp,tcl=tcl,labels=prettyNum(tidyPlotBreaks[limitsIndex],digits=digits,format="G"))

} else if(legendLabels=="none"){
  #to get no labels need to change ac to whether horizontal
  if ( horizontal ) axis.args=list(xaxt="n")
  else axis.args=list(yaxt="n")
  
} else if(legendLabels=="all"){
  axis.args=list(at=colourBarBreaks,cex.axis=labelFontSize,mgp=mgp,tcl=tcl,labels=prettyNum(tidyPlotBreaks,digits=digits,format="G"))
}

#The actual legend plotting command
#**BEWARE this can cause mfcol to be set back to mfrow and seems to be no way to fix
image.plot(zlim=zlim,legend.only=TRUE,horizontal=horizontal,legend.args=legendArgs,legend.mar=legendMar,col=colourVector,breaks=colourBarBreaks,axis.args=axis.args,legend.width=legendWidth,legend.shrink=legendShrink)

#image.plot(zlim=zlim,legend.only=TRUE,graphics.reset=TRUE,horizontal=horizontal,legend.args=legendArgs,legend.mar=legendMar,col=coloursUsed,breaks=colourBarBreaks,axis.args=axis.args,legend.width=legendWidth,legend.shrink=legendShrink)

#par(oldPar)
#print('test')

}


