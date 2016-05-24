#' dcscatter adds a scatter plot to the dashboard
#'
#' \code{dcscatter} generates a scatter plot
#' 
#' @docType methods
#' @param x column name of a single numeric in data frame \code{data} for drawing a scatter plot
#' @param y column name of a single numeric in data frame \code{data} for drawing a scatter plot
#' @param title character for the title of the generated scatter plot
#' @param spansize integer between 1 to 12 for the width of the element in the web page
#' @param width integer for the width (in pixels) of the element in the web page
#' @param height integer for the height (in pixels) of the element in the web page
#' @param symbolesize integer for adjusting the symbole size
#' @param symboletype character for defining the symbole type
#' @param clipaddingsize integer for adjusting the clipadding size
#' @param hlightedsize integer for adjusting the highlighted size
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5])
#' dcscatter(x=names(iris)[1], y=names(iris)[3] )
#' dcbarchart(x=names(iris)[1] , gap=75)
#' dcpiechart(x=names(iris)[2])
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files and server is not launched
#'  

dcscatter <- function(
x,
y,
title=paste(x," * ",y),
spansize=4,
width=dashboard.env$outerwidth*spansize/12,
height=250, 
symbolesize=2,
symboletype="circle", #diamond
clipaddingsize =10,
hlightedsize=4
){
	# dimension name is used by convention (DC.js)
	#dimension = x

	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
		# close the previous row, write in a row then start a new one 
		linebreak()
	}
	
	dashboard.env$numeltgraph <- dashboard.env$numeltgraph + 1
	dashboard.env$nbelementinrow <- dashboard.env$nbelementinrow + 1	
	dashboard.env$sumspaninrow <- dashboard.env$sumspaninrow + spansize
		
	dashboard.env$temp <- paste(dashboard.env$temp,"<div class='span",spansize,"' id='dc-chart-",dashboard.env$numeltgraph,"'><h4>",title," <span><a class='reset' href='javascript:dcChart",dashboard.env$numeltgraph,".filterAll();dc.redrawAll();' style='display: none;'> reset </a> </span> </h4> </div>", sep="")	
	
			
	# js temp head : Create the dc.js chart objects & link to div		
	dashboard.env$jstemphead <- paste(dashboard.env$jstemphead, " var dcChart",dashboard.env$numeltgraph,"= dc.scatterPlot('#dc-chart-",dashboard.env$numeltgraph,"'); ",sep="")
		
	# js temp core	
	dashboard.env$jstempcore <- paste(dashboard.env$jstempcore, " 
	
	var dcChartDim", dashboard.env$numeltgraph," = records.dimension(function (d) {
		 return [d.",dashboard.env$matchlistcolnb [match(x, dashboard.env$listcol)],", d.",dashboard.env$matchlistcolnb [match(y, dashboard.env$listcol)],"]; 
		 }), 
	dcChartGroup", dashboard.env$numeltgraph," = dcChartDim", dashboard.env$numeltgraph,".group(); 

	var	xMin = d3.min(dataset, function(d) {return d.",dashboard.env$matchlistcolnb [match(x, dashboard.env$listcol)],";}),
		xMax = d3.max(dataset, function(d) {return d.",dashboard.env$matchlistcolnb [match(x, dashboard.env$listcol)],";});

	var	yMin = d3.min(dataset, function(d) {return d.",dashboard.env$matchlistcolnb [match(y, dashboard.env$listcol)],";}),
		yMax = d3.max(dataset, function(d) {return d.",dashboard.env$matchlistcolnb [match(y, dashboard.env$listcol)],";});

	dcChart", dashboard.env$numeltgraph,".width(",width,").height(",height,") 
	.dimension(dcChartDim", dashboard.env$numeltgraph,")
	.group(dcChartGroup", dashboard.env$numeltgraph,")
	.x(d3.scale.linear().domain([xMin- (xMax -xMin)/10, xMax+(xMax -xMin)/10]))
	.y(d3.scale.linear().domain([yMin- (yMax -yMin)/10, yMax+(yMax -yMin)/10]))
    .brushOn(true)
    .clipPadding(",clipaddingsize,")
    .xAxisLabel('",x,"')
    .yAxisLabel('",y,"')
	.symbol('", symboletype ,"')
	.symbolSize(", symbolesize ,")
	.highlightedSize(", hlightedsize ,")
	.title(function(d){return d.value;});
	
	",sep="")

			
	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
		# close the previous row, write in a row then start a new one 
		linebreak()
	}			
	
}


