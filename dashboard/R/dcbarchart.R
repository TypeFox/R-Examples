#' dcbarchart adds a bar chart to the dashboard
#'
#' \code{dcbarchart} generates a bar chart
#' 
#' @docType methods
#' @param x column name of data frame \code{data} for drawing a bar chart
#' @param title character for the title of the generated bar chart
#' @param spansize integer between 1 to 12 for the width of the element in the row
#' @param width integer for the width (in pixel) of the element in the web page
#' @param height integer for the height (in pixel) of the element in the web page
#' @param gap integer for adjusting the width of each bar
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5])
#' dcbarchart(x=names(iris)[1] , gap=75)
#' dcpiechart(x=names(iris)[2])
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files. Server is not launched
#'  

dcbarchart <- function(
x,
title=paste(x,"bar chart"),
spansize=4,
width=dashboard.env$outerwidth*spansize/12,
height=250,
gap=width/2
){
	# dimension name is used by convention (DC.js)
	dimension = x
	
	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
	# close the previous row, write in a row then start a new one 
		linebreak()
	}
	
	dashboard.env$numeltgraph <- dashboard.env$numeltgraph + 1
	dashboard.env$nbelementinrow <- dashboard.env$nbelementinrow + 1	
	dashboard.env$sumspaninrow <- dashboard.env$sumspaninrow + spansize
	
	dashboard.env$temp <- paste(dashboard.env$temp,"<div class='span",spansize,"' id='dc-chart-",dashboard.env$numeltgraph,"'><h4>",title," <span><a class='reset' href='javascript:dcChart",dashboard.env$numeltgraph,".filterAll();dc.redrawAll();' style='display: none;'> reset </a> </span> </h4> </div>", sep="")	
	
	
	# js temp head : Create the dc.js chart objects & link to div		
	dashboard.env$jstemphead <- paste(dashboard.env$jstemphead, " var dcChart",dashboard.env$numeltgraph,"= dc.barChart('#dc-chart-",dashboard.env$numeltgraph,"'); ",sep="")
	
	# js temp core	
	dashboard.env$jstempcore <- paste(dashboard.env$jstempcore, " 
	
	var dcChartDim", dashboard.env$numeltgraph," = records.dimension(function (d) { 
		return d.",dashboard.env$matchlistcolnb [match(dimension, dashboard.env$listcol)],";
		}), 
	dcChartGroup", dashboard.env$numeltgraph," = dcChartDim", dashboard.env$numeltgraph,".group().reduceCount(function (d) { 
		return d.",dashboard.env$matchlistcolnb [match(dimension, dashboard.env$listcol)],"; 
	}), 
	xMin = d3.min(dataset, 	function (d) { return d.",dashboard.env$matchlistcolnb [match(dimension, dashboard.env$listcol)],"; }),
	xMax = d3.max(dataset, 	function (d) { return d.",dashboard.env$matchlistcolnb [match(dimension, dashboard.env$listcol)],"; });
	dcChart", dashboard.env$numeltgraph,".width(",width,").height(",height,") .margins(margin) 
	.dimension(dcChartDim", dashboard.env$numeltgraph,")
	.group(dcChartGroup", dashboard.env$numeltgraph,")
	.transitionDuration(500)	
	.gap(",gap,")
	.x(d3.scale.linear().domain([xMin- (xMax -xMin)/10, xMax+(xMax -xMin)/10]))
	.elasticY(true)
	.xAxis().tickFormat(); 
	
	",sep="")
	
	
	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
	# close the previous row, write in a row then start a new one 
		linebreak()
	}			
	
}


