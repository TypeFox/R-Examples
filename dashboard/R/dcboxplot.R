#' dcboxplot adds a box plot to the dashboard
#'
#' \code{dcboxplot} generates a box plot
#' 
#' @docType methods
#' @param x column name of categorical variable of data frame \code{data} for drawing a box plot. One box is drawed for each distinct value of \code{x}
#' @param val column name of a single numerical variable in the data frame \code{data} for computing the size of the box 
#' @param title character for the title of the generated box plot
#' @param spansize integer between 1 to 12 for the width of the element in the row
#' @param width integer for the width in pixel of the element in the web page
#' @param height integer for the height in pixel of the element in the web page
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5]) 
#' dcboxplot(x=names(iris)[5], val=names(iris)[3] )
#' dcbarchart(x=names(iris)[1] , gap=75)
#' dcpiechart(x=names(iris)[2])
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files. Server is not launched
#'  

dcboxplot <- function(
x,
val,
title=paste(x,"boxplot: ",val),
spansize=4,
width=dashboard.env$outerwidth*spansize/12,
height=250
){
	# dimension name is used by convention (DC.js)
	dimension = x
	groupval = val

	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
		# close the previous row, write in a row then start a new one 
		linebreak()
	}
	
	dashboard.env$numeltgraph <- dashboard.env$numeltgraph + 1
	dashboard.env$nbelementinrow <- dashboard.env$nbelementinrow + 1	
	dashboard.env$sumspaninrow <- dashboard.env$sumspaninrow + spansize
		
	dashboard.env$temp <- paste(dashboard.env$temp,"<div class='span",spansize,"' id='dc-chart-",dashboard.env$numeltgraph,"'><h4>",title," <span><a class='reset' href='javascript:dcChart",dashboard.env$numeltgraph,".filterAll();dc.redrawAll();' style='display: none;'> reset </a> </span> </h4> </div>", sep="")	
	
			
	# js temp head : Create the dc.js chart objects & link to div		
	dashboard.env$jstemphead <- paste(dashboard.env$jstemphead, " var dcChart",dashboard.env$numeltgraph,"= dc.boxPlot('#dc-chart-",dashboard.env$numeltgraph,"'); ",sep="")
		
	# js temp core	
	dashboard.env$jstempcore <- paste(dashboard.env$jstempcore, " 
	
	var dcChartDim", dashboard.env$numeltgraph," = records.dimension(function (d) {
		 return d.",dashboard.env$matchlistcolnb [match(dimension, dashboard.env$listcol)],"; 
		 }), 
	dcChartGroup", dashboard.env$numeltgraph," = dcChartDim", dashboard.env$numeltgraph,".group().reduce(
		function(p,v){
			p.push(v.",dashboard.env$matchlistcolnb [match(groupval, dashboard.env$listcol)],");
			return p;
		},
		function(p,v){
			p.slice(p.indexOf(v.",dashboard.env$matchlistcolnb [match(groupval, dashboard.env$listcol)],"),1);
			return p;
		},
		function(){
			return[];
		}
	); 
	var	yMin = d3.min(dataset, 	function (d) { return d.",dashboard.env$matchlistcolnb [match(groupval, dashboard.env$listcol)],"; }),
		yMax = d3.max(dataset, 	function (d) { return d.",dashboard.env$matchlistcolnb [match(groupval, dashboard.env$listcol)],"; });

	dcChart", dashboard.env$numeltgraph,".width(",width,").height(",height,").margins({top: 10, right:50, bottom:30, left:50})
	.dimension(dcChartDim", dashboard.env$numeltgraph,")
	.group(dcChartGroup", dashboard.env$numeltgraph,")
	.y(d3.scale.linear().domain([yMin- (yMax -yMin)/10, yMax+(yMax -yMin)/10]));
	",sep="")
			
	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
		# close the previous row, write in a row then start a new one 
		linebreak()
	}			
	
}


