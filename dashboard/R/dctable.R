#' dctable add a data table to the dashboard
#'
#' \code{dctable} displays a data table
#' 
#' @docType methods
#' @param index column name of data frame \code{data} for indexing the data table
#' @param listcol vector of column names of data frame \code{data} to display in the data table
#' @param title character for the title of the generated data table
#' @param spansize integer between 1 to 12 for the width of the element in the web page
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5])
#' dcbarchart(x=names(iris)[1] , gap=75)
#' dcpiechart(x=names(iris)[2])
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files. Server is not launched
#'  

dctable <- function(
index= listcol[1],
listcol=dashboard.env$listcol,
title="data list table",
spansize=12
){
	# dimension name is used by convention (DC.js)
	dimension =index 
	
	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
		# close the previous row, write in a row then start a new one 
		linebreak()
	}
	
	dashboard.env$numelttable <- dashboard.env$numelttable + 1
	dashboard.env$nbelementinrow <- dashboard.env$nbelementinrow + 1	
	
	dashboard.env$sumspaninrow <- dashboard.env$sumspaninrow + spansize

if (dashboard.env$numelttable >1){
	cat("Maybe several table in the dashboard would crash")
}	
	dashboard.env$temp <- paste(dashboard.env$temp,"<div class='span",spansize,"'><table class='table table-hover list table-striped table-bordered' id='dc-data-table",dashboard.env$numelttable,"'><thead><tr class='header'>", capture.output(cat(paste(rep("<th data-dynatable-column='",length(listcol)),dashboard.env$matchlistcolnb[match(listcol, dashboard.env$listcol)],rep("'>",length(listcol)),listcol,rep("</th>",length(listcol)), sep=""))), "</tr></thead></table></div>", sep="")
							
							
	# Create the dc.js chart objects & link to div		
	dashboard.env$jstemphead <- paste(dashboard.env$jstemphead, " var dataList",dashboard.env$numelttable,"; ",sep="")

	
	dashboard.env$jstempcore <- paste(dashboard.env$jstempcore, " var dataListDim", dashboard.env$numelttable," = records.dimension(function (d) { return d.",dashboard.env$matchlistcolnb [match(dimension, dashboard.env$listcol)],"; }); 
	
	dataList", dashboard.env$numelttable,"= $('#dc-data-table",dashboard.env$numelttable,"').dynatable({features:{pushState: false},
												  dataset: {records: dataListDim", dashboard.env$numelttable,".top(Infinity), perPageDefault: 10, perPageOptions: [10,20,50,100,200,500] }}).data('dynatable');",sep="")


if (dashboard.env$numelttable ==1){
	
dashboard.env$jstempendt1 <- "function RefreshTable() {
		dc.events.trigger(function () {"
			

dashboard.env$jstempendt3 <- "			
							   });	
			 
			 var filtered = all.value();
			 var total = records.size();
			 
			 if(filtered==total) {
				 $('a.dc-data-count-reset-all').text(' ');
				 $('span.filter-count' ).text(' ');
				 $('span.filter-count-text' ).text('All selected out of');
				 $('span.total-count' ).text(total);
				 
			 } else{
				 $('a.dc-data-count-reset-all').text('Reset All');
				 $('span.filter-count' ).text(filtered);
				 $('span.total-count' ).text(total);
				 $('span.filter-count-text' ).text('selected out of');				 
			 }
			 
		 };
		 
		 for (var i = 0; i < dc.chartRegistry.list().length; i++) {
			 var chartI = dc.chartRegistry.list()[i];
			 chartI.on('filtered', RefreshTable);
		 }
		 
		 RefreshTable();	
"

}



dashboard.env$jstempend <- paste(dashboard.env$jstempend,
"
			dataList", dashboard.env$numelttable,".settings.dataset.originalRecords = dataListDim", dashboard.env$numelttable,".top(Infinity);
			dataList", dashboard.env$numelttable,".process();",sep="")
			

			
			
			
			
	if (spansize==12 || (dashboard.env$sumspaninrow + spansize) > 12) {
		# close the previous row, write in a row then start a new one 
		linebreak()
	}			
}


