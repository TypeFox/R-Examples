#' dashboard_open innitializes a new dashboard
#'
#' \code{dashboard_open} sets the configuration for the web page
#' 
#' @docType methods
#' @param data data frame used for drawing a dashboard
#' @param title character for the title of the generated dashboard
#' @param filename character for the name of the generated html file
#' @param pathoutput character for the output path of generated files
#' @param outerwidth integer for the outer width (in pixel) of the web page 
#' @param outerheight integer for the outer height (in pixel) of the web page
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5])
#' dcbarchart(x=names(iris)[1] , gap=75)
#' dcpiechart(x=names(iris)[2])
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files. Server is not launched
#'  

dashboard_open <- function( 
data, 
title="Dashbord test",
filename = "dashboard",
pathoutput=tempdir(),
outerwidth = 960,
outerheight = 700) {
	
	# Dependency on Rook package for runing a local server
	#require(Rook)

	# Here to propose to use a predefine output folder path
	folder_output <- paste(filename,as.integer(Sys.time()),sep="_")

	path_folder_output <- paste(pathoutput, folder_output, sep="/")			

	dir.create(paste(pathoutput, folder_output, sep="/"), showWarnings=FALSE)

	dir.create(paste(pathoutput, folder_output, "lib", sep="/"), showWarnings=FALSE)
		
    listnumeric <- sapply(data, is.numeric)
	listvarfloat<-names(listnumeric[listnumeric==TRUE])
	listchar <- names(listnumeric[listnumeric==FALSE])
	listcol<-names(data)

	# Need to remap names due to . or special character in the column name	
	matchlistcolnb <-paste(rep("vard",length(listcol)),c(1:length(listcol)),sep="")
	#matchlistcolnb [match(listvarfloat, listcol)]
	
	# For optimisation keep a list of variable really used
	listvarused <- c()

	# May need to change this comand to the dashboard_end in order to optimize it : remove not needed columns
    write.csv(data, file = paste(path_folder_output,"/", filename,".csv",sep="") ,row.names=FALSE, na="")
		
	# Generate the core javascript file
	jstempcore <-paste("d3.csv('",filename,".csv', function(error, data) {if (error) { console.log(error);} else { data.forEach(function(d) { ",sep="")

	# Define the numeric variables		
	if(length(listvarfloat)>0){
		jstempcore <-paste(jstempcore, capture.output(cat(paste(rep("d.",length(listvarfloat)), matchlistcolnb [match(listvarfloat, listcol)],rep("= +d['",length(listvarfloat)), listvarfloat,rep("']; ",length(listvarfloat)), sep=""))),sep="")
 	}
 
 	# Define the character variables
 	if(length(listchar)>0){
jstempcore <-paste(jstempcore, capture.output(cat(paste(rep("d.",length(listchar)), matchlistcolnb [match(listchar, listcol)],rep("=d['",length(listchar)), listchar,rep("']; ",length(listchar)), sep=""))),sep="")
 	}

	# End of the core javacript file for the setting
	jstempcore <-paste(jstempcore, paste("});dataset = data;displayData();}}); function displayData() { var margin = {top: 10, right: 10, bottom: 20, left: 40}, padding = {top: 2, right: 2, bottom: 2, left: 2}, outerWidth =", outerwidth,", outerHeight =", outerheight ,", width = outerWidth - margin.left - margin.right, height = outerHeight - margin.top - margin.bottom; var records = crossfilter(dataset); var all = records.groupAll(); dc.dataCount('.dc-data-count').dimension(records).group(all);",sep=""),sep="")
	
	#save all parameters in global variable
	dashboard.env <<- new.env()
	
	# Assign global variable in dashboard package environement
	dashboard.env$filename <- filename
	dashboard.env$folder_output <- folder_output
	dashboard.env$path_folder_output <- path_folder_output
	dashboard.env$pathoutput <- pathoutput
	dashboard.env$listvarfloat <- listvarfloat
	dashboard.env$listchar <- listchar
	dashboard.env$outerwidth <- outerwidth
	dashboard.env$outerheight <- outerheight	
	dashboard.env$numelttable <- 0
	dashboard.env$numeltgraph <- 0
	dashboard.env$nbelementinrow <- 0
	dashboard.env$rowopen <- TRUE 
	dashboard.env$listcol <- listcol
	dashboard.env$sumspaninrow <- 0
	dashboard.env$jstemphead <- "var dataset;"
	dashboard.env$jstempcore <- jstempcore
	dashboard.env$jstempend <-	" "
	dashboard.env$jstempendt1 <-	" "
	dashboard.env$jstempendt3 <-	" "
	dashboard.env$title <- title 
	dashboard.env$matchlistcolnb <- matchlistcolnb
	
	"setting ok"
}



