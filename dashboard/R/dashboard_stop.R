#' dashboard_stop stops the local Rook server
#'
#' \code{dashboard_stop} stops the Rook server running. Not needed in linux, but required in unix environment
#' 
#' @docType methods
#' @param dashboard.env name of the global environment variable used across the dashboard package
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5])
#' dcbarchart(x=names(iris)[1] , gap=75)
#' dcpiechart(x=names(iris)[2])
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files. Server is not launched
#' dashboard_stop(dashboard.env) # should have a server running
#'  

dashboard_stop <-function(dashboard.env=dashboard.env){
	# stop the server 
	if (interactive()){
	    dashboard.env$s$stop()	
	}
}

