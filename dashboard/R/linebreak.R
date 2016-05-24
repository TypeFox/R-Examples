#' linebreak adds a line break in the dashboard
#'
#' \code{linebreak} generates a line break in the dashboard  
#' 
#' @docType methods
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5])
#' linebreak()
#' dcbarchart(x=names(iris)[1] , gap=75)
#' linebreak()
#' dcpiechart(x=names(iris)[2])
#' linebreak()
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files. Server is not launched
#'  
 
linebreak <- function(){
	if (dashboard.env$rowopen == TRUE && dashboard.env$nbelementinrow == 0) {
		cat("row already ready")
		dashboard.env$sumspaninrow <- 0	
	} else if (dashboard.env$rowopen == TRUE)  {
		dashboard.env$temp <- paste(dashboard.env$temp, "</div><div class='row'>", sep="")
		dashboard.env$nbelementinrow <- 0
		dashboard.env$sumspaninrow <- 0			
	}
	else {
		dashboard.env$temp <- paste(dashboard.env$temp, "<div class='row'>", sep=" ")
		dashboard.env$rowopen <- TRUE
		dashboard.env$nbelementinrow <- 0
		dashboard.env$sumspaninrow <- 0		
	}
}

