

#' Run the interAdapt shiny application
#'
#' \code{runInterAdapt} Runs the interactive shiny application
#' 
#' @export
#' @import shiny RCurl mvtnorm knitr knitcitations
runInterAdapt<-function(){
	cat('starting application...\n')
	rtn <- shiny::runApp(system.file('shinyInterAdaptApp', package='interAdapt'))
	if(rtn > 0)
		cat('\napplication stopped\n')
}

