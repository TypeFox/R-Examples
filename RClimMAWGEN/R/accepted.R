NULL
#'
#' This functions lists the realizations which pass successfully Ks or Wilcoxon test.  
#' 
#' @title Which generations pass the tests with success? 
#' 
#' @param tests list of objects returned by \code{\link{wilcox.test}} and \code{\link{ks.test}}
#' @param significance significance for statistical tests (maximum accepted \code{p-Value}). Default is 0.05. 
#' @return Vector with names of successful realizations.
#' 
#' @export
#' @seealso \code{\link{climdex.data.frame}},\code{\link{ks.test}},\code{\link{ks.test.climdex.data.frame}},\code{\link{wilcox.test}}
#' @examples
#' 
#' # See the example of 'climdex.data.frame' function
#'  




accepted <- function(tests,significance=0.05) {
	
	out <- names(tests)
	
	accepted <- array(FALSE,length(out))

	for (i in 1:length(tests)) {
		
		  if (tests[[i]]$p.value>significance) accepted[i] <- TRUE
		
	}
	
	out <- out[accepted]
	
	return(out)
	
}