#' Prints the summary results of the logistic regression models.
#' 
#' The \code{\link{trips}} function estimates three separate logistic regression
#' models for each pair of groups. This function will print a combined table
#' of the three summaries.
#' 
#' @param object the results of \code{\link{trips}}.
#' @param ... currently unused.
#' @method summary triangle.psa
#' @export
summary.triangle.psa <- function(object, ...) {
	m1 <- attr(object, 'model1')
	m2 <- attr(object, 'model2')
	m3 <- attr(object, 'model3')
	
	m1sum <- summary(m1)
	m2sum <- summary(m2)
	m3sum <- summary(m3)
	
	ls(m1sum)
	print(
	cbind(as.data.frame(m1sum$coefficients),
		  'm1'=star(m1sum$coefficients[,4]),
		  as.data.frame(m2sum$coefficients),
		  'm2'=star(m2sum$coefficients[,4]),
		  as.data.frame(m3sum$coefficients), 
		  'm3'=star(m3sum$coefficients[,4])),
	digits=3)

	
}

# star <- function(pvals) {
# 	ifelse(pvals < .001, '***', 
# 		   ifelse(pvals < .01, '**',
# 		   	   ifelse(pvals < .05, '*',
# 		   	   	   ifelse(pvals < .1, '.', ''))))
# }
