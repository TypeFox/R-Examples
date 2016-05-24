## Linear Predictor Score formula builders
## Author : Sylvain Mareschal <maressyl@gmail.com>

LPS.p.form <- function(responseName, p.values, threshold=0.05, method="fdr") {
	# p-value threshold
	selection <- p.adjust(p.values, method=method) <= threshold
	selection <- names(p.values)[ selection ]
	if(length(selection) == 0) stop("No gene selected")
	
	# Formula building
	form <- sprintf("%s~%s", responseName, paste(sprintf("`%s`", selection), collapse="+"))
	form <- as.formula(form)
	
	return(form)
}

LPS.k.form <- function(responseName, coeff, k) {
	# k best statistics
	if(as.integer(k) <= 0)            stop("No gene selected")
	if(as.integer(k) > length(coeff)) stop("'k' is larger than available feature count")
	selection <- head(order(abs(coeff), na.last=TRUE, decreasing=TRUE), as.integer(k))
	selection <- names(coeff)[ selection ]
	
	# Formula building
	form <- sprintf("%s~%s", responseName, paste(sprintf("`%s`", selection), collapse="+"))
	form <- as.formula(form)
	
	return(form)
}
