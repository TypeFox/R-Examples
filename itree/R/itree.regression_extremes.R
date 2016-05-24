#alg 3/3/2012
# Mostly a copy from rpart.anova.s
# but also does checking to see if we're looking for
# high means or low means. If nothing specified,
# use high means. 1 <==> high mean, -1 <==> low mean
itree.regression_extremes <- function(y, offset, parms, wt) {

    if (!is.null(offset)){
		y <- y-offset
	}

	#figure out if we want high or low means. First check if anything has been passed
	# and if not, assume high means.
	meantype <- ''
	if(missing(parms)){
		warning("Nothing specified, so assuming high means.")
		meantype <- 1
	}
    else {
		parms <- as.list(parms)
		#allow a user to pass "high" or "low" or list with 1st element = high or low, but nothing else

		# check only one item was passed
		if(length(parms) > 1){
			stop("For method='extremes' with regression must specify parms=1 or parms=-1 for \\
			high and low means respectively")
		}

		# we coerced to a list and that list has length=1, so get 1st element
		meantype <- parms[[1]]

		# print error if meantype isn't correctly specified.
		if(meantype != 1 && meantype != -1){
			stop("For method='extremes' must specify parms=1 or parms=-1 for \\
			high and low means respectively")
		}
	}

	#now return a list with required entries. highmeans = +1 / lowmeans = -1
   list(y=y, parms=meantype, numresp=1L, numy=1L,
	 summary= function(yval, dev, wt, ylevel, digits ) {
	     paste("  mean=", formatg(yval, digits),
		   ", MSE=" , formatg(dev/wt, digits),
		   sep='')
	     },
	 text= function(yval, dev, wt, ylevel, digits, n, use.n ) {
	     if(use.n) {paste(formatg(yval,digits),"\nn=", n,sep="")} else
	               {paste(formatg(yval,digits))}}
	     )
}
