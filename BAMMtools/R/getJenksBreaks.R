###############################
#
#	getJenksBreaks <- function(...)
#	Function returns the Jenks natural breaks for a vector of values
#
#		var = numeric vector
#		k = number of categories
#		subset = if not NULL, this is the number of regularly sampled values to be taken from var



getJenksBreaks <- function(var, k, subset = NULL) {
	k <- k - 1;
	
	#if more breaks than unique values, segfault, so avoid
	if (k > length(unique(var))) {
		k <- length(unique(var));
	}
	brks <- rep(1, k + 1);
	
	#if requested, regularly sample subset values
	if (!is.null(subset)) {
		if (length(var) > subset) {
			ind <- c(seq(from=1, to=length(var), by=floor(length(var)/subset)), length(var));
			var <- var[ind];
		}
	}
	
	d <- sort(var);
	length_d <- length(d);
	return(.C("jenksBrks", as.double(d), as.integer(k), as.integer(length_d), as.double(brks))[[4]]);
}
