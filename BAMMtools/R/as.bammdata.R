as.bammdata <- function(x, ...) {
	if (length(class(x) == 1) && class(x) == "bammdata") {
		return(x);
	}
	UseMethod("as.bammdata");
}

as.bammdata.credibleshiftset <- function(x, ...) {
	obj <- x;
	class(obj) <- "bammdata";
	return(obj);
}
