#' @name pkg_prep
#' @title PkgName Prep Functions
#' @keywords internal
new.prelimstats <- function() {
	new <- matrix()
	class(new) <- c("drsmooth", "matrix")
	return(new)
}
#' @rdname pkg_prep
#' @keywords internal
validate.prelimstats <- function(x) {
    allowed_types <- list("numeric", "integer")
    
	for (i in 1:ncol(x)) {
       if (!isTRUE(class(x[,i]) %in% allowed_types)) stop ("All input data other than column headers must be of type numeric or integer.")

	}
}

