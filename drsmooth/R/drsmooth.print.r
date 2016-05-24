#' @name drsmooth.print
#' @title Print Formatted Matrices for drsmooth
#' @param m  test output matrices
#' @keywords internal

drsmooth.print <- function(m) {
	dimnames(m) <- list(rep("", dim(m)[1]), rep("", dim(m)[2]))
	print(m, justify="right", na.print="", quote=F)
}

