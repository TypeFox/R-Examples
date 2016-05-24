# genepop.to.genind

#' Imports a \code{.txt} file in \code{GenePop} format into an object of type \code{genind}

#' 
#' The main work is done by the function \code{adegenet::read.genepop}. However, that function requires text files with an extension of \code{.gen}, whereas such files usually have extension \code{.txt}. The sole purpose of this function is to work around the \dQuote{.gen} requirement.
#' 
#' @param name the name of a file in \code{GenePop} format
#' @param quiet whether a conversion message should be printed
#' @param ncode Set to the number of characters per allele name
#' 
#' @return an object of class \code{genind}
#' 
#' @export
genepop.to.genind <- function(name, quiet = TRUE, ncode = 3) {
	if (requireNamespace("adegenet")) {
		tempfile <- file(name)
		tmp <- readLines(tempfile)
		writeLines(tmp, "tempgenepop.gen")
		ind <- adegenet::read.genepop("tempgenepop.gen", quiet = quiet, ncode = ncode)
		unlink("tempgenepop.gen")
		return(ind)
	} else {
		return(0)
	}
}