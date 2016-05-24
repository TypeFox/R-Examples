#' Wrapper for constructor of Bed class
#' 
#' @param bed_path character. Path to bed file.
#' @param bytes_snp integer. Bytes per SNP.
#' @param nindiv integer. Number of individuals.
#' @examples 
#' ## do not run
#' # rbed_obj = rBed("test.bed")
#' # geno = rbed_obj$readBed()
#' # geno = getJArray(geno)
#' # print(geno)
#' @return jobjRef object.
#' @import rJava
#' @author Kaiyin Zhong
#' @export
rBed = function(bed_path, bytes_snp = NULL, nindiv = NULL) {
	if(is.null(bytes_snp) && is.null(nindiv)) {
		.jnew("vu/co/kaiyin/Bed", bed_path)
	} else {
		.jnew("vu/co/kaiyin/Bed", bed_path, bytes_snp, nindiv)
	}
}



#' Truncate n bytes from end of file
#' 
#' @param filename character. Filename.
#' @param len numeric. Number of bytes to truncate
#' @examples 
#' \dontrun{
#' fn = tempfile()
#' f = file(fn, "wb")
#' writeBin("a", f)
#' writeBin("b", f)
#' writeBin("c", f)
#' close(f)
#' file.info(fn)$size == 6
#' truncateEndOfFile(fn, 1)
#' file.info(fn)$size == 5
#' }
#' 
#' @author Kaiyin Zhong
#' @export
truncateEndOfFile = function(filename, len) {
	stopifnot(is.character(filename)
					&& (length(filename) == 1)
					&& (is.numeric(len))
					&& (length(len) == 1))
	# "V" for void, return signature
	.jcall("vu/co/kaiyin/Utils", "V", "truncateFromEnd", filename, as.integer(len))
}

#' Import Java array into R 
#' 
#' A thin wrapper around \code{rJava::.jevalArray}
#' 
#' @param mat_ref Reference object of the Java array
#' @param na_vals NA code. Default to -9.
#' 
#' @author Kaiyin Zhong
#' @export
getJArray <- function(mat_ref, na_vals = -9) {
	res = .jevalArray(mat_ref, simplify = TRUE)
	res[res %in% na_vals] = NA
	res = as.data.frame(res)
	res
}






