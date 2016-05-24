#' @name gtypes2loci
#' @title Convert Between \code{gtypes} And \code{loci} objects.
#' @description Convert a \code{gtypes} object to a \code{\link[pegas]{loci}} object.
#' 
#' @param x a \linkS4class{gtypes} or \code{loci} formatted object.
#' @param description a label for the \code{gtypes} object (optional).
#'  
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, \link{sequence2gtypes}, 
#'   \link{gtypes2df}, \link{gtypes2genind}
#' 
#' @importFrom pegas as.loci
#' @export
#' 
gtypes2loci <- function(x) {
  mat <- as.matrix(x, one.col = TRUE, sep = "/")
  mat <- data.frame(cbind(as.character(strata(x)), mat))
  mat <- data.frame(lapply(mat, factor))
  as.loci(mat, col.pop = 1)
}

#' @rdname gtypes2loci
#' @export
#' 
loci2gtypes <- function(x, description = NULL) {
  lc <- attr(x, "locicol")
  x <- as.data.frame(x)
  loci <- do.call(cbind, lapply(colnames(x)[lc], function(l.name) {
    locus <- x[, l.name]
    locus <- strsplit(as.character(locus), split = "/")
    max.len <- max(sapply(locus, length))
    locus <- lapply(locus, function(l) if(length(l == max.len)) l else rep(l, max.len))
    locus <- do.call(rbind, locus)
    colnames(locus) <- paste(l.name, 1:ncol(locus), sep = ".")
    locus
  }))
  gen.mat <- cbind(ids = rownames(x), x[, 1], loci)
  df2gtypes(gen.mat, ploidy = 2, description = description)
}
  