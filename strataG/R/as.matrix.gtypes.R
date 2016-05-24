#' @title Convert \code{gtypes} To \code{matrix}
#' @description Convert genotypic data in the \code{@@loci} slot of a 
#'   \linkS4class{gtypes} object to a \code{matrix}.
#'   
#' @param x a \linkS4class{gtypes} object.
#' @param one.col logical. If \code{TRUE}, then result has one column per 
#'   locus.
#' @param sep character to use to separate alleles if \code{one.col} is 
#'   \code{TRUE}.
#' @param ... additional arguments ot be passed to or from methods.
#'   
#' @return A \code{matrix} or \code{data.frame} with one row per sample.
#' 
#' @author Eric Archer \email{eric.archer@@noa.gov}
#' 
#' @seealso \link{df2gtypes}
#' 
#' @aliases as.matrix.gtypes
#' 
#' @export
#' 
setMethod("as.matrix", "gtypes", function(x, one.col = FALSE, sep = "/", ...) {
  # create matrix identifying which rows (columns) each sample (row) is in 
  #   in the @loci slot
  id.mat <- matrix(1:nrow(x@loci), ncol = ploidy(x))
  # loop through each locus
  gen.mat <- do.call(cbind, lapply(locNames(x), function(locus) {
    # extract alleles for this locus as a matrix with one allele per colum
    this.loc <- x@loci[[locus]]
    this.loc <- apply(id.mat, 1, function(i) as.character(this.loc[i]))
    this.loc <- if(is.vector(this.loc)) cbind(this.loc) else t(this.loc)  
    # collapse alleles to create a one column matrix if one.col == TRUE
    if(one.col & x@ploidy > 1) {
      this.loc <- cbind(apply(this.loc, 1, function(alleles) {
        if(any(is.na(alleles))) NA else paste(alleles, collapse = sep)
      }))
    }
    # assign column names
    colnames(this.loc) <- if(ncol(this.loc) > 1) {
      paste(locus, 1:ncol(this.loc), sep = ".")
    } else {
      locus
    }
    this.loc
  }))
  rownames(gen.mat) <- indNames(x)
  gen.mat    
})