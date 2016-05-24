#' @name gtypes2genind
#' @title Convert Between \code{gtypes} And \code{genind} objects.
#' @description Convert a \code{gtypes} object to a \code{genind} object 
#'   and vice-versa.
#' 
#' @param x either a \linkS4class{gtypes} or \linkS4class{genind} object
#'   to convert from.
#' @param type a character string indicating the type of marker for 
#'   \linkS4class{genind} objects: 'codom' stands for 'codominant' 
#'   (e.g. microstallites, allozymes); 'PA' stands for 'presence/absence' 
#'   markers (e.g. AFLP, RAPD).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, \link{sequence2gtypes}, 
#'   \link{gtypes2df}, \link{gtypes2loci}
#' 
#' @export
#' 
gtypes2genind <- function(x, type = c("codom", "PA")) {
  df2genind(X = as.matrix(x, one.col = TRUE, sep = "/"),
            sep = "/", 
            pop = strata(x),
            NA.char = NA,
            ploidy = ploidy(x),
            type = match.arg(type)
  )
}


#' @rdname gtypes2genind
#' @export
#' 
genind2gtypes <- function(x) {
  gen.mat <- genind2df(x, usepop = TRUE, oneColPerAll = TRUE)
  gen.mat[gen.mat == "NA"] <- NA
  has.pop <- !is.null(x@pop)
  df2gtypes(x = gen.mat,
            ploidy = x@ploidy[1],
            id.col = NULL,
            strata.col = if(has.pop) 1 else NULL,
            loc.col = if(has.pop) 2 else 1,
            other = other(x)
  )  
}
