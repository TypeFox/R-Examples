#' @name TiTvRatio
#' @title Transition / Transversion Ratio
#' @description Calculate transition/transversion ratio. Test substitution 
#'   type of two bases.
#' 
#' @param x a \linkS4class{gtypes} object with aligned sequences or a list of 
#'   aligned DNA sequences.
#' @param b1,b2 two bases to be compared.
#' 
#' @return \code{TiTvRatio}: a vector providing the number of 
#'   transitions (\code{Ti}), transversions (\code{Tv}), and the 
#'   transition/transversion ratio (\code{Ti.Tv.ratio}).\cr
#' \code{subType}: either "ti" for transition, or "tv" for transversion.\cr
#' \code{isTi} and \code{isTv}: a logical identifying whether 
#'   the \code{b1} to \code{b2} is a transition or transversion.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.seqs)
#' 
#' TiTvRatio(dolph.seqs)
#' 
#' subType("a", "c")
#' 
#' isTi("a", "c")
#' 
#' isTv("a", "c")
#' 
#' @importFrom utils combn
#' @export
#' 
TiTvRatio <- function(x) {
  site.freqs <- baseFreqs(x, c("a", "c", "g", "t"))$site.freqs
  ti.tv <- apply(site.freqs, 2, function(freqs) {
    freqs <- freqs[freqs > 0]
    if(length(freqs) < 2) {
      c(Ti = 0, Tv = 0)
    } else {
      pairs <- combn(names(freqs), 2)
      c(Ti = sum(apply(pairs, 2, function(bases) isTi(bases[1], bases[2]))),
        Tv = sum(apply(pairs, 2, function(bases) isTv(bases[1], bases[2])))
      )
    }
  })
  
  Ti = sum(ti.tv["Ti", ])
  Tv = sum(ti.tv["Tv", ])
  c(Ti = Ti, Tv = Tv, Ti.Tv.ratio = Ti / Tv)
}

#' @rdname TiTvRatio
#' @export
#' 
subType <- function(b1, b2) {
  b1 <- tolower(b1)
  b2 <- tolower(b2)
  if(!(all(c(b1, b2) %in% colnames(ti.tv.mat)))) return(NA)
  ti.tv.mat[b1, b2]
}

#' @rdname TiTvRatio
#' @export
#' 
isTi <- function(b1, b2) {
  x <- subType(b1, b2) 
  if(is.na(x)) return(FALSE)
  x == "ti"
}

#' @rdname TiTvRatio
#' @export
#' 
isTv <- function(b1, b2) {
  x <- subType(b1, b2)
  if(is.na(x)) return(FALSE)
  x == "tv"
}
  
