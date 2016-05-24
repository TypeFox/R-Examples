#' @title Nei's Da
#' @description Calcuate frequency-based Nei's Da for haploid or diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @details Returns Nei's Da for each pair of strata.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Nei et al 1983 Accuracy of Estimated Phylogenetic Trees from 
#'   Molecular Data. J Mol Evol 19:153-170 (eqn 7)\cr
#'   Nei, M., and S. Kumar (2000) Molecular Evolution and Phylogenetics. 
#'   Oxford University Press, Oxford. (pp. 268, eqn 13.6)
#'
#' @examples
#' data(msats.g)
#' 
#' neiDa(msats.g)
#'  
#' @export
#' 
neiDa <- function(g) {
  st.pairs <- as.matrix(.strataPairs(g))
  st.col <- rep(strata(g), ploidy(g))

  Da <- apply(st.pairs, 1, function(sp) {
    loc.sum <- sapply(1:nLoc(g), function(i) {
      locus <- g@loci[, i]
      to.use <- st.col %in% sp & !is.na(st.col) & !is.na(locus)
      freqs <- prop.table(table(locus[to.use], droplevels(st.col[to.use])))
      sum(apply(freqs, 1, function(f) {
        if(all(f == 0)) NA else sqrt(prod(f))
      }), na.rm = TRUE)
    })
    1 - sum(loc.sum, na.rm = TRUE) / sum(!is.na(loc.sum))
  })
  cbind(data.frame(st.pairs), Da = Da)
}