#' @title Duplicate Genotypes
#' @description Identify duplicate or very similar genotypes.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param num.shared either number of loci or percentage of loci two 
#'   individuals must share to be considered duplicate individuals.
#' @param num.cores number of CPU cores to use.
#' 
#' @return a data.frame with the following columns:
#' \tabular{ll}{
#'   \code{ids.1, ids.2} \tab sample ids.\cr
#'   \code{strata1, strata2} \tab sample strata.\cr
#'   \code{num.loci.genotyped} \tab number of loci genotyped for both 
#'     samples.\cr
#'   \code{num.loci.shared} \tab number of loci shared between both samples.\cr
#'   \code{prop.loci.shared} \tab proportion of loci genotyped for both samples 
#'     that are shared.\cr
#'   \code{mismatch.loci} \tab loci where the two samples do not match.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
dupGenotypes <- function(g, num.shared = 0.8, num.cores = 1) {
  #if not already, convert num.shared to %
  if(num.shared > 1) num.shared <- num.shared / nLoc(g) 
    
  shared.locs <- propSharedLoci(g, type = "ids", num.cores = num.cores)
  dup.df <- shared.locs[shared.locs[, "prop.same"] >= num.shared, ]
  if(nrow(dup.df) > 0) {
    dup.df$pct.loci.shared <- dup.df$num.same / dup.df$num.not.missing
    dup.df$strata.1 <- strata(g)[dup.df$ids.1]
    dup.df$strata.2 <- strata(g)[dup.df$ids.2]
    dup.df$mismatch.loci <- sapply(1:nrow(dup.df), function(i) {
      is.same <- as.matrix(dup.df[i, locNames(g)])[1, ] < 1
      paste(names(is.same)[is.same], collapse = ", ")
    })
    colnames(dup.df)[c(3:5)] <- c(
      "num.loci.shared", "num.loci.genotyped", "prop.loci.shared"
    )
    dup.df <- dup.df[, c("ids.1", "ids.2", "strata.1", "strata.2", 
                         "num.loci.genotyped", "num.loci.shared", 
                         "prop.loci.shared", "mismatch.loci")]
  } 
  
  if(nrow(dup.df) > 0) {
    sort.order <- order(dup.df$prop.loci.shared, dup.df$num.loci.shared, 
                        rev(dup.df$ids.1), rev(dup.df$ids.2), decreasing = TRUE
    )
    dup.df <- dup.df[sort.order, ]
    rownames(dup.df) <- NULL
  } else dup.df <- NULL
  
  dup.df
}
