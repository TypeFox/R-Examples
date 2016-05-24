#' Calculate Jost's D 
#'
#' This function calculates Jost's D from a genind object
#'
#' Takes a genind object with population information and calculates Jost's D
#' Returns a list with values for each locus as well as two global estimates.
#' 'global.het' uses the averages of Hs and Ht across all loci while 
#' 'global.harm_mean' takes the harmonic mean of all loci.
#' 
#' Because estimators of Hs and Ht are used, its possible to have negative
#' estimates of D. You should treat these as numbers close to zero.
#'
#' @return per.locus values for each D for each locus in the dataset
#' @return global estimtes for D based on overall heterozygosity or the harmonic
#' mean of values for each locus
#' @param hsht_mean The type of mean to use to calculate values of Hs and Ht
#' for a global estimate. (Default is teh airthmetic mean, can also be set to
#' the harmonic mean).
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' D_Jost(nancycats)
#' D_Jost(nancycats, hsht_mean= "arithmetic")
#' @references
#'  Jost, L. (2008), GST and its relatives do not measure differentiation. Molecular Ecology, 17: 4015-4026. 
#' @family diffstat
#' @family D

D_Jost <- function(x, hsht_mean = "arithmetic"){
  mean_type <- match.arg(hsht_mean, c("arithmetic", "harmonic"))
  mean_f <- if(mean_type == "arithmetic") mean else harmonic_mean
  gn <- length(unique(pop(x))) 
  loci <- t(sapply(seploc(x), D.per.locus))
  global_Hs <- mean_f(loci[,1], na.rm=T)
  global_Ht <- mean_f(loci[,2], na.rm=T)
  global_D <-  (global_Ht - global_Hs)/(1 - global_Hs ) * (gn/(gn-1))
  harm_D <- harmonic_mean(loci[,3])
  return(list("per.locus"=loci[,3],
              "global.het"=global_D,
              "global.harm_mean" = harm_D
        ))

}

D.per.locus <- function(g) {
    hets <- HsHt(g)
    if(all(is.na(hets))){
       return(hets)
    }
    Ht_est <- hets[[1]]
    Hs_est <- hets[[2]]
    n <- hets[[3]]
    D <- (Ht_est-Hs_est)/(1-Hs_est) * (n/(n-1))
    return(c(Hs_est, Ht_est, D))
  }

