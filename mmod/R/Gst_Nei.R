#' Calculate Nei's Gst using estimators for Hs and Ht
#'
#' This function calculates Gst following Nei's method and using
#' Nei and Chesser's estimators for Hs and Ht
#'
#' @param x genind object (from package adegenet)
#' @return per.locus estimates of Gst for each locus in the dataset
#' @return per.locus estimates of Gst for across all loci
#' @references
#'  Nei M. (1973) Analysis of gene diversity in subdivided populations. PNAS: 3321-3323. 
#' @references
#'  Nei M, Chesser RK. (1983). Estimation of fixation indices and gene diversities. Annals of Human Genetics. 47: 253-259.
#' @family diffstat
#' @family Nei
#' @export
#' @examples
#' 
#' data(nancycats)
#' Gst_Nei(nancycats)

Gst_Nei <- function(x){
  
  Gst.per.locus <- function(g) { 
    hets <- HsHt(g) #A private function form mmod
    Ht_est <- hets[1]
    Hs_est <- hets[2]
    G_est <- (Ht_est-Hs_est)/Ht_est
    return(c(Hs_est, Ht_est, G_est))
  }
  
 loci <- t(sapply(seploc(x), Gst.per.locus))
  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_G_est <-  (global_Ht - global_Hs)/global_Ht
  
  return(list("per.locus"=loci[,3], "global"=global_G_est))

}

