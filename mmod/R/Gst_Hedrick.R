#' Calculate Nei's Gst using estimators for Hs and Ht
#'
#' This function calculates Hedrick's G'st from a genind object
#'
#' Takes a genind object with population information and calculates Hedrick's 
#' G''st.
#' 
#' Because estimators of Hs and Ht are used, it's possible to have negative
#' estimates of G''st. You should treat such results as zeros (or 
#' an attempt to estimate a very low number with some error which might push it
#' below zero)
#' 
#'
#' @param x genind object (from package adegenet)
#' @export
#' @return per.locus values for each G''st for each locus in the dataset
#' @return global estimtes for G''st based on overall heterozygosity 
#' @references
#'  Hedrick, PW. (2005), A Standardized Genetic Differentiation Measure. Evolution 59: 1633-1638. 
#' @references
#' Meirmans PG, Hedrick PW (2011), Assessing population structure: FST and related measures. Molecular Ecology Resources, 11:5-18
#' @family diffstat
#' @family Hedrick
#' @examples 
#' data(nancycats) 
#' Gst_Hedrick(nancycats)

Gst_Hedrick <- function(x){
  gn <- length(unique(pop(x)))  

  Gst.per.locus <- function(g) {
    hets <- HsHt(g) #A private function form mmod
    Ht_est <- hets[[1]]
    Hs_est <- hets[[2]]
    n <- hets[[3]]
    Gprime_st <- n * (Ht_est - Hs_est) / ((n * Ht_est - Hs_est) * (1 - Hs_est))
    return(c(Hs_est, Ht_est, Gprime_st))
  }
 loci <- t(sapply(seploc(x), Gst.per.locus))
  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_GstH <-  gn * (global_Ht - global_Hs) / ((gn * global_Ht - global_Hs)*(1-global_Hs))
  return(list("per.locus"=loci[,3], "global"=global_GstH))

}

