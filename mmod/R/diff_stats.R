#' Calculate differentiation statistics for a genind object
#'
#' By default this function calculates three different statistics of 
#' differentiation for a genetic dataset. Nei's Gst, Hedrick's G''st and
#' Jost's D. Optionally, it can also calculate Phi'st, which is not calculated 
#' by default as it can take somewhat more time to run.
#' 
#' See individual functions (listed below) for more details.
#'
#' @param x genind object (from package adegenet)
#' @param phi_st Boolean Calculate Phi_st (default is FALSE)
#' @export
#'
#' @return per.locus values for each statistic for each locus in the dataset
#' @return global estimtes for these statistics across all loci in the dataset 
#' 
#' @references
#'  Hedrick, PW. (2005), A Standardized Genetic Differentiation Measure. Evolution 59: 1633-1638. 
#' @references
#'  Jost, L. (2008), GST and its relatives do not measure differentiation. Molecular Ecology, 17: 4015-4026.
#' @references
#' Meirmans PG, Hedrick PW (2011), Assessing population structure: FST and related measures. Molecular Ecology Resources, 11:5-18
#' @references
#'  Nei M. (1973) Analysis of gene diversity in subdivided populations. PNAS: 3321-3323. 
#' @references
#'  Nei M, Chesser RK. (1983). Estimation of fixation indices and gene diversities. Annals of Human Genetics. 47: 253-259.
#' @references
#'  Meirmans, PW. (2005), Using the AMOVA framework to estimate a standardized genetic differentiation measure. Evolution 60: 2399-402.
#' @references
#'  Excoffier, L., Smouse, P., Quattro, J. (1992), Analysis of molecular variance inferred from metric distances among DNA haplotypes: application to human mitochondrial DNA restriction data. Genetics 131: 479-91
#' @family diffstat
#' @examples
#' data(nancycats)
#' diff_stats(nancycats)

diff_stats <- function(x, phi_st=FALSE){  
  #Global mean population sizes for global stats (can vary per locus due to NAs)
  gn <- length(unique(pop(x))) 
  per.locus <- function(g) {
    hets <- HsHt(g) #A private function form mmod
    Ht_est <- hets[[1]]
    Hs_est <- hets[[2]]
    n <- hets[[3]]
    G_est <- (Ht_est-Hs_est)/Ht_est
    D <- (Ht_est-Hs_est)/(1-Hs_est) * (n/(n-1))
    Gprime_st <- n * (Ht_est - Hs_est) / ((n * Ht_est - Hs_est) * (1 - Hs_est))
    #And the results formated as list
    result <- c("Hs" = Hs_est, 
                "Ht" = Ht_est, 
                "Gst"=G_est, 
                "Gprime_st" = Gprime_st,
                "D" = D)
    return(result)
  }
 loci <- t(sapply(seploc(x), per.locus))

  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_G_est <- (global_Ht - global_Hs)/global_Ht
  
  if(phi_st){
    phi_stat <- Phi_st_Meirmans(x)
    loci <- cbind(loci, Phi_st=phi_stat$per.locus)

  }
  
  global <- c(Hs = global_Hs, 
              Ht = global_Ht, 
              Gst_est = global_G_est, 
              "Gprime_st"= gn * (global_Ht - global_Hs) / ((gn * global_Ht - global_Hs)*(1-global_Hs)),
              "D_het" = (global_Ht - global_Hs)/(1 - global_Hs ) * (gn/(gn-1)),
              "D_mean"= harmonic_mean(loci[,5]))
  
  if(phi_st){
    global <- c(global, Phi_st=phi_stat$global)
  }
  
  return(list("per.locus"=loci, global=global))
}

