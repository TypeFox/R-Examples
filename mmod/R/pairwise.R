#Let's DRY out the pairwise functions but writing a general case that 
#all others (below) can make use of.

pairwise_fxn <- function(x, linearized=FALSE, FXN, VAL, ...){
    n.pops <- nPop(x)
    #all combinations 
    res <- utils::combn(1:n.pops, 2, function(p) FXN(x[pop=p], ...)[[VAL]])
    attributes(res) <- list(class="dist", Diag=FALSE, Upper=FALSE, 
                            Labels=popNames(x),Size=n.pops)
    if(linearized){
        return(res/(1-res))
    }
    res
}

#' Calculates pairwise values of Jost's D
#'
#' This function calculates Jost's D, a measure of genetic
#' differentiation, between all combinations of populaitons
#' in a genind object.
#' @param hsht_mean type of mean to use for the global estimates of Hs and Ht 
#' default it "arithmetic", can also be set to "harmonic". 
#' @param x genind object (from package adegenet)
#' @param linearized logical, if TRUE will turned linearized D (1/1-D)
#' @export
#' @examples
#' 
#' data(nancycats)
#' pairwise_D(nancycats[1:26,])
#' @return  A distance matrix with between-population values of D
#' @references
#'  Jost, L. (2008), GST and its relatives do not measure differentiation. Molecular Ecology, 17: 4015-4026. 
#' @family pairwise
#' @family D

pairwise_D <- function(x, linearized=FALSE, hsht_mean="arithmetic") {
    pairwise_fxn(x, linearized, D_Jost, "global.het", hsht_mean=hsht_mean)
}

#' Calculates pairwise values of Hedrick's G'st
#'
#' This function calculates Hedrick's G'st, a measure of genetic
#' differentiation, between all combinations of populaitons
#' in a genind object.
#'
#' @return A distance matrix with between-population values of G''st
#' @param x genind object (from package adegenet)
#' @param linearized logical, if TRUE will turned linearized G'st (1/()1-G'st))
#' @export
#' @examples
#' 
#' data(nancycats)
#' pairwise_Gst_Hedrick(nancycats[1:26,])
#' @references
#'  Hedrick, PW. (2005), A Standardized Genetic Differentiation Measure. Evolution 59: 1633-1638. 
#' @family pairwise
#' @family Hedrick

pairwise_Gst_Hedrick<- function(x, linearized=FALSE) {
    pairwise_fxn(x, linearized, Gst_Hedrick, "global")
}

#' Calculates pairwise values of Nei's Gst
#'
#' This function calculates Nei's Gst, a measure of genetic
#' differentiation, between all combinations of populaitons
#' in a genind object.
#'
#' @return dist A distance matrix with between-population values of Gst
#' @param x genind object (from package adegenet)
#' @param linearized logical, if TRUE will turned linearized Gst (1/(1-Gst))
#' @export
#' @examples
#' 
#' data(nancycats)
#' pairwise_Gst_Nei(nancycats[1:26,])
#' @references
#'  Nei M. (1973) Analysis of gene diversity in subdivided populations. PNAS: 3321-3323. 
#' @references
#'  Nei M, Chesser RK. (1983). Estimation of fixation indices and gene diversities. Annals of Human Genetics. 47: 253-259.
#' @family pairwise
#' @family Nei

pairwise_Gst_Nei <- function(x, linearized=FALSE) {
    pairwise_fxn(x, linearized, Gst_Nei, "global")
}
