#' Obtention of the multilocus Loiselle's Fij matrix
#' 
#' @param eco Object of class ecogen.
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' loiselle <- eco.kin.loiselle(eco)
#' loiselle[1:5, 1:5]
#' 
#' }
#' 
#' @references
#' Kalisz, S., J. Nason, F.M. Handazawa, and S. Tonsor. 2001. 
#' Spatial population genetic structure in Trillium grandiflorum: 
#' the roles of dispersal, mating, history, and selection. 
#' Evolution 55: 1560-1568.
#' 
#' Loiselle, B., V. Sork, J. Nason, and C. Graham. 1995. 
#' Spatial genetic structure of a tropical understory shrub, 
#' Psychotria officinalis (Rubiaceae). 
#' American Journal of Botany 1420-1425.
#' 
#' @export

eco.kin.loiselle <- function(eco) {
  geno <- eco@A
  locus <- as.numeric(eco@INT@loc.fac)
  nloc <- max(locus)
  nal <- ncol(geno)
  N <- nrow(geno)
  alelos <- nrow(geno)
  p.la <- apply(geno, 2, sum, na.rm = TRUE)
  p.locus <- tapply(p.la,locus, sum, na.rm = TRUE)
  p.la <- p.la / p.locus[locus]
  
  p.la.cen <- t(replicate(N, p.la))
  p.center <- geno - p.la.cen
  
  n1 <- 2 * p.locus
  n1 <- n1[locus]
  
  
  ###SUPERIOR COEFFICENT###
  
  #RIGHT TERM	
  T2  <- p.la * (1 - p.la) / (n1 - 2)
  
  #LEFT TERM
  out <- list()
  for(i in 1:nal) {
    out[[i]] <- outer(p.center[,i], p.center[,i])
  }
  
  numer <- matrix(0, N , N)
  for(i in 1:nal){
    numer <- numer + ifelse(is.na(out[[i]]), 0, out[[i]] + matrix(T2[i], N, N))
  }
  
  ###INFERIOR COEFFICENT###
  
  denom <- sum(p.la * (1 - p.la))
  
  ###Fij COMPUTATION###
  Fij <- (numer / denom)
  diag(Fij) <- NA
  rownames(Fij) <- rownames(eco@G)
  colnames(Fij) <- rownames(eco@G)
  
  Fij
  
}
