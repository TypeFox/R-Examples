#' @title Computes a Genetic Relationship Matrix from a GenABEL object
#'
#' @description
#' Two alternative methods for GRM computations implemented.
#'
#' @param gen.data The GenABEL object.
#' @param method Method to be used.
#' 
#' @author Lars Ronnegard
#' 
compute.GRM <-
function(gen.data, method = "GenABEL") {
  if (method == "GenABEL"){
    gkins <- ibs(gen.data, weight = "freq")
    GRM <- gkins
    GRM[t(lower.tri(gkins))] <- 0
    rm(gkins)
    GRM <- 2*(GRM + t(GRM) - diag(diag(GRM)))
  }
  if (method == "ZZt"){
    SNP.matrix <- scale(as.double(gen.data))
    m <- ncol(SNP.matrix)
    markers.to.fit <- (1:m)[!is.na(colMeans(SNP.matrix))] 
    GRM <- tcrossprod(SNP.matrix[,markers.to.fit])/length(markers.to.fit)
    rm(SNP.matrix)
  }
  return(GRM)
}
