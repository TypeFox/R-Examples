#' Pairwise Euclidean distances between landmarks 
#'
#' This function computes all possible \emph{n(n-1)/2} pairwise Euclidean distances between \emph{n} landmarks.
#' @param x a two-column matrix for the xy landmark coordinates
#' @param average if TRUE, the pairwise distances for both left and right anchors are averaged. For data quality checks, set to FALSE
#' @details Comparison of the pairs of pairwise Euclidean distance from the left and right form of anchors is 
#' the basis of the quality control procedure implemented in \code{Qscore}. In addition, pairwise Euclidean distances 
#' provide length variables useful for the analysis of size variation (Lele & Richtsmeier, 2001).
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Lele SR,  Richtsmeier  JT. 2001. An Invariant Approach to Statistical Analysis of Shape. Boca Raton: Chapman and Hall.
#' @examples
#' library(cluster)
#'
#' data(ligophorus_tpsdata)
#'
#' #There are 11 landmarks for the ventral and dorsal anchors, 
#' #yielding 110 pairwise Euclidean distances
#' #The indices for the pairwise Euclidean distances map to the upper triangle 
#' #of the pairwise distance matrix by row
#' pwdist(ligophorus_tpsdata$bantingensis[[1]],average=FALSE)
#'

pwdist <- function(x,average=TRUE){

result <- vector("list",4)
for(j in 1:4){
coordinate_x <- x[ (11*j - 10) : (11*j), ]
length_x <- as.numeric (daisy(coordinate_x))
result[[j]] <- length_x
}

if(average == TRUE){
ventral <- (result[[1]]+result[[2]] )/2
dorsal <- (result[[3]]+result[[4]] )/2
return(c(ventral,dorsal))
}

else return(cbind(c(result[[1]],result[[3]]),c(result[[2]],result[[4]])))

}
