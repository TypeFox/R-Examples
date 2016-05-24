
#' Recombination frequencies computation
#' 
#' Calculate recombination frequencies for a whole matrix 
#'  
#' @param input Matrix of genotypes. Rows represent markers. 
#' Columns represent samples.
#' @param ploidy Ploidy level of the organism.
#' @param sparse Logical, if the matrix is a sparse matrix or not.
#' @param ... arguments are forwarded to \code{pairwRF}.
#' @return Matrix of pairwise recombination frequencies.
#' @examples
#' data(simTetra)
#' simTetraGen <- bases2genotypes(simTetra, ploidy = 4)
#' calcRec(simTetraGen, 4)
#' @export
calcRec <- function(input, ploidy, sparse = FALSE, ...){
  nMark <- nrow(input)
  com <- utils::combn(1:nMark, 2)
  vec <- apply(com, MARGIN = 2, FUN = function(x) pairwRF(input[x, ], ploidy, ...))  
  if(sparse){
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix needed for this function to work. Please install it.",
           call. = FALSE)
    }
    mat <- Matrix::Matrix(0, nrow = nMark, ncol = nMark, sparse = TRUE) 
    mat[lower.tri(mat)] <- 0.5 - vec
  }else{
    mat <- matrix(0.5, nrow = nMark, ncol = nMark)    
    ind <- upper.tri(mat)
    mat[lower.tri(mat)] <- vec    
    mat[ind] <- t(mat)[ind]
    diag(mat) <- 0
  }
  rownames(mat) <- colnames(mat) <- rownames(input)
  return(mat)
}

#
#' Pairwise recombination frequency calculation
#' 
#' Calculates the pairwise recombination frequencies for two given markers.
#' Heuristic approach, assuming minimal recombination between markers.  
#'  
#' @param input Matrix with two rows of genotypes. Each row is one marker.
#' @param ploidy Ploidy level of the organism.
#' @return Numeric value of recomination between the two markers.
#' @keywords internal
pairwRF <- function(input, ploidy = 4, na.penalty = 0){
  n <- ncol(input)
  a1 <- input[1, ]
  a2 <- input[2, ]
  nas <- sum(is.na(a1 + a2))
  aa <- sum(abs(a1 - a2), na.rm = T) + nas * na.penalty
  ab <- sum(abs(a1 - (ploidy - a2)), na.rm = T) + nas * na.penalty  
  if(!(aa + ab) == 0){
    r <- min(aa, ab)
    theta <- r / (aa + ab)
    return(theta)
  }else{
    return(NA)
  }
}