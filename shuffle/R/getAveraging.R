getAveraging <-
function(des){
  #########
  # Convert a design (either permutation pi or matrix X)
  # into the averaging mat format:
  # D = res$B * res$ns
  # B = res$B
  # B-G = res$B - res$G
  #########
  
  # If design is a vector, convert it into a matrix
  if (is.null(dim(des))) {   
    des = designVec2Mat(des)
  }
  D = des %*% t(des)
  B = D/colSums(D)
  avgmat = list(B= B, ns = colSums(D), G = 1/nrow(D),m=ncol(des))
  # B_G = B-1/nrow(D)
  # D = B*ns
  return(avgmat)
}
