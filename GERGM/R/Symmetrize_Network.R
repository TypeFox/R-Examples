Symmetrize_Network <- function(network){
  if(isSymmetric(network) == FALSE){
    lower_to_upper <- function(m) {
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
      m
    }

    upper_to_lower <- function(m) {
      m[lower.tri(m)] <- t(m)[lower.tri(m)]
      m
    }

    if(sum(upper.tri(network)) > sum(lower.tri(network))){
      network <- upper_to_lower(network)
    }else{
      network <- lower_to_upper(network)
    }
    return(network)
  }else{
    return(network)
  }
}
