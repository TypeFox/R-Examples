# Function to join the duplicated species (i.e. with more than one polygon/shapefile) in the presence-absence matrix
# Bruno Vilela

.unicas <- function(resu) {
  nomes <- colnames(resu)
  
  if(!any(duplicated(nomes))) {
    return(resu)
  } else {
    n <- ncol(resu)
    for(i in 1:(n - 1)){  
      nome1 <- nomes[i]
      for(j in 1:n){
        nome2 <- nomes[j]    
        if (nome1 == nome2) {
          divid <- which((resu[, i] != 0 & resu[, j] != 0))    
          resu[,i] <- resu[, i] + resu[, j]
          resu[divid, i] <- resu[divid, i]/2
        }
      }
    }
    pos <- duplicated(nomes)
    resu <- resu[, !pos, drop = FALSE]
    if (is.vector(resu)) {
      nomes <- names(resu)
      resu <- matrix(resu, ncol = length(resu))
      colnames(resu) <- nomes                
    }
    return(resu)
  }
}

