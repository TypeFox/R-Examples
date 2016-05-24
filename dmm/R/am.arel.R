am.arel <-
function(df){
  # am.arel()
  #compute additive relationship matrix from a data frame
  m <- length(df$Id)   # no of animals in data frame
  arel <- matrix(0,m,m)
  for(i in 1:m) {
    pat <- df$SId[i]
    mat <- df$DId[i]
    # diagonal elements
    if(is.na(pat) || is.na(mat)) {
      arel[i,i] <- 1.0
    }
    else {
      arel[i,i] <- 1.0 + 0.5 * arel[pat,mat]
    }
    # off diagonal elements
    if(is.na(pat) && is.na(mat)) {
      # base animal - sire and dam unknown
      if(i > 1) {
        for(j in 1:(i-1)) {
  arel[j,i] <- 0.0
  arel[i,j] <- 0.0
        }
      }
    }
    else {
      if (i > 1) {
        for(j in 1:(i-1)) {
            if(is.na(pat)) {
              apat <- 0
            }
            else {
              apat <- arel[j,pat]
            }
            if(is.na(mat)) {
              amat <- 0
            }
            else {
              amat <- arel[j,mat]
            }
      arel[j,i] <- 0.5 * (apat + amat)
    arel[i,j] <- arel[j,i]
        }
      }
    }
  }
  return(arel)
}
