TraFromMpl <- function(mpl){
  tra <- matrix(nrow = length(mpl), ncol = length(mpl), FALSE)
  for(k in 1:length(mpl)){
    if(!is.null(mpl[[k]]$all.to)){
      tra[k, mpl[[k]]$all.to] <- TRUE
    }
  }
  return(tra)
}