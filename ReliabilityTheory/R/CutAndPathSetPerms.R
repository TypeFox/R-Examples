cutAndPathSetPerms <- function(stsep, n) {
  p <- permn(2:(n+1))
  #p <- permn(sort(unique(unlist(stsep))))
  l <- length(stsep)
  
  permute <- function(perm) {
    for(i in 1:l) {
      for(j in 1:length(stsep[[i]])) {
        stsep[[i]][j] <- perm[stsep[[i]][j]-1]
      }
      stsep[[i]] <- sort(stsep[[i]])
    }
    stsep
  }
  stsepPerms <- lapply(p, permute)
  stsepPerms
}
