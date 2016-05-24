duplicMat <- function(m, n){
  dm <-NULL
  for(i in 1:n){
   dm <- cbind(dm, m)
  }
  dm
}

duplicSumDirect <- function(m, n){
  m <- as.matrix(m)
  dm <-m
  if(n> 1){
    for(i in 2:n){
      dm <- dm %s% m
    }
  }
  dm
}


StackDiag <- function(vect, n){
  dm <-NULL
  for(i in 1:length(vect)){
   dm <- rbind(dm, diag(vect[i],n))
  }
  dm
}

