jm.conv <- function(a){ 
  v1 <- seq(1, dim(a)[2], by=2)
  v2 <- v1+1
  a2 <- as.matrix(a)
  
  a3 <- apply(a2, 2, round)
  
  alls <- list(NA)
  for(i in 1:length(v1)){
    s1 <- v1[i]
    s2 <- v2[i]
    x <- num.to.lett(as.matrix(a3[,s1:s2]))
    alls[[i]] <- apply(x,2, letter.to.jm)
  }
  res <- data.frame(matrix(unlist(alls), ncol=length(alls)))
  names(res) <- names(a)[v1]
  return(res)
}