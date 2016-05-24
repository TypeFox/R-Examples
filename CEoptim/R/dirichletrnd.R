dirichletrnd <- function(a,n) { 
  out <- matrix(nrow=n,ncol=length(a))  
  for (i in 1:n){
    y <- rgamma(length(a), a, 1)    
    out[i,] <- y / sum(y)
  }
  return(out)
}