bernstein <- function(v,x,n) {
  return((choose(n,v)*x^v*(1-x)^(n-v))*(n+1))
}

int.bernstein <- function(x,n) {
  base.int.h <- apply(matrix(0:n),1,bernstein,x,n=n)

  dim.base <- dim(base.int.h)

  base.int <- matrix(0,nrow=dim.base[1],ncol=(dim.base[2]-1))

  for(j in 1:(dim.base[2]-1)) {
    for(i in 1:dim.base[1]) {
      base.int[i,j] <- sum(base.int.h[i,((j+1):dim.base[2])])/(n+1)
    }
  }
  return(base.int)
}
