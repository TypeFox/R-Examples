patternMatch <- function(X.discret,unique.pat)
{
  N <- nrow(X.discret)
  n.entity <- ncol(X.discret)
# filter full zero lines
  X.discret <- X.discret[which(apply(X.discret,1,sum) > 0),]
  X.string <- patternToString(X.discret)

  n.random <- array(0,length(unique.pat))

  unique.X <- patternmatrix(unique.pat,n.entity)
  for ( i in 1:length(unique.pat))
  {
    comp.called <- which(unique.X[i,] == 1)
    subX <- as.matrix(X.discret[,comp.called])
    n.random[i] <- length(which(apply(subX,1,sum)==length(comp.called)))
  }

 return(n.random)

}
