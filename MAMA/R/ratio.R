# calculate the ration of co-significant: expectd/observed
##
ratio <- function(X.discret)
{
  N <- nrow(X.discret)
  n.entity <- ncol(X.discret)
# filter full zero lines
  X.discret <- X.discret[which(apply(X.discret,1,sum) > 0),]
  X.string <- patternToString(X.discret)

  unique.X <- unique(X.discret)
  unique.pat <- patternToString(unique.X)

  n.pat <- length(unique.pat)
  p.soft <- array(0,n.pat)
  p.strong <- p.soft
  n <- apply(X.discret,2,sum)
  p <- n/N
  q <- 1-p
  for(i in 1:n.pat)
 {
   prob <- 1
   for (j in 1:n.entity)
   {
     if  (unique.X[i,j]==1)   prob <- prob * p[j]
  }
   p.soft[i] <- prob
 }

  for(i in 1:n.pat)
 {
   prob <- 1
   for (j in 1:n.entity)
   {
     if  (unique.X[i,j]==1)   prob <- prob * p[j] else  prob <- prob * q[j]
  }
   p.strong[i] <- prob
 }

 return(list(n=n,X.string=X.string,p.strong=p.strong,p.soft=p.soft))
}

