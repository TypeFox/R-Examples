`lambda.crown.ms01` <-
function(n, tb, eps=0) #n <- extant lineages; tb<-time of basal divergence
{
  res<-list()
  x1<- .5*n*(1-eps^2) + 2*eps
  x2 <- .5*(1-eps)*sqrt(n*((n*eps^2)- 8*eps + 2*n*eps + n))
  x <- (log(x1 + x2)-log(2))/tb
  
  res$r <- x
  res$lambda <- x/(1-eps) 
  res
}

