`lambda.stem.ms01` <-
function(n, tb, eps=0)
{
  res<-list()
  x<- (1/tb)*log(n*(1-eps)+eps)
  res$r <- x
  res$lambda <- x/(1-eps) 
  res
}

