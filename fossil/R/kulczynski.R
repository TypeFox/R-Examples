`kulczynski` <-
function(x,y) 
{
  ai <- (y[y>0]-x[y>0])!=y[y>0]
  a <- length(ai[ai==TRUE])
  b<-length(x[x>0])-a
  c<-length(y[y>0])-a
  kul <- (a/(a + b) + a/(a + c))/2
  return(kul)
}

