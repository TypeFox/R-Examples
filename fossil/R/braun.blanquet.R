`braun.blanquet` <-
function(x,y) 
{
ai<-(y[y>0]-x[y>0])!=y[y>0]
a<-length(ai[ai==TRUE])
if (a!= 0) {
b <- length(x[x > 0])
c <- length(y[y > 0])
simp <- a/max(b, c)
}
else simp <- 0
return(simp)
}

