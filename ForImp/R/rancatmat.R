rancatmat <-
function(n,m,cat=3)
{
mat<-matrix(0,n,m)
if(length(cat)<m)
cat<-rep(cat,length.out=m)
for (j in 1:m)
{
mat[,j]<-floor(runif(n)*cat[j])+1
}
mat
}

