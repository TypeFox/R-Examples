transfmatcat <-
function(mat,cat=3)
{
n<-dim(mat)[1]
m<-dim(mat)[2]
cat<-rep(cat,length.out=m)
for (j in 1:m)
{
mat[,j]<-as.integer(cut(mat[,j],breaks=c(min(mat[,j])-1, qnorm(seq(1/cat[j],1,1/cat[j]),mean(mat[,j]),sd(mat[,j])),max(mat[,j])+1)))
}
mat
}

