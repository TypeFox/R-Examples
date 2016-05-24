medianimp <-
function(mat)
{
m<-ncol(mat)
for (i in 1:m)
{
median<-quantile(mat[which(is.na(mat[,i])==FALSE),i], prob=0.5, type=1)[[1]]
mat[which(is.na(mat)[,i]),i]<-median
}
mat
}

