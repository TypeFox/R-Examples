meanimp <-
function(mat)
{
m<-ncol(mat)
for (i in 1:m)
{
mat[which(is.na(mat)[,i]),i]<-mean(mat[which(is.na(mat[,i])==FALSE),i])
}
mat
}

