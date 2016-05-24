modeimp <-
function(mat)
{
m<-ncol(mat)
for (i in 1:m)
{
mode<-which.max(tabulate(mat[which(is.na(mat[,i])==FALSE),i]))
mat[which(is.na(mat)[,i]),i]<-mode
}
mat
}

