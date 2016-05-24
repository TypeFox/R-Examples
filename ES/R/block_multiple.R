block_multiple <-
function(A, B, ii, sortcol, sortcol2)
{
result <- NULL
object <- NULL
for (i in 1:ncol(A))
{
if(sum(sortcol[ii]==i)==0)
object <- rep(0,nrow(A))
if(sum(sortcol[ii]==i)>0)
object <- A[,sortcol2[ii][which(sortcol[ii]==i)]]%*%as.matrix(B[which(sortcol[ii]==i),])
result <- c(result,object)
}
return(as.vector(result))
}
