cl.memb.t <-
function(U,t)
{
if (missing(U))
stop("The membership degree matrix U must be given")
if (is.null(U))
stop("The membership degree matrix U is empty")
U=as.matrix(U)
if (any(is.na(U)))
stop("The membership degree matrix U must not contain NA values")
if (!is.numeric(U)) 
stop("The membership degree matrix U must be numeric")
if (missing(t))
t=0
if (!is.numeric(t))
{
t=0
cat("The threshold t is not numeric: the default value t=0 will be used ",fill=TRUE)
}
if ((t<0) || (t>1))
{
t=0
cat("The threshold t must be in the interval [0,1]: the default value t=0 will be used ",fill=TRUE)
}
n=nrow(U)
info.U=cbind(max.col(U),apply(U,1,max))
for (i in 1:n)
{
if (info.U[i,2]<t)
{
info.U[i,1]=0
info.U[i,2]=NA
}
}
if (is.null(rownames(U)))
rownames(info.U)=paste("Obj",1:n,sep=" ")
else
rownames(info.U)=rownames(U)
colnames(info.U)=c("Cluster","Membership degree")
return(info.U)
}