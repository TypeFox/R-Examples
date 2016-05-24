cl.memb <-
function(U)
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
n=nrow(U)
info.U=cbind(max.col(U),apply(U,1,max))
if (is.null(rownames(U)))
rownames(info.U)=paste("Obj",1:n,sep=" ")
else
rownames(info.U)=rownames(U)
colnames(info.U)=c("Cluster","Membership degree")
return(info.U)
}
