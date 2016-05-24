summary.fclust <-
function (object,...)  
{
fclust.obj=object  
if ((missing(fclust.obj)) || (!inherits(fclust.obj, "fclust")))
stop("An object of class fclust must be given")
X=fclust.obj$X
U=fclust.obj$U 
k=fclust.obj$k
H=fclust.obj$H
n=nrow(X)
info.U=cl.memb(U)
cat("\n Fuzzy clustering object of class 'fclust' ")
cat("\n ")
cat("\n Number of objects: \n", n); 
cat("\n ")
cat("\n Number of clusters: \n", k); 
cat("\n ")
cat("\n Cluster sizes: \n"); 
print(cl.size(U))
cat("\n Closest hard clustering partition: \n")
print(info.U[,1])
cat("\n Cluster memberships:")
cat("\n ")
for (c in 1:k)
{  
cat(" Clus", c)
objs=rownames(info.U[info.U[,1]==c,])
if (length(objs)<51)
{
cat("\n ")
print(objs)
}
else
{
cat(" (First 50 objects) \n")
print(objs[1:50])
}
}
noua=sum(info.U[,2]<0.5)
cat("\n Number of objects with unclear assignment (maximal membership degree <0.5): \n", noua)
cat("\n ")
if (noua>0)
{
cat("\n Objects with unclear assignment: \n")
print(rownames(info.U[info.U[,2]<0.5,]))
cat("\n Cluster sizes (without unclear assignments): \n") 
print(cl.size.H(U))
}
cat("\n Membership degree matrix (rounded): \n")
print(round(U,2))
cat("\n Cluster summary: \n")
minU=rep(0,k);
names(minU)=names(cl.size(U))
unasU=minU
maxU=minU
meanU=minU
for (c in 1:k)
{  
minU[c]=round(min(info.U[info.U[,1]==c,2]),2)
maxU[c]=round(max(info.U[info.U[,1]==c,2]),2)
meanU[c]=round(mean(info.U[info.U[,1]==c,2]),2)
unasU[c]=sum(info.U[info.U[,1]==c,2]<0.5)
}
summU=cbind(cl.size(U),minU,maxU,meanU,unasU)
colnames(summU)=c("Cl.size", "Min.memb.deg.", "Max.memb.deg.", "Av.memb.deg.", "N.uncl.assignm.")
print(summU)
cat("\n Euclidean distance matrix for the prototypes (rounded): \n")
print(round(dist(H),2))
cat("\n Available components: \n")
print(names(fclust.obj))
cat("\n ")
invisible(fclust.obj)
}