#重心法求距离：
centdis<-function(C,num,Meth){
C<-t(C);
A<-as.matrix(C[1:num,]);
B<-as.matrix(C[(1+num):nrow(C),]);
if (Meth=="stat")
{
C<-rbind(A,B);
cr<-nrow(C);
cc<-ncol(C);
for(n in 1:cc){
A[,n]<-A[,n]*(1/sd(C[,n]));
B[,n]<-B[,n]*(1/sd(C[,n]));
}
AA<-apply(A,2,mean);#矩阵A的重心（均值为AA）
BB<-apply(B,2,mean);#矩阵B的重心（均值为BB）
C<-rbind(AA,BB)
D<-as.matrix(dist(C,method="euclidean",diag=TRUE,upper=TRUE))
}
else
{AA<-apply(A,2,mean);#矩阵A的重心（均值为AA）
BB<-apply(B,2,mean);#矩阵B的重心（均值为BB）
C<-rbind(AA,BB)
D<-as.matrix(dist(C,method=Meth,diag=TRUE,upper=TRUE));}##算距离阵的下三角阵并且显示对角线上的元素
return (D[1,2])
}
