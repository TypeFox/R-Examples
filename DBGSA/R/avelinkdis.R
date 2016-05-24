#用类平均法求:
avelinkdis<-function(C,num,Meth){
#例如：Meth为"euclidean"
C<-t(C);
A<-as.matrix(C[1:num,]);
B<-as.matrix(C[(1+num):nrow(C),]);
rowA<-nrow(A);
rowB<-nrow(B);
if (Meth=="stat")
{
cr<-nrow(C);
cc<-ncol(C);
for(n in 1:cc){
C[,n]<-C[,n]*(1/sd(C[,n]))
}
D<-dist(C,method="euclidean",diag=TRUE,upper=TRUE)
}
else
{
D<-dist(C,method=Meth,diag=TRUE,upper=TRUE);##算距离阵的下三角阵并且显示对角线上的元素
}
D<-as.matrix(D);#强制将D的dist数据类型转换为矩阵格式
Out<-D[(rowA+1):nrow(D),1:rowA];
return(1/(rowA*rowB)*sum(Out));#根据公式求出的距离
 }
