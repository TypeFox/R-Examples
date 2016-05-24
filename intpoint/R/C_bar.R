C_bar <-
function(c,A,Dx){
At<-t(A)
if(nrow(At)<=40){
Dx2<-Dx^2
ec<-A%*%Dx2%*%At
aux <- svd(ec)
Positive <- aux$d > max(1e-08 * aux$d[1L], 0)
if (all(Positive)) 
inv<-aux$v %*% (1/aux$d * t(aux$u))
else if (!any(Positive)) 
inv<-array(0, dim(ec)[2L:1L])
else inv<-aux$v[, Positive, drop = FALSE] %*% ((1/aux$d[Positive]) * 
t(aux$u[, Positive, drop = FALSE]))
C<-c-At%*%inv%*%A%*%Dx2%*%c
return(C)
}
else {
Dx2<-Dx^2
B<-A%*%Dx2%*%At
y<-A%*%Dx2%*%c
w<-solve(B,y)
C<-c-At%*%w
return(C)}
}
