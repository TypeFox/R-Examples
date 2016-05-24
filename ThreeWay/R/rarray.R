rarray <-
function(Xa,n,m,p){

Xa=as.matrix(Xa)
X=array(0,dim=c(n,m,p))
k=1
while (k<=p){
    j=1
    while (j<=m){
        X[,j,k]=Xa[, (k-1)*m+j]
        j=j+1
    }
    k=k+1
}
return(X)
}
