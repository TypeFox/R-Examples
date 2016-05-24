permnew <-
function(X,n,m,p){

X=as.matrix(X)
Y=t(X)		
y=as.vector(Y)		
Y=matrix(0,m,n*p)
Y=matrix(y,nrow(Y),ncol(Y))
return(Y)
}
