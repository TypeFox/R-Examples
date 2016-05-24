estmf0 <-
function(x,p,alphau)
{
x <- as.matrix(x)
n <- nrow(x)
y <- x[(p+1):n,1]- alphau*x[p:(n-1),1]
{
if( alphau == 1){
    xmat <- numeric()
    if(p > 1){
    z<- x[2:n,1]-x[(1:n-1),1]
    index <- (p-1):(n-2)
    for(i in 1:(p-1)){
    xmat <- cbind(xmat,z[index])
    index <- index -1}
    }
    btem <- rbind( solve( t(xmat) %*% xmat) %*% t(xmat) %*% y , 0) 
    }
else{
    xmat <- matrix(1,nrow=n-p)
    if(p > 1){
    z<- x[2:n,1]-x[(1:n-1),1]
    index <- (p-1):(n-2)
    for(i in 1:(p-1)){
    xmat <- cbind(xmat,z[index])
    index <- index -1}
    xmat <- cbind(xmat[,2:ncol(xmat)],xmat[,1])
    }
    btem <-  solve( t(xmat) %*% xmat) %*% t(xmat) %*% y     
    }
}
return(btem)
}
