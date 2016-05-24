### function to compute a symmetric square root of a matrix
### internal function
###

`mat.sqrt` <-
function(A) 
    {
    eigen.A<-eigen(A)
    sqrt.A<-eigen.A$vectors %*% tcrossprod(diag(eigen.A$values^0.5), eigen.A$vectors)
    return(sqrt.A)
    }
