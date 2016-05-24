"InvertQ" <- 
function(coef){
 stopifnot(class(coef)=="numeric"||class(coef)=="matrix"||(class(coef)=="array" && (dim(coef)[1]==dim(coef)[2])))
     if (class(coef) == "numeric")
       coef <- array(coef,dim=c(1,1,length(coef)))
     if (class(coef) == "matrix")
       coef <- array(coef,dim=c(NROW(coef),NROW(coef),1))
     k <- dim(coef)[1]
     order <- dim(coef)[3]
     if (order==1)
       ans <- eigen(coef[,,1], symmetric=FALSE, only.values =TRUE)$value 
     else{
        blockMat <- matrix(numeric((k*order)^2),k*order,k*order)
        blockMat[1:k,] <- coef
        Imat <- diag((order-1)*k)
        blockMat[((k+1):(k*order)),1:((order-1)*k)] <- Imat
        ans <- eigen(blockMat, symmetric=FALSE, only.values =TRUE)$value 
     }
   MaxEigenvalue <- max(Mod(ans))
   if (MaxEigenvalue >= 1) 
     return( warning("check stationary/invertibility condition !"))
}