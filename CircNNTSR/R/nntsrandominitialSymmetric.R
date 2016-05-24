nntsrandominitialSymmetric <-
function(M){

res <- rep(0,2*M)
aux <- rnorm(M+1)
aux <- sqrt(1/(2*pi))*(aux/sqrt(sum(aux^2))) 
res[1:M] <- aux[1:M]^2
res[(M+1):(2*M)] <- runif(M,0,2*pi-0.00000000000001)
return(res)
}

