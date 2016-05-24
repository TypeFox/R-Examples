mnntsrandominitial<- function(M=1,R=1) 
{
    if (sum(M)==0)
        return(1/(2*pi)^R)
    if (R != length(M)) 
        return("Error: Length of M and number of dimensions are not equal")
    res <- complex(prod(M + 1))
    aux <- rnorm(2*prod(M + 1))
    aux <- sqrt((1/(2 * pi))^R) * (aux/sqrt(sum(aux^2)))
#    aux <- (aux/sqrt(sum(aux^2)))
    res[1]<-Mod(aux[1] + 1i*aux[prod(M+1)+1])
    for (k in 2:prod(M+1)){
	res[k] <- aux[k] + 1i*aux[prod(M+1)+k]
    }
    return(res)
}