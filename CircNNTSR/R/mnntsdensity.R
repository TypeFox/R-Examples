mnntsdensity<-function (data, cpars = 1/sqrt(2 * pi), M = 0, R=1) 
{
    data <- as.matrix(data)
    n <- nrow(data)



    if (R != length(M)) 
        return("Error: Length of M and number of dimensions are not equal")
    size <- length(cpars)
    if (size != prod(M + 1)) 
        return("Error: Length of cpars must be equal to prod(M+1)")
    if (abs(sum(Mod(cpars)^2) - ((1/(2 * pi))^R)) > 1e-10) 
        return("Error: Sum of the squared norm of componentes greater than condition")
    if (sum(M) == 0) 
        return(rep((1/(2 * pi))^R, n))


    sec<-list(R)

    for (k in 1:R){
	sec[[k]] <- 0:M[k]
    }

    ind<-expand.grid(sec,KEEP.OUT.ATTRS = FALSE)
    ind<-as.matrix(ind)

    statisticsmatrix <- matrix(0, nrow = prod(M + 1), ncol = n)
    statisticsmatrix <- exp((ind %*% t(data))*complex(real=0,imaginary=1))
    aux <- t(as.matrix(cpars)) %*% statisticsmatrix

    res <- Re(aux*Conj(aux))
    res <- c(t(res))
    return(res) 
}