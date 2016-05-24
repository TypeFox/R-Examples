#Funcion auxiliar para la funcion de grafica bivariada 


mnntsdensity2d<-function (x,y, cpars,M) 
{
    R<-2
    data<-matrix(c(x,y),ncol=2)
    data <- as.matrix(data)
    n <- nrow(data)



#    if (R != length(M)) 
#        return("Error: Dimensions of M and vector of observations are not equal")
#    if (abs(sum(Mod(cpars)^2) - ((1/(2 * pi))^R)) > 1e-10) 
#        return("sum of the squared norm of componentes greater than condition")
#    if (sum(M) == 0) 
#        return(rep((1/(2 * pi))^R, n))


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