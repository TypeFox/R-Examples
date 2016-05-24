mnntsmanifoldnewtonestimation<-function (data, M = 0, R=1, iter = 1000, initialpoint = FALSE, cinitial) 
{
    data <- as.matrix(data)
    n <- nrow(data)
    if (R != length(M)) 
        return("Error: Length of M and number of dimensions are not equal")
    
    sec<-list(R)

    for (k in 1:R){
	sec[[k]] <- 0:M[k]
    }

    ind<-expand.grid(sec,KEEP.OUT.ATTRS = FALSE)
    ind<-as.matrix(ind)

    statisticsmatrix <- matrix(0, nrow = prod(M + 1), ncol = n)
    statisticsmatrix <- exp(-(ind %*% t(data))*complex(real=0,imaginary=1))


    if (initialpoint)
    { 
        size <- length(cinitial)
    	if (size != prod(M + 1)) 
        	return("Error: Length of cinitial must be equal to prod(M+1)")
        c0 <- cinitial
    }
    else c0 <- apply(statisticsmatrix, 1, mean)
    c0 <- c0/sqrt(sum(Mod(c0)^2))


    eta <- matrix(0, nrow = prod(M + 1), ncol = 1)
    for (k in 1:n) 
	eta <- eta + (1/n) * (1/(t(Conj(c0)) %*% statisticsmatrix[, k])) * statisticsmatrix[, k]
    eta <- eta - c0
    newtonmanifold <- (c0 + eta)
    newtonmanifold <- newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
    newtonmanifold <- newtonmanifold * exp(-(0+1i) * Arg(newtonmanifold[1]))
    newtonmanifoldprevious <- newtonmanifold

    for (j in 1:iter) {
        eta <- matrix(0, nrow = prod(M + 1), ncol = 1)
        for (k in 1:n) {
            eta <- eta + (1/n) * (1/(t(Conj(newtonmanifold)) %*% 
                statisticsmatrix[, k])) * statisticsmatrix[, 
                k]
        }
        eta <- eta - newtonmanifold
        newtonmanifold <- newtonmanifold + eta
        newtonmanifold <- newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
        newtonmanifold <- newtonmanifold * exp(-(0+1i) * Arg(newtonmanifold[1]))
        if (j == iter) 
            normsequence <- (sqrt(sum(Mod(newtonmanifold - newtonmanifoldprevious)^2)))
        newtonmanifoldprevious <- newtonmanifold
    }
    newtonmanifold <- newtonmanifold/sqrt( (2 * pi)^R)
    loglik <- mnntsloglik(data, newtonmanifold, M, R)
    AIC <- -2 * loglik + 2 * (2*(prod(M+1)-1))
    BIC <- -2 * loglik + (2*(prod(M+1)-1)) * log(n)
    gradnormerror <- normsequence


    cestimatesarray <- data.frame(cbind(ind, newtonmanifold))
    cestimatesarray[, 1:R] <- as.integer(Re(as.matrix(cestimatesarray[, 
        1:R])))
    names(cestimatesarray)[1:R] <- 1:R
    names(cestimatesarray)[R+1] <- "cestimates"
    res <- list(cestimates = cestimatesarray, loglik = loglik, 
        AIC = AIC, BIC = BIC, gradnormerror = gradnormerror)
    return(res)
}

#cestimatearray es la salida de mnntsmanifoldnewtonestimation $cestimation (las tres columnas)