snntsmanifoldnewtonestimation <- function (data, M = c(0,0), iter = 1000, initialpoint = FALSE, cinitial) 
{
    data <- as.matrix(data)
    n <- nrow(data)
    R <- 2
    if (R != length(M)) 
        return("Error: Dimensions of M and vector of observations are not equal")
    auxcond1 <- sum(data[,1] > 2*pi) + sum(data[,1] < 0)
    auxcond2 <- sum(data[,2] > pi) + sum(data[,2] < 0)
    if (auxcond1>0)
	return("First column of the data matrix must have values between 0 and 2*pi")
    if (auxcond2>0)
	return("Second column of the data matrix must have values between 0 and pi")

    statisticsmatrix <- matrix(0, nrow = (M[1] + 1)*(M[2] + 1), ncol = n)

    for (j in 1:nrow(data)){
    	for (k1 in 0:M[1]){
		for (k2 in 0:M[2]){ 
			statisticsmatrix[k1*(M[2]+1) + k2 + 1,j] <- sqrt(sin(data[j,2])) * Conj(exp((0+1i) * (k1 * data[j,1] + k2 * data[j,2])))
		}
	}
    }

    A <- matrix(0,nrow=M[2]+1,ncol=M[2]+1)

    for (k2 in 0:M[2]){
    	for (m2 in 0:M[2]){
        	if (abs(k2 - m2) != 1){
            		A[k2+1,m2+1] <- (2*pi)*((1 + cos((k2-m2)*pi))/(1 - ((k2-m2)^2)));
 		}
	}
    }

    Ac<-chol(A)

    Acinv<-solve(Ac)	

    for (j in 1:nrow(data)){
    	for (k1 in 0:M[1]){
		statisticsmatrix[(k1*(M[2]+1)+1):((k1+1)*(M[2]+1)),j] <- t(Acinv) %*% statisticsmatrix[(k1*(M[2]+1)+1):((k1+1)*(M[2]+1)),j]
	}
    }

    if (initialpoint) 
        c0 <- cinitial
    else c0 <- apply(statisticsmatrix, 1, mean)
    c0 <- c0/sqrt(sum(Mod(c0)^2))

    c0<-Conj(c0)

    eta <- matrix(0, nrow = (M[1] + 1)*(M[2] + 1), ncol = 1)

    for (k in 1:n) 
    	eta <- eta + (1/n) * (1/(t(Conj(c0)) %*% statisticsmatrix[, k])) * statisticsmatrix[, k]
    eta <- eta - c0
    newtonmanifold <- (c0 + eta)
    newtonmanifold <- newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
    newtonmanifold <- newtonmanifold * exp(-(0+1i) * Arg(newtonmanifold[1]))
    newtonmanifoldprevious <- newtonmanifold
    for (j in 1:iter) {
        eta <- matrix(0, nrow = (M[1] + 1)*(M[2] + 1), ncol = 1)
        for (k in 1:n) {
            eta <- eta + (1/n) * (1/(t(Conj(newtonmanifold)) %*% 
                statisticsmatrix[, k])) * statisticsmatrix[, k]
        }
        eta <- eta - newtonmanifold
        newtonmanifold <- newtonmanifold + eta
        newtonmanifold <- newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
        newtonmanifold <- newtonmanifold * exp(-(0+1i) * Arg(newtonmanifold[1]))
        if (j == iter) 
            normsequence <- (sqrt(sum(Mod(newtonmanifold - newtonmanifoldprevious)^2)))
        newtonmanifoldprevious <- newtonmanifold
    }

    loglik <- snntsloglik(data, newtonmanifold, M)
    AIC <- -2 * loglik + 2 * (2 * (M[1]+1)*(M[2]+1) - 2)
    BIC <- -2 * loglik + (2 * (M[1]+1)*(M[2]+1) - 2) * log(n)
    gradnormerror <- normsequence
    
    index<-expand.grid(0:M[1],0:M[2])
    index<-index[order(index[,1]),]

    cestimatesarray <- data.frame(cbind(index, newtonmanifold))
    cestimatesarray[, 1:2] <- as.integer(Re(as.matrix(cestimatesarray[,1:2])))
    names(cestimatesarray)[1] <- "k1"
    names(cestimatesarray)[2] <- "k2"
    names(cestimatesarray)[3] <- "cestimates"

    res <- list(cestimates = cestimatesarray, loglik = loglik, 
        AIC = AIC, BIC = BIC, gradnormerror = gradnormerror)
    return(res)
}