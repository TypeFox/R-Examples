snntsdensity <- function (data, cpars = 1, M = c(0,0)) 
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

    if (abs(sum(Mod(cpars)^2) -  1) > 1e-10) 
        return("cpars must lie on the surface of the unit hypersphere")
    if (sum(M) == 0)
	return(sin(data[,2])/(4*pi))

    statisticsmatrix <- matrix(0, nrow = (M[1] + 1)*(M[2] + 1), ncol = n)

    for (j in 1:nrow(data)){
    	for (k1 in 0:M[1]){
		for (k2 in 0:M[2]){ 
			statisticsmatrix[k1*(M[2]+1) + k2 + 1,j] <- sqrt(sin(data[j,2])) * exp((0+1i) * (k1 * data[j,1] + k2 * data[j,2]))
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

    aux <- t(cpars)%*%statisticsmatrix
    res <- aux * Conj(aux)
    return(t(Re(res)))
}
