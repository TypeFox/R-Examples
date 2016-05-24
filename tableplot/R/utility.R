## last modified 1-22-10 MF: fixed bug in transpose

# Transpose an array of 2 or 3 dimensions,
# where transposition is carried out on the first two dimensions:
transpose <- function(x) {
	# TODO: check for class matrix or array or something else suitable
	if (length(dim(x))==2) result = t(x)
	else if (length(dim(x))==3) result = aperm(x, c(2,1,3))
	else stop("Can't transpose")
	result
}

## Function to calculate a congruence coefficient
congruence.coef <- function(a,b){
	(t(a)%*%b) / sqrt( (t(a)%*%a) * (t(b) %*% b) )
	}
cg.cf <- congruence.coef

identity.coef <- function(a,b){
	2 * (t(a)%*% b) / ( (t(a)%*%a) + (t(b)%*%b) )
	}
id.cf <- identity.coef

## Function to facilitate the comparison of two or more factor patterns
## stored as a 3-dimensional array; if only two patterns, a matrix of 
## differences are calculated; if three of more patterns, a matrix of
## standard deviations are calculated

compare <- function(X){
	
	d1 <- dim(X)[1]
	d2 <- dim(X)[2]
	d3 <- dim(X)[3]
	
	if (d3 == 2) {
		
		Y <- array(NA, c(d1+1, (d2*2)+3, d3))
		Y[1:d1,1:d2,1:d3] <- X
		
		Y[1:d1, (d2+3):(dim(Y)[2]-1), 1] <- X[,,1]-X[,,2]
		for (i in 1:d1) { Y[i,d2+1,1] <- congruence.coef(X[i,,1], X[i,,2])*100 }
		for (j in 1:d2) { Y[d1+1,j,1] <- congruence.coef(X[,j,1], X[,j,2])*100 }
		Y[d1+1,d2+1,1] <- congruence.coef(as.vector(X[,,1]),as.vector(X[,,2]))*100
		
		for(i in 1:d1) { Y[i,dim(Y)[2],1] <- mean(abs(Y[i,(d2+3):(dim(Y)[2]-1),1]))}
		for(j in (d2+3):(dim(Y)[2]-1)) { Y[dim(Y)[1],j,1] <- mean(abs(Y[1:d1,j,1]))}
		Y[d1+1,(d2*2)+3,1] <- mean(abs(Y[1:d1,(d2+3):((d2*2)+2),1]))
		}
	
	if (d3 >= 3) {
		
		Y <- array(NA, c(d1+1, (d2*2)+2, d3))
		Y[1:d1,1:d2,1:d3] <- X
		
		for (i in 1:d1){
			for (j in 1:d2){
				Y[i,j+d2+1,1] <- sd(X[i,j,])
				}
			}
		
		for (i in 1:d1) { Y[i,dim(Y)[2],1] <- mean(Y[i,(d2+2):(d2+1+d2),1]) }
		for (j in (d2+2):(dim(Y)[2]-1)) { Y[dim(Y)[1],j,1] <- mean(Y[1:d1,j,1])}
		Y[d1+1,(d2*2)+2,1] <- mean(abs(Y[1:d1,(d2+2):((d2*2)+1),1]))
		}
		
	round(Y)
	
	}