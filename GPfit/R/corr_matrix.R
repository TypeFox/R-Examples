## Creating the correlation matrix based on hyperparameters 
## beta and a data matrix x
##
## May 8th, 2012

corr_matrix <- function(X,beta,corr=list(type="exponential",power=1.95)){
## Checking to make sure the data is a matrix, and sets it as one 
## if it is not
if (is.matrix(X) == FALSE){
	X = as.matrix(X)
}
d = ncol(X)
n = nrow(X)
## Checking the dimensions between the two inputs 
if (d != length(beta)){
	stop("The dimensions of beta and X do not match. \n")
}
R = matrix(0,n,n)
Points = expand.grid((1:n),(1:n))
Points = cbind(Points[,2],Points[,1])
#Points = Points[Points[,2]>Points[,1],]
Points = matrix(Points[Points[, 2] > Points[, 1], ],,2)
junk = abs(X[Points[,2],] - X[Points[,1],])
if (corr$type == "exponential"){
	power = corr$power
	if (is.null(power)){
		power = 1.95
	}
	dim(beta) = c(d,1)
	junk = junk^power
	Beta = matrix(beta,nrow = (length(junk)/d),ncol=d,byrow=TRUE)
	Rtemp=(10^Beta)*junk
	Rtemp=rowSums(Rtemp)
	R[Points]=Rtemp
	R = R+t(R)
	R=exp(-R)
}
if (corr$type == "matern"){
	nu = corr$nu
	if (is.null(nu)){
		nu = 2.5
	}
	dim(beta) = c(1,d)
	temp = 10^beta
	temp = matrix(temp,ncol=d,nrow=(length(junk)/d),byrow=TRUE)
	junk = 2*sqrt(nu)*junk*(temp)
	ID = which(junk==0)
	Rtemp = 1/(gamma(nu)*2^(nu-1))*(junk)^nu*besselK(junk,nu)
	Rtemp[ID]=1;

	Rtemp = apply(Rtemp,1,prod)
	R[Points] = Rtemp
	R = R + t(R)
	diag(R) = 1
}
return(R)
}

