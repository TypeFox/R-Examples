"nipals" <-
function(X,a,it=10,tol=1e-4) 
#fct nipals calculates the principal components of a given data matrix X according to 
#the NIPALS algorithm (Wold).
#X...data matrix, a...number of components, 
#it...maximal number of iterations per component,
#tol...precision tolerance for calculation of components
{
Xh <- scale(X,center=TRUE,scale=FALSE)		#mean-centering of data matrix X
nr <- 0
T <- NULL
P <- NULL
for (h in 1:a){
	th <- Xh[,1]		#starting value for th is 1st column of Xh
	ende <- FALSE
	#3 inner steps of NIPALS algorithm
	while (!ende){
		nr <- nr+1
		ph <- t((t(th)%*%Xh) * as.vector(1/(t(th)%*%th)))	#LS regression for ph
		ph <- ph * as.vector(1/sqrt(t(ph)%*%ph))		#normalization of ph
		thnew <- t(t(ph)%*%t(Xh) * as.vector(1/(t(ph)%*%ph)))	#LS regression for th
		prec <- t(th-thnew)%*%(th-thnew)	#calculate precision
		cat("actual precision: ",sqrt(prec),"\n")
		th <- thnew	#refresh th in any case
		#check convergence of th
		if (prec <= (tol^2)) {
			ende <- TRUE
		}
		else if (it <= nr) {	#too many iterations
			ende <- TRUE
			cat("\nWARNING! Iteration stop in h=",h," without convergence!\n\n")
		}
	}
	Xh <- Xh-(th%*%t(ph))	#calculate new Xh
	T <- cbind(T,th)	#build matrix T
	P <- cbind(P,ph)	#build matrix P
	nr <- 0
}
list(T=T,P=P)
}

