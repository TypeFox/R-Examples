hbn <- function(X,lambda1,lambda2=100000,lambda3=100000,convergence=1e-8,maxiter=1000,start="cold",var.init=NULL,trace=FALSE){

# Checking arguments	
if(lambda1<0) stop("tuning parameter cannot be negative!")
if(lambda2<0) stop("tuning parameter cannot be negative!")
if(lambda3<0) stop("tuning parameter cannot be negative!")
if(lambda1==0 & lambda2==0 & lambda3==0) stop("At least one tuning parameter needs to be nonzero!")
if(start!="cold" && start!="warm") stop("Invalid argument")

  cl <- match.call()
  hubind <- NULL
  S <- t(X)%*%X
  n <- nrow(X)
  p <- ncol(X)
  i <- 1
  rho <- 5
  
# Variables initialization
  if(start=="cold"){
  oldTheta <- Theta <- V <- Z <- diag(rep(1,p))
  W1 <- W2 <- W3 <- Gamma <- matrix(0,p,p)
  tildeTheta <- tildeV <- tildeZ <- matrix(0,p,p)
  }
  
  else if(start=="warm"){
  Theta <- diag(rep(1,p))	
  oldTheta <- var.init$Theta
  V <- var.init$V
  Z <- var.init$Z  
  W1 <- var.init$W1 
  W2 <- var.init$W2
  W3 <- var.init$W3
  Gamma <- var.init$Gamma  
  tildeTheta <- var.init$tildeTheta    
  tildeV <- var.init$tildeV
  tildeZ <- var.init$tildeZ
  }

  criteria <- 1e10
  while(criteria > convergence && i <= maxiter){
	if(trace==TRUE && i%%10==0) print(paste("iteration = ",i,sep=""))
  	
    Theta <- BB_logistic(X,tildeTheta-W1,rho)    

    Z<-updateZ(tildeZ,W3,lambda1,rho)	 

    V<-updateV(tildeV,W2,lambda2,lambda3,rho)

    Gamma<-updateGamma(Theta,V,Z,W1,W2,W3,rho)

    tildeTheta<-updatetildeTheta(Theta,W1,Gamma,rho)

    tildeV<-updatetildeV(V,W2,Gamma,rho)

    tildeZ<-updatetildeZ(Z,W3,Gamma,rho)
    
    W1<-W1+Theta-tildeTheta
    W2<-W2+V-tildeV
    W3<-W3+Z-tildeZ

	i <- i+1
	Theta <- Z+V+t(V)	
	  
	criteria <- sum((Theta-oldTheta)^2)/sum((oldTheta)^2)
	oldTheta<-Theta

  } 	

  if(i>maxiter){
  	 warning("The algorithm has not converged by the specified maximum number of iteration")
  }

  Theta<-Z+V+t(V)


	tempV<-V
	diag(tempV) <- 0
    tempV[which(abs(tempV)<1e-5)] <- 0
  	hubind <- which(apply(abs(tempV),2,sum)!=0)

# return output
	out <- list()
	class(out) <- "hbn"
	out$Theta <- Theta
	out$V <- V
	out$Z <- Z
	out$W1 <- W1
	out$W2 <- W2  
	out$W3 <- W3  
	out$Gamma <- Gamma  
	out$tildeTheta <- tildeTheta
	out$tildeV <- tildeV
	out$tildeZ <- tildeZ    
	out$cl <- cl  
	out$hubind <- hubind
	out$lambda1 <- lambda1
	out$lambda2 <- lambda2
	out$lambda3 <- lambda3	
  	out$p <- p
  	out$n <- n	
	return(out)
		
}