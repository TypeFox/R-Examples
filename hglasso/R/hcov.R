hcov <-
function(S,lambda1,lambda2=100000,lambda3=100000,convergence=1e-10,maxiter=1000,start="cold",var.init=NULL,trace=FALSE){

# Checking arguments
if(isSymmetric(S)==FALSE) stop("Input matrix is not symmetric!")
if(lambda1<0) stop("tuning parameter cannot be negative!")
if(lambda2<0) stop("tuning parameter cannot be negative!")
if(lambda3<0) stop("tuning parameter cannot be negative!")
if(lambda1==0 & lambda2==0 & lambda3==0) stop("At least one tuning parameter needs to be nonzero!")
if(start!="cold" && start!="warm") stop("Invalid argument")

	
cl <- match.call()
  hubind <- NULL
  n <- nrow(S)
  p <- nrow(S)
  rho <- 2.5
  i <- 1
# Variables initialization
  if(start=="cold"){
  oldSigma <- Sigma <- V <- Z <- diag(rep(1,p))
  W1 <- W2 <- W3 <- Gamma <- matrix(0,p,p)
  tildeSigma <- tildeV <- tildeZ <- matrix(0,p,p)
  }
  
  else if(start=="warm"){
  Sigma <- diag(rep(1,p))	
  oldSigma <- var.init$Sigma
  V <- var.init$V
  Z <- var.init$Z  
  W1 <- var.init$W1 
  W2 <- var.init$W2
  W3 <- var.init$W3
  Gamma <- var.init$Gamma  
  tildeSigma <- var.init$tildeSigma    
  tildeV <- var.init$tildeV
  tildeZ <- var.init$tildeZ
  }
  
  criteria <- 1e10
  	
# While loop for the iterations
  while(criteria > convergence && i<=maxiter){
	if(trace==TRUE && i%%10==0) print(paste("iteration = ",i,sep=""))

    Sigma <- updateSigma(tildeSigma,W1,S,rho,0.001)

    Z <- updateZ(tildeZ,W3,lambda1,rho)	 

    V <- updateV(tildeV,W2,lambda2,lambda3,rho)

    Gamma <- updateGamma(Sigma,V,Z,W1,W2,W3,rho)

    tildeSigma <- updatetildeTheta(Sigma,W1,Gamma,rho)

    tildeV <- updatetildeV(V,W2,Gamma,rho)

    tildeZ <- updatetildeZ(Z,W3,Gamma,rho)
    
    W1 <- W1+Sigma-tildeSigma
    W2 <- W2+V-tildeV
    W3 <- W3+Z-tildeZ

	i <- i+1
	Sigma <- Z+V+t(V)	
	  
	criteria <- sum((Sigma-oldSigma)^2)/sum((oldSigma)^2)
	oldSigma <- Sigma

  } 	

  if(i>maxiter){
  	 warning("The algorithm has not converged by the specified maximum number of iteration")
  }

  Sigma<-Z+V+t(V)
  Obj <- HcovObj(Sigma,S,Z,V,lambda1,lambda2,lambda3)

	tempV<-V
	diag(tempV) <- 0
    tempV[which(abs(tempV)<1e-5)] <- 0
  	hubind <- which(apply(abs(tempV),2,sum)!=0)
# return output
	out <- list()
	class(out) <- "hcov"
	out$Sigma <- Sigma
	out$V <- V
	out$Z <- Z
	out$W1 <- W1
	out$W2 <- W2  
	out$W3 <- W3  
	out$Gamma <- Gamma  
	out$tildeTheta <- tildeSigma
	out$tildeV <- tildeV
	out$tildeZ <- tildeZ    
	out$objective <- Obj
	out$cl <- cl  
	out$hubind <- hubind
	out$lambda1 <- lambda1
	out$lambda2 <- lambda2
	out$lambda3 <- lambda3	
  	out$p <- p
  	out$n <- n	
	return(out)
	
}
