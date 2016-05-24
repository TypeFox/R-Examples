hglasso <-
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
  p <- ncol(S)
  A <- 1*(abs(S)>=min(lambda1,lambda2/2))
  diag(A) <- 1
  g1 <- graph.adjacency(A)
  cout <- clusters(g1)
  unconnected <- NULL
  blocklist <- list()
  rho <- 2.5

  for(i in 1:cout$no){
    if(sum(cout$membership==i)==1){ 
      unconnected <- c(unconnected,which(cout$membership==i))
    }
    if(sum(cout$membership==i)>1){
      blocklist[[length(blocklist)+1]]<- which(cout$membership==i)
    }
  }

# Initialize the variables   
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
    
diag(Theta)[unconnected] <- 1 
    
  if(length(blocklist)>0){
    for(i in 1:length(blocklist)){
      bl <- blocklist[[i]]
      temp <- hgl(S[bl,bl],lambda1,lambda2,lambda3,convergence=convergence,maxiter=maxiter,start=start,var.init=var.init,trace=trace)
      Theta[bl,bl] <- temp$Theta
      V[bl,bl] <- temp$V
      Z[bl,bl] <- temp$Z
      W1[bl,bl] <- temp$W1
      W2[bl,bl] <- temp$W2
      W3[bl,bl] <- temp$W3
      Gamma[bl,bl] <- temp$Gamma      
      tildeTheta[bl,bl] <- temp$tildeTheta      
      tildeV[bl,bl] <- temp$tildeV      
      tildeZ[bl,bl] <- temp$tildeZ      
    }
  }

  Obj <- HGLObj(Theta,S,Z,V,lambda1,lambda2,lambda3)
  tempV<-V
  diag(tempV) <- 0
  tempV[which(abs(tempV)<1e-5)] <- 0
  hubind <- which(apply(abs(tempV),2,sum)!=0)
# return output
  out <- list()
  class(out) <- "hglasso"
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
  out$blocklist <- blocklist
  out$unconnected <- unconnected
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
