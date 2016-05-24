#' For given frequency domain operator, compute its eigenvalues and eigendirections
#' such that for close frequencies eigendirection matrices are close to each other
#' ('quasi-continuity').
#'
#' @title Eigendevompose a frequency domain operator at each frequency
#' @param S frequency domain operator
#' @return Rotated matrix \code{M}
#' @importFrom graphics par plot title
#' @importFrom stats optim rnorm
#' @export
freqdom.eigen = function(S){
  if (!is.freqdom(S))
    stop("S must be a freqdom object")
  
  op = S$operators[1,,]
  if (dim(op)[1] != dim(op)[2])
    stop("S$operators[,,theta] must be a square matrix")
  
  thetas = S$freq
  E = list()
  nbasis = dim(S$operators)[2]
  T = length(thetas)
  
  E$freq = S$freq
  E$vectors = array(0,c(T,nbasis,nbasis))
  E$values = array(0,c(T,nbasis))
  Prev = diag(nbasis)
  
  for (theta in 1:T)
  {
    Eg = close.eigen(S$operators[theta,,],Prev)
    Prev = Eg$vectors
    
    ## TAKE EIGEN WITHOUT(!) HEURISTIC
    #Eg = eigen(S$operators[theta,,])
    
    E$vectors[theta,,] = Eg$vectors
    E$values[theta,] = Eg$values
  }
  
  E
#  OneSideFilter(E,10,thetas)
}


OneSideFilter = function(E,L,thetas,side="left"){
  D = dim(E$vectors)
  T = length(thetas)
  
  if (side == "left") side = 1
  else side = -1
  
  weights = c(sqrt(1:(L)),-sqrt((L+1):1))
  for (k in 1:1){
    
    fr = function(param){
      RES = t(exp(1i*param) * (E$vectors[,k,])) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) #inv fourier
      -sum(weights*abs(t(RES))^2)
#      - sum(abs(sqrt((L+1):1)*t(RES[,L + 1:(L+1)]))^2) - sum(abs(sqrt(1:L)*t(RES[,1:L]))^2)
      
    }

start = rep(0,T)
start = rnorm(T)
Opt = optim(start,fr,method = "L-BFGS-B",control = list(maxit = 100, temp = 1000, trace=TRUE))
#Opt = optim(rep(0,T),fr,method = "SANN",control = list(maxit = 30000, temp = 1000, trace=TRUE),
    #  lower=-pi,upper=pi)
    
    param = Opt$par
    # plot(exp(1i*param)*1:(2*T+1),t='l')
    ORG = t(E$vectors[,k,]) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) / T
    RES = t(exp(1i*param) * (E$vectors[,k,])) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) / T
    
    E$vectors[,k,] = (exp(1i*param) * (E$vectors[,k,]))
    
print("Weights prev i<0");
print(sum(abs(t(ORG[,1:L]))^2))
print("Weights prev i>=0");
print(sum(abs(t(ORG[,L+1:(L+1)]))^2))
    
print("Weights i<0");
print(sum(abs(t(RES[,1:L]))^2))
print("Weights i>=0");
print(sum(abs(t(RES[,L+1:(L+1)]))^2))
    
par(mfrow=c(2,2))
plot(thetas, param, t='l', xlab="frequency", ylab="magnitude")
title("Rotation")
plot(0:(2*L)-L,weights,xlab="l",ylab="L_2 norm",type='h');
  title("Rotation")
    
plot(0:(2*L)-L,colMeans(abs(ORG[,])^2),xlab="l",ylab="L_2 norm",type='h');
title("Original filters' magnitudes")
plot(0:(2*L)-L,colMeans(abs(RES[,])^2),xlab="l",ylab="L_2 norm",type='h');
title("Rotated filters' magnitudes")
par(mfrow=c(1,1))

  }
  E
}
