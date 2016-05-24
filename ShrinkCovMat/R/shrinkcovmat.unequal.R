shrinkcovmat.unequal <-
function(data,centered=FALSE){
if(!is.matrix(data)) data <- as.matrix(data)
p <- nrow(data)
N <- ncol(data)
centered <- as.logical(centered)
if(centered!=TRUE && centered!=FALSE)
  stop("'centered' must be either 'TRUE' or 'FALSE'")
if(!centered){
if(N < 4) stop("the number of columns should be greater than 3")
datacen <- data-rowMeans(data)
Sigmasam <- tcrossprod(datacen)/(N-1)
trSigmahat <- sum(diag(Sigmasam))
Q <- sum(colSums(datacen^2)^2)/(N-1)
trSigma2hat <- (N-1)/(N*(N-2)*(N-3))*((N-1)*(N-2)*sum(Sigmasam^2)+(trSigmahat)^2-N*Q)
help1 <- help21 <- help22 <- help3 <- rep(0,p)
for(i in 1:(N-1)){
data2 <- matrix(data[,(i+1):N],p,N-i)
help1 <- rowSums(data[,i]*data2)+help1
help21 <- rowSums(data[,i]^3*data2)+help21
help22 <- rowSums(data2^3*data[,i])+help22
help3 <- rowSums(data[,i]^2*data2^2)+help3
}
trD21 <- 2*sum(help3)/N/(N-1)
trD22 <- 2*(sum(help1*rowSums(data^2))-sum(help21+help22))
trD23 <- 4*(sum(help1^2)-sum(help3)-trD22)
trD22 <- trD22/N/(N-1)/(N-2)
trD23 <- trD23/N/(N-1)/(N-2)/(N-3)
trDs2 <- trD21-2*trD22+trD23   
lambdahat <- (trSigmahat^2+trSigma2hat-2*trDs2)/(N*trSigma2hat+trSigmahat^2-(N+1)*trDs2)
lambdahat <- max(0,min(lambdahat,1))
} else {
if(N < 2) stop("the number of columns should be greater than 1")   
Sigmasam <- tcrossprod(data)/N
trSigmahat <- sum(diag(Sigmasam))
trSigma2hat <- trDs2 <- 0
for(i in 1:(N-1)) {
trSigma2hat <- sum(crossprod(data[,i],data[,(i+1):N])^2) + trSigma2hat              
trDs2 <- sum((data[,i]*data[,(i+1):N])^2) + trDs2
}
trSigma2hat <- 2*trSigma2hat/N/(N-1)
trDs2 <- 2*trDs2/N/(N-1)
lambdahat <- (trSigmahat^2+trSigma2hat-2*trDs2)/((N+1)*trSigma2hat+trSigmahat^2-(N+2)*trDs2)
lambdahat <- max(0,min(lambdahat,1))
}
diagSam <- diag(Sigmasam)
if(lambdahat<1) {
Sigmahat <- (1-lambdahat)*Sigmasam
diag(Sigmahat) <- lambdahat+diagSam
} else Sigmahat <- diag(lambdahat*diagSam,p)
Target <- diag(diagSam,p)
ans <- list(Sigmahat=Sigmahat,lambdahat=lambdahat,Sigmasample=Sigmasam,Target=Target,centered=centered)
class(ans) <- "covmat"
ans
}
