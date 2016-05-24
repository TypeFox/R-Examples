shrinkcovmat.equal <-
function(data,centered=FALSE){
if(!is.matrix(data)) data <- as.matrix(data)
p <- nrow(data)
N <- ncol(data)
centered <- as.logical(centered)
if(centered!=TRUE && centered!=FALSE)
  stop("'centered' must be either 'TRUE' or 'FALSE'")
if(!centered){
if(N < 4) stop("The number of columns should be greater than 3")
datacen <- data-rowMeans(data)
Sigmasam <- tcrossprod(datacen)/(N-1)
trSigmahat <- sum(diag(Sigmasam))
nuhat <- trSigmahat/p
Q <- sum(colSums(datacen^2)^2)/(N-1)
trSigma2hat <- (N-1)/(N*(N-2)*(N-3))*((N-1)*(N-2)*sum(Sigmasam^2)+(trSigmahat)^2-N*Q)
lambdahat <- (trSigmahat^2+trSigma2hat)/(N*trSigma2hat+(p-N+1)/p*trSigmahat^2)
lambdahat <- min(lambdahat,1)
            } else {
if(N < 2) stop("The number of columns should be greater than 1")              
Sigmasam <- tcrossprod(data)/N
trSigmahat <- sum(diag(Sigmasam))
nuhat <- trSigmahat/p
trSigma2hat <- 0
for(i in 1:(N-1)) trSigma2hat <- sum(crossprod(data[,i],data[,(i+1):N])^2) + trSigma2hat              
trSigma2hat <- 2*trSigma2hat/N/(N-1)
lambdahat <- (trSigmahat^2+trSigma2hat)/((N+1)*trSigma2hat+(p-N)/p*trSigmahat^2)
lambdahat <- min(lambdahat,1)
                   }
if(lambdahat<1) {
Sigmahat <- (1-lambdahat)*Sigmasam
diag(Sigmahat) <- nuhat*lambdahat+diag(Sigmahat)
} else Sigmahat <- diag(lambdahat*nuhat,p)
Target <- diag(nuhat,p)
ans <- list(Sigmahat=Sigmahat,lambdahat=lambdahat,Sigmasample=Sigmasam,Target=Target,centered=centered)
class(ans) <- "covmat"
ans
}
