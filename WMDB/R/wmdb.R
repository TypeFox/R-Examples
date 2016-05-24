wmahalanobis <- function (x, center, cov, weight) 
{
    if (is.vector(x)) 
        x=matrix(x, ncol = length(x))
    else x=as.matrix(x)
    x <- sweep(x, 2, center)
    cov <- weight%*%solve(cov)
    retval <- diag(x%*%cov%*%t(x))
    retval
}


wmd <- function
(TrnX, TrnG, Tweight=NULL, TstX = NULL, var.equal = F){
if ( is.factor(TrnG) == FALSE){
mx <-nrow(TrnX); mg <- nrow(TrnG)
TrnX<-rbind(TrnX, TrnG)
TrnG<-factor(rep(1:2, c(mx, mg)))
}
if (is.null(TstX) == TRUE) TstX <- TrnX
if (is.vector(TstX) == TRUE) TstX <- t(as.matrix(TstX))
else if (is.matrix(TstX) != TRUE)
TstX <- as.matrix(TstX)
if (is.matrix(TrnX) != TRUE) TrnX <- as.matrix(TrnX)
if (is.null(Tweight)==TRUE) Tweight=cor(TstX)
nx <- nrow(TstX)
blong <- matrix(rep(0, nx), nrow=1,
dimnames=list("blong", 1:nx))
g <- length(levels(TrnG))
mu <- matrix(0, nrow=g, ncol=ncol(TrnX))
for (i in 1:g)
mu[i,] <- colMeans(TrnX[TrnG==i,])
D <-matrix(0, nrow=g, ncol=nx)
if (var.equal == TRUE || var.equal == T){
for (i in 1:g)
D[i,]<- wmahalanobis(TstX, mu[i,], var(TrnX),Tweight)
}
else{
for (i in 1:g)
D[i,]<-wmahalanobis(TstX,mu[i,],var(TrnX[TrnG==i,]),Tweight)
}
for (j in 1:nx){
dmin <- Inf
for (i in 1:g)
if (D[i,j] < dmin){
dmin <- D[i,j]; blong[j] <- i
}
}
print(blong)
print("num of wrong judgement")
print(which(blong!=TrnG))
print("samples divided to")
print(blong[which(blong!=TrnG)])
print("samples actually belongs to")
print(TrnG[which(blong!=TrnG)])
print("percent of right judgement")
print(1-length(which(blong!=TrnG))/length(blong))
}


dbayes <- function
(TrnX, TrnG, p = rep(1, length(levels(TrnG))),
TstX = NULL, var.equal = FALSE){
if ( is.factor(TrnG) == FALSE){
mx <- nrow(TrnX); mg <- nrow(TrnG)
TrnX <- rbind(TrnX, TrnG)
TrnG <- factor(rep(1:2, c(mx, mg)))
}
if (is.null(TstX) == TRUE) TstX <- TrnX
if (is.vector(TstX) == TRUE) TstX <- t(as.matrix(TstX))
else if (is.matrix(TstX) != TRUE)
TstX <- as.matrix(TstX)
if (is.matrix(TrnX) != TRUE) TrnX <- as.matrix(TrnX)
nx <- nrow(TstX)
blong <- matrix(rep(0, nx), nrow=1,
dimnames=list("blong", 1:nx))
g <- length(levels(TrnG))
mu <- matrix(0, nrow=g, ncol=ncol(TrnX))
for (i in 1:g)
mu[i,] <- colMeans(TrnX[TrnG==i,])
D <- matrix(0, nrow=g, ncol=nx)
if (var.equal == TRUE || var.equal == T){
for (i in 1:g){
d2 <- mahalanobis(TstX, mu[i,], var(TrnX))
D[i,] <- d2 - 2*log(p[i])
}
}
else{
for (i in 1:g){
S <- var(TrnX[TrnG==i,])
d2 <- mahalanobis(TstX, mu[i,], S)
D[i,] <- d2 - 2*log(p[i])-log(det(S))
}
}
for (j in 1:nx){
dmin <- Inf
for (i in 1:g)
if (D[i,j] < dmin){
dmin <- D[i,j]; blong[j] <- i
}
}
print(blong)
print("num of wrong judgement")
print(which(blong!=TrnG))
print("samples divided to")
print(blong[which(blong!=TrnG)])
print("samples actually belongs to")
print(TrnG[which(blong!=TrnG)])
print("percent of right judgement")
print(1-length(which(blong!=TrnG))/length(blong))
}


