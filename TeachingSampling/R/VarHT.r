VarHT<-function(y, N, n, p){
Ind <- Ik(N,n)
pi1 <- as.matrix(Pik(p, Ind))
pi2 <- Pikl(N,n,p)
Delta <- Deltakl(N,n,p)
y <- t(as.matrix(y))
ykylexp <- t(y/pi1)%*%(y/pi1)
A <- (Delta)*(ykylexp)
Var <- sum(A)
return(Var)
}