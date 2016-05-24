Deltakl <- function(N, n, p){
Ind <- Ik(N,n)
P1 <- as.matrix(Pik(p, Ind))
Delta <-Pikl(N,n,p)-(t(P1)%*%P1)
return(Delta)
}
