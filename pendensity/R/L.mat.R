#Generates the difference matrix of order m.
L.mat <- function(penden.env) {
  m <- get("m",penden.env)
  K <- get("K",penden.env)
  M <- get("M",penden.env)
  if(m==1) {
  L <- diag(1,K)
  L.1 <- diag(-1,K,K-1)
  L.2 <- matrix(0,K,1)
  L1 <- cbind(L.2,L.1)
  L <- L+L1
  #M <- floor(K/2) + 1
  L <- cbind(L[,1:(M-1)],L[,(M+1):K])
  L <- L[1:(K-1),]
  return(L)
 }
 if(m==2) {
  L <- diag(1,K,K)
  L1 <- diag(-2,K,(K-1))
  L2 <- diag(1,K,(K-2))
  L.1 <- matrix(0,K,1)
  L1 <- cbind(L.1,L1)
  L2 <- cbind(L.1,L.1,L2)
  L <- L+L1+L2
  L <- cbind(L[,1:(M-1)],L[,(M+1):K])
  L <- L[1:(K-2),]
  return(L)
 }
 if(m==3) {
  L <- diag(1,(K-3),(K-3))
  L.help <- matrix(0,(K-3),1)
  L1 <- diag(-3,(K-3),(K-3))
  M1 <- cbind(L,L.help,L.help,L.help)
  M2 <- cbind(L.help,L1,L.help,L.help)
  M3 <- cbind(L.help,L.help,-L1,L.help)
  M4 <- cbind(L.help,L.help,L.help,-L)
  L <- (M1+M2+M3+M4)
  L <- cbind(L[,1:(M-1)],L[,(M+1):K])
  return(L)
  }
 if(m==4) {
  L <- diag(1,(K-4),(K-4))
  L.help <- matrix(0,(K-4),1)
  L1 <- diag(-4,(K-4),(K-4))
  L2 <- diag(6,(K-4),(K-4))
  M1 <- cbind(L,L.help,L.help,L.help,L.help)
  M2 <- cbind(L.help,L1,L.help,L.help,L.help)
  M3 <- cbind(L.help,L.help,L2,L.help,L.help)
  M4 <- cbind(L.help,L.help,L.help,L1,L.help)
  M5 <- cbind(L.help,L.help,L.help,L.help,L)
  L <- (M1+M2+M3+M4+M5)
  L <- cbind(L[,1:(M-1)],L[,(M+1):K])
  return(L)
 }
}
