context("rotation - source")

theta <- 0.225406088569612
q <- c(1, 4, 9, 31, 71, 244, 1047, 5479, 22963, 74368, 171699, 761164)
a <- c(4, 2, 3, 2, 3, 4, 5, 4, 3, 2, 4)
#
mu1 <- c(1-a[1]*theta, a[1]*theta)
mu2 <- c((a[2]+1)*theta, 1-(a[2]+1)*theta)
# noyau de 2 à 1
P21_0 <- as.bigq(c(q[1], q[3]-q[1]), q[3])
P21_1 <- as.bigq(c(0, 1), 1)
P21 <- matrix(c(P21_0, P21_1), byrow=TRUE, nrow=2)
# noyau de 3 à 2
P32_0 <- as.bigq(c(1, 0), 1)
P32_1 <- as.bigq(c(q[4]-q[2], q[2]), q[4])
P32 <- matrix(c(P32_0, P32_1), byrow=TRUE, nrow=2)
#
rho1 <- function(x,y) discrete(x, y, gmp=TRUE)
kanto2 <- kantorovich(as.vector.bigq(P21[1,]), as.vector.bigq(P21[2,]), dist=rho1)
#
rho <- function(x,y){
  if(x==y) return(0) else return(as.character(kanto2))
}
kanto3 <- kantorovich(as.vector.bigq(P32[1,]), as.vector.bigq(P32[2,]), dist=rho)

