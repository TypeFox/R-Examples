library("RUnit")
library("aucm")

test.srauc <- function() {

RNGkind("Mersenne-Twister", "Inversion")
#RNGkind("Marsaglia-Multicarry", "Kinderman-Ramage") 

tol=1e-6


dat = sim.dat.1(n=200,seed=1, add.outliers=TRUE)
fit=srauc(y~x1+x2, dat, lambda=1, method="tron")
checkEqualsNumeric(coef(fit), c(3.2934696, 0.3524196), tol=tol)


## the following gives different results on 32 bit windows and 64 bit linux
#dat = sim.disc (n=200,seed=1, add.outliers=TRUE) 
#fit=srauc(y~x1+x2, dat, kernel="rbf", para=0.363, lambda=11.3, method="tron")
#checkEqualsNumeric(coef(fit)[1:2], c(3794.167, -1857.790), tol=tol)# 32 bit windows results


}
