

library(LatticeKrig)
options(echo = FALSE)
test.for.zero.flag <- 1
###########################
######## Inverse model or passed "X"  and "U" matrices
##########################

set.seed(122)
lambda <- 0.1
x <- matrix(runif(20 * 2), 20, 2)
nObs <- nrow(x)
LKinfo <- LKrigSetup(x, NC = 4, nlevel = 1, alpha = 1, lambda = lambda, 
	a.wght = 5, NC.buffer = 1)
W <- LKrig.basis(x, LKinfo)
W<- as.matrix(W)
x0<- matrix(runif( 5* 2), 5, 2)
W0 <- LKrig.basis(x0, LKinfo)
W0<- as.matrix(W0)
T.matrix <- do.call(LKinfo$fixedFunction, c(list(x = x, distance.type = LKinfo$distance.type), 
	LKinfo$fixedFunctionArgs))
T0 <- do.call(LKinfo$fixedFunction, c(list(x = x0, distance.type = LKinfo$distance.type), 
	LKinfo$fixedFunctionArgs))
	
A <- matrix(runif(nObs^2), nObs, nObs)
#   A<- diag( 1, nObs)
X <- as.spam(A %*% W)
U <- A %*% T.matrix
c.true <- runif(LKinfo$latticeInfo$m)
y <- X %*% c.true

obj <- LKrig(x, y, LKinfo = LKinfo, X = X, U = U)
rho<- obj$rho.MLE

Q <- LKrig.precision(LKinfo)
Sigma <- solve(Q)
Mlambda <- solve(as.matrix(X %*% Sigma %*% t(X) + lambda * diag(1, nObs)))
PMatrix <- solve(t(U) %*% Mlambda %*% U) %*% t(U) %*% Mlambda
d.check <- PMatrix %*% y
test.for.zero(d.check, obj$d.coef, tag = "checking d coef with inverse X")

findCoef <- solve(t(X) %*% X + lambda * Q) %*% t(X) %*% (diag(1, nObs) - 
	U %*% PMatrix)

c.check <- findCoef %*% y
test.for.zero(c.check, obj$c.coef, tag = "checking c coef with inverse X")

AMatrix<-  W0 %*% findCoef + T0 %*% PMatrix
hold1<- AMatrix%*%y
hold2<- predict( obj, xnew=x0)
test.for.zero(hold1, hold2, tag = "checking A matrix coef with inverse X")

rho<- obj$rho.MLE
sigma<- obj$sigma.MLE

Q <- LKrig.precision(LKinfo)
S <- rho*solve(Q)
temp0<- W0%*% (S)%*%t(W0)
temp1<- W0%*% (S)%*% t( AMatrix%*%X)
temp2<- AMatrix%*%( X%*%S%*%t(X) + diag( sigma^2 ,nObs) )%*%t(AMatrix)
covTest<- temp0 - temp1 - t(temp1) + temp2

test1.se <- predictSE(obj, xnew = x0)
test.for.zero( sqrt(diag( covTest)), test1.se, tag="SE inverse from first formulas" )

# set.seed( 233)
# out<- LKrig.sim.conditional( obj,x.grid=x0, M=2000 )
# test0.sim<- var( t(out$g.draw))
# test.for.zero( sqrt(diag(test0.sim)), test1.se, tol=.06)

