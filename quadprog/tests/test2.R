library(quadprog)

Dmat <- matrix(c(4,-2,-2,4),2,2)
dvec <- c(-6,0)
Amat <- matrix(c(1,0,0,1,1,1),2,3)
bvec <- c(0,0,2)
res<-solve.QP(Dmat,dvec,Amat,bvec=bvec)
print(res)
res<-solve.QP(solve(chol(Dmat)),dvec,Amat,bvec=bvec,fac=T)
print(res)

print(crv1 <- (crossprod(res$unc, Dmat)/2-dvec)%*%res$unc)
print(crv2 <- (crossprod(res$solution, Dmat)/2-dvec)%*%res$solution)
print(res$value)
print(crv2 >= crv1)
print(all(crossprod(Amat, res$solution) >= bvec))


Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat  <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec  <- c(-8,2,0)
res<-solve.QP(Dmat,dvec,Amat,bvec=bvec)
print(res)
res<-solve.QP(solve(chol(Dmat)),dvec,Amat,bvec=bvec,fac=T)
print(res)

print(crv1 <- (crossprod(res$unc, Dmat)/2-dvec)%*%res$unc)
print(crv2 <- (crossprod(res$solution, Dmat)/2-dvec)%*%res$solution)
print(res$value)
print(crv2 >= crv1)
print(all(crossprod(Amat, res$solution) >= bvec))

