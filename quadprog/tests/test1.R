library(quadprog)

FullAmat <- function(Dmat, Amat, Aind){
  res <- matrix(0, nrow=NROW(Dmat), ncol=NCOL(Amat))
  ii <- as.vector(Aind[-1,])
  ii <- cbind(ii, rep(1:NCOL(Aind), each=NROW(Aind)-1))
  ind <- ii[,1] != 0
  res[ii[ind,]] <- as.vector(Amat)[ind]
  res
}

Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Aind <- c(2,2,2)
Aind <- rbind(Aind,c(1,1,2))
Aind <- rbind(Aind,c(2,2,3))
Amat <- c(-4,2,-2)
Amat <- rbind(Amat,c(-3,1,1))
bvec  <- c(-8,2,0)
res<-solve.QP.compact(Dmat,dvec,Amat,Aind,bvec=bvec)
print(res)
res<-solve.QP.compact(solve(chol(Dmat)),dvec,Amat,Aind,bvec=bvec,fac=T)
print(res)

print(crv1 <- (crossprod(res$unc, Dmat)/2-dvec)%*%res$unc)
print(crv2 <- (crossprod(res$solution, Dmat)/2-dvec)%*%res$solution)
print(res$value)
print(crv2 >= crv1)
Amat <- FullAmat(Dmat, Amat, Aind)
print(all(crossprod(Amat, res$solution) >= bvec))

Dmat <- matrix(c(4,-2,-2,4),2,2)
dvec <- c(-6,0)
Amat <- c(1,1,1)
Amat <- rbind(Amat,Amat)
Aind <- c(1,1,2)
Aind <- rbind(Aind,c(1,2,1))
Aind <- rbind(Aind,c(0,0,2))
bvec <- c(0,0,2)
res<-solve.QP.compact(Dmat,dvec,Amat,Aind,bvec=bvec)
print(res)
res<-solve.QP.compact(solve(chol(Dmat)),dvec,Amat,Aind,bvec=bvec,fac=T)
print(res)

print(crv1 <- (crossprod(res$unc, Dmat)/2-dvec)%*%res$unc)
print(crv2 <- (crossprod(res$solution, Dmat)/2-dvec)%*%res$solution)
print(res$value)
print(crv2 >= crv1)
Amat <- FullAmat(Dmat, Amat, Aind)
print(all(crossprod(Amat, res$solution) >= bvec))

