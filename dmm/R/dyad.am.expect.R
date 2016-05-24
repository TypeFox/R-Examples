dyad.am.expect <-
function(am, gls,dmeopt){
# dyad.am.expect()
# evaluate parts of dyadic model equation  including emat and emat.qr
# am is an ante-model object
# am$z is a list of Z matrices
# am$rel is a list of (pseudo) relationship matrices
# returns a list object containing emat, emat.qr, and cnames(col names for emat)
#
# Note: This code sets the order of components in siga[,]
#
# setup m matrix
  m <- diag(am$n) - am$x %*% ginv(am$x)
# cat("m:\n")
# print(m)
#
  cnames <- "VarE(I)"
  emat <- matrix(0,am$n * am$n, am$v)
  vmat <- matrix(0,am$n * am$n, am$v)
  zaz <- matrix(0,am$n,am$n)
  icol <- 1

  if(any(am$components == "VarE(I)")){
    zaz <- am$z$i %*% am$rel$e %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarE(I)"
    icol <- icol + 1
  }

  if(any(am$components == "VarG(Ia)")) {
    zaz <- am$z$i %*% am$rel$a %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Ia)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Id)")) {
    zaz <- am$z$i %*% am$rel$d %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Id)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Ia:a)")) {
    zaz < am$z$i %*% am$rel$aa %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Ia:a)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Ia:d)")) {
    zaz <- am$z$i %*% am$rel$ad %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Ia:d)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Id:d)")) {
    zaz <- am$z$i %*% am$rel$dd %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Id:d)"
    icol <- icol + 1
  }

  if(any(am$components == "VarGs(Ia)")) {
    zaz <- am$z$i %*% am$rel$s %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarGs(Ia)"
    icol <- icol + 1
  }
  
  if(any(am$components == "VarE(M)")) {
    zaz <- am$z$m %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarE(M)"
    icol <- icol + 1
  }
  if(any(am$components == "VarE(C)")) {
    zaz <- am$z$c %*% t(am$z$c)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarE(C)"
    icol <- icol + 1
  }

  if(any(am$components == "VarE(M&!C)")) {
    zaz <- ((am$z$m %*% t(am$z$m)) & !(am$z$c %*% t(am$z$c)) + 0)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarE(M&!C)"
    icol <- icol + 1
  }

  if(any(am$components == "VarE(M&C)")) {
    zaz <- ((am$z$m %*% t(am$z$m)) & (am$z$c %*% t(am$z$c)) + 0)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarE(M&C)"
    icol <- icol + 1
  }

  if(any(am$components == "VarG(Ma)")) {
    zaz <- am$z$m %*% am$rel$a %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Ma)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Md)")) {
    zaz <- am$z$m %*% am$rel$d %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Md)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Ma:a)")) {
    zaz <- am$z$m %*% am$rel$aa %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Ma:a)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Ma:d)")) {
    zaz <- am$z$m %*% am$rel$ad %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Ma:d)"
    icol <- icol + 1
  }
  if(any(am$components == "VarG(Md:d)")) {
    zaz <- am$z$m %*% am$rel$dd %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarG(Md:d)"
    icol <- icol + 1
  }

  if(any(am$components == "VarGs(Ma)")) {
    zaz <- am$z$m %*% am$rel$s %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "VarGs(Ma)"
    icol <- icol + 1
  }

  if(any(am$components == "CovE(I,M)")) {
    zaz <- am$z$i %*% am$rel$e %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m) 
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovE(I,M)"
    icol <- icol + 1
  }

  if(any(am$components == "CovE(M,I)")) {
    zaz <- am$z$m %*% am$rel$e %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovE(M,I)"
    icol <- icol + 1
  }

  if(any(am$components == "CovE(I,M&!C)")) {
    zaz <- ((am$z$i %*% t(am$z$m)) & !(am$z$c %*% t(am$z$c)) + 0)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovE(I,M&!C)"
    icol <- icol + 1
  }

  if(any(am$components == "CovE(M&!C,I)")) {
    zaz <- ((am$z$m %*% t(am$z$i)) & !(am$z$c %*% t(am$z$c)) + 0)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovE(M&!C,I)"
    icol <- icol + 1
  }

  if(any(am$components == "CovE(I,M&C)")) {
    zaz <- ((am$z$i %*% t(am$z$m)) & (am$z$c %*% t(am$z$c)) + 0)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovE(I,M&C)"
    icol <- icol + 1
  }

  if(any(am$components == "CovE(M&C,I)")) {
    zaz <- ((am$z$m %*% t(am$z$i)) & (am$z$c %*% t(am$z$c)) + 0)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovE(M&C,I)"
    icol <- icol + 1
  }

  if(any(am$components == "CovG(Ia,Ma)")) {
    zaz <- am$z$i %*% am$rel$a %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m) 
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Ia,Ma)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Ma,Ia)")) {
    zaz <- am$z$m %*% am$rel$a %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m)
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Ma,Ia)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Id,Md)")) {
    zaz <- am$z$i %*% am$rel$d %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m ) 
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Id,Md)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Md,Id)")) {
    zaz <- am$z$m %*% am$rel$d %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m )
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Md,Id)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Ia:a,Ma:a)")) {
    zaz <- am$z$i %*% am$rel$aa %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m ) 
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Ia:a,Ma:a)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Ma:a,Ia:a)")) {
    zaz <- am$z$m %*% am$rel$aa %*% t(am$z$i)
    emat[,icol] <- as.vector( m %*% zaz %*% m )
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Ma:a,Ia:a)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Ia:d,Ma:d)")) {
    zaz <- am$z$i %*% am$rel$ad %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m ) 
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Ia:d,Ma:d)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Ma:d,Ia:d)")) {
    zaz <- am$z$m %*% am$rel$ad %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m )
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Ma:d,Ia:d)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Id:d,Md:d)")) {
    zaz <- am$z$i %*% am$rel$dd %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m ) 
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Id:d,Md:d)"
    icol <- icol + 1
  }
  if(any(am$components == "CovG(Md:d,Id:d)")) {
    zaz <- am$z$m %*% am$rel$dd %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*% zaz %*% m )
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovG(Md:d,Id:d)"
    icol <- icol + 1
  }

  if(any(am$components == "CovGs(Ia,Ma)")) {
    zaz <- am$z$i %*% am$rel$s %*% t(am$z$m)
    emat[,icol] <- as.vector(m %*% zaz %*% m ) 
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovGs(Ia,Ma)"
    icol <- icol + 1
  }
  if(any(am$components == "CovGs(Ma,Ia)")) {
    zaz <- am$z$m %*% am$rel$s %*% t(am$z$i)
    emat[,icol] <- as.vector(m %*%  zaz %*% m )
    vmat[,icol] <- as.vector(zaz)
    cnames[icol] <- "CovGs(Ma,Ia)"
    icol <- icol + 1
  }

  dimnames(emat) <- list(NULL,cnames)
  dimnames(vmat) <- list(NULL,cnames)
#   cat("Checking dyadic model equations:\n")
#   cat("emat:\n")
#   print(emat)
#   emat.sum <- apply(emat,2,sum)
#   cat("column.sums:\n")
#   print(emat.sum)
    emat.mean <- apply(emat,2,mean)
#   cat("column.means:\n")
#   print(emat.mean)
    emat.var <- apply(emat,2,var)
#   cat("column.variances:\n")
#   print(emat.var)
    emat.cor <- cor(emat,emat)
#   cat("column.correlations:\n")
#   print(emat.cor)

#  do QR transform on emat
  emat.qr <- qr(emat)

# check emat for E'E positive definite - ie emat.qr$rank should be am$v
  if(emat.qr$rank != am$v) {
    cat(" Rank of DME matrix = ",emat.qr$rank," no of components = ",am$v,"\n")
    if(dmeopt != "pcr") stop("Dyadic model equations not of full rank: either omit some components or try dmeopt='pcr' \n")
  }

# do qr on vmat  - only need if gls=T
  if(gls) {
    vmat.qr <- qr(vmat)
  }
  else {
    vmat.qr <- NULL
  }

  explist <- list(emat=emat, emat.qr=emat.qr, cnames=cnames, vmat=vmat, vmat.qr=vmat.qr,
                  emat.mean=emat.mean,emat.var=emat.var,emat.cor=emat.cor) 
  return(explist)
}
