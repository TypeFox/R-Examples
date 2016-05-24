isKerAinKerB <- function(A,B, tol = .Machine$double.eps){
A <- as.matrix(A)
A.svd <- svd(A, nu = 0)
if (sum(A.svd$d > tol * max(A.svd$d))>0)
    {kerA.perp <- A.svd$v[,A.svd$d > tol * max(A.svd$d)]
     ## projector to ker A perp
     Pi.kerA.perp <- kerA.perp%*%solve(t(kerA.perp)%*%kerA.perp, tol = tol)%*%t(kerA.perp)
}else{Pi.kerA.perp <- 0}

B <- as.matrix(B)
B.svd <- svd(B, nu = 0)
if (sum(B.svd$d > tol * max(B.svd$d))>0)
    {kerB.perp <- B.svd$v[,B.svd$d > tol * max(B.svd$d)]
     ## projector to ker B perp
     Pi.kerB.perp <- kerB.perp%*%solve(t(kerB.perp)%*%kerB.perp, tol = tol)%*%t(kerB.perp)
}else{Pi.kerB.perp <- 0}
isTRUE(all.equal(Pi.kerB.perp%*%Pi.kerA.perp, Pi.kerB.perp, tolerance = tol ))
}

setMethod("existsPIC", "L2ParamFamily", function(object, warning = TRUE, tol = .Machine$double.eps){
if(!isKerAinKerB(object@FisherInfo, trafo(object), tol = tol))
  {if(warning)
      warning("trafo of parameter is not (locally) identifyable for parameter theta.")
      return(FALSE)}
  return(TRUE)})

