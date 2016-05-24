###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
glsfit <-
function(Xlist, ylist, Slist, nalist, Psi, onlycoef=TRUE)  {
#
################################################################################
# FUNCTION TO COMPUTE THE GLS ESTIMATE + OPTIONAL INTERMEDIATE PRODUCTS
#
  Sigmalist <- mapply(function(S,na) S+Psi[!na,!na,drop=FALSE],
    Slist,nalist,SIMPLIFY=FALSE)
  Ulist <- lapply(Sigmalist,chol)
  invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
  invtUXlist <- mapply(function(invU,X) crossprod(invU,X),
    invUlist,Xlist,SIMPLIFY=FALSE)
  invtUylist <- mapply(function(invU,y) crossprod(invU,y),
    invUlist,ylist,SIMPLIFY=FALSE)
  invtUX <- do.call("rbind",invtUXlist)
  invtUy <- do.call("rbind",invtUylist)
  coef <- as.numeric(qr.solve(invtUX,invtUy))
#
  if(onlycoef) return(coef)
#
  list(coef=coef,Sigmalist=Sigmalist,Ulist=Ulist,invUlist=invUlist,
    invtUXlist=invtUXlist,invtUX=invtUX,invtUy=invtUy)
}

#
