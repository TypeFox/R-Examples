"vma.sim" <-
function(psi,a){
    if (!is.null(psi) && class(psi)!="array")
	stop("Psi must be enterd as NULL or array with dimension (k*k*Q)")
      if (is.null(psi)) return(a)
      n <- NROW(a)
      k <- NCOL(a)
      dimPsi <- dim(psi)
      if (length(dimPsi) !=3 || dimPsi[1] !=dimPsi[2] || dimPsi[1] !=k)
        stop("invalid dimensions for Psi")
      q <-  dimPsi[3]
      extend.psi <- array(numeric(n*k*k), dim=c(k, k, n))
      extend.psi[,,1:q] <- psi
      sim.vma <- matrix(numeric(0),nrow=n,ncol=k)
	  for (i in 1:n){
	    out <- 0
	     for (j in 1:i){
	      out=out+crossprod(t(extend.psi[,,j]),a[i-j+1,])
	     }
          sim.vma[i,] <- out
        }
      Q <- q-1
    a <- matrix(sim.vma[-(1:Q),],ncol=k)
  return(matrix(a[1:(n-Q),],ncol=k))
}