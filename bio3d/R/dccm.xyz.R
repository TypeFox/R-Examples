`dccm.xyz` <-
function(x, reference=NULL, grpby=NULL, ncore=1, nseg.scale=1, ... ) {
  xyz <- x
  # Parallelized by parallel package (Wed Dec 12 18:36:39 EST 2012)
  ncore <- setup.ncore(ncore)

  if(is.null(reference)) ref = colMeans(xyz)
  else ref = reference
  dxyz  <- sweep(xyz, 2, ref)

  covmat <- cov(dxyz)
  
  if(!is.null(reference)) {
     # moment instead of covariance
     mxyz <- colMeans(dxyz)
     covmat <- covmat + outer(mxyz, mxyz)
  }

  ccmat <- cov2dccm(covmat, ncore = ncore)
 
  if(is.null(grpby)) {
    return(ccmat)
  } else {
    ##- Group by concetive numbers in 'grpby'
    if( ncol(xyz) != (length(grpby)*3) )
      stop("dimension miss-match in 'xyz' and 'grpby', check lengths")
    
    ##- Bounds of 'grpby' numbers
    inds <- bounds(grpby, dup.inds=TRUE)
    nres <- nrow(inds)

    ##- Per-residue matrix
    m <- matrix(, ncol=nres, nrow=nres)
    ij <- pairwise(nres)

    ##- Max (absolute value) per residue
    for(k in 1 : nrow(ij) ) {
      m[ij[k,1],ij[k,2]] <-
        min( ccmat[ (inds[ij[k,1],"start"]:inds[ij[k,1],"end"]),
                  (inds[ij[k,2],"start"]:inds[ij[k,2],"end"])],
            na.rm=TRUE )
      tmax <- max( ccmat[ (inds[ij[k,1],"start"]:inds[ij[k,1],"end"]),
                  (inds[ij[k,2],"start"]:inds[ij[k,2],"end"])],
            na.rm=TRUE )
      if(tmax > abs(m[ij[k,1],ij[k,2]])) m[ij[k,1],ij[k,2]] = tmax 
    }
#    if( !mask.lower )
    m[lower.tri(m)] = t(m)[lower.tri(m)]
    diag(m) <- 1

    class(m)=c("dccm","matrix")
    return(m)
  }
}

