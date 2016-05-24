`dist.xyz` <-
function(a, b=NULL, all.pairs=TRUE, ncore=1, nseg.scale=1){

  ## if 'a' is a vector (or matrix) and 
  ## 'b' is a matrix
  ## compare (each row of) 'a' to all rows in 'b'

  ## if 'a' is a matrix and 'b' is NULL
  ## call 'dist' on 'a'

  ## if 'a' is a vector and 'b' is NULL
  ## make 'a' a 3 col matrix and call 'dist' 

  # Parallelized by parallel package (Fri Jul  5 19:58:32 EDT 2013)
  ncore <- setup.ncore(ncore)
  if(ncore > 1) {
     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }

  if(is.vector(a)) {
    a <- matrix(a, ncol=3, byrow=TRUE)
  } else {
    a <- as.matrix(a)
  }
  
  if(is.null(b)) {
    return(as.matrix(dist(a)))
  } else {
    if(is.vector(b)) {
      b <- matrix(b, ncol=3, byrow=TRUE)
    } else {
      b <- as.matrix(b)
    }
  }

  dima <- ncol(a)
  dimb <- ncol(b)
  if(dima != dimb)
    stop("Dimension miss-match of input 'a' and 'b'")
    
  if(dima != 3) {
    warning(paste("input does not have three columns: assuming you want",
                  dima, "dimensional distances")) 
  }

  if(!all.pairs) {
    ## distance between coresponding rows
    d <- rep( NA, max(nrow(a), nrow(b)) )
    ind <- 1:min(nrow(a), nrow(b))
    d[ind] <- sqrt( rowSums((a[ind,] - b[ind,])^2) )
    ##    return( sqrt( rowSums((a - b)^2) ) )
    return(d)
  } else {
    
    if(ncore > 1) {
       RLIMIT = floor(R_NCELL_LIMIT / nrow(b))
       nDataSeg = floor((nrow(a)-1)/RLIMIT) + 1
       nDataSeg = floor(nDataSeg * nseg.scale)
       lenSeg = floor(nrow(a)/nDataSeg)
       d.l <- NULL
       for(i in 1:nDataSeg) {
          istart = (i-1)*lenSeg + 1
          iend = if(i<nDataSeg) i*lenSeg else nrow(a)
          d.l <- c(d.l, mclapply(istart:iend, function(j) {
             sqrt(colSums((a[j,] - t(b))^2)) } ) )
       }   
       d <- do.call(rbind, d.l)
    } else {
       d <- matrix(0, nrow=nrow(a), ncol=nrow(b))
       for(i in 1:nrow(a)){
         d[i,] <- sqrt(colSums((a[i,] - t(b))^2))
       }
    }
    return(d)
  }
}

