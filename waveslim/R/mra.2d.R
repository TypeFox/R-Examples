mra.2d <-
  function(x, wf="la8", J=4, method="modwt", boundary="periodic")
{
  m <- dim(x)[1]
  n <- dim(x)[2]

  switch(boundary,
    "periodic" = invisible(),
    stop("Invalid boundary rule in mra"))

  if(method == "modwt") {
    x.wt <- modwt.2d(x, wf, J, "periodic")
  } else {
    x.wt <- dwt.2d(x, wf, J, "periodic")
  }
  
  x.mra <- vector("list", 3*J+1)

  ## Smooth
  zero <- vector("list", 3*J+1)
  names(zero) <-
    c(matrix(rbind(paste("LH", 1:J, sep=""), paste("HL", 1:J, sep=""),
                   paste("HH", 1:J, sep="")), nrow=1), paste("LL", J, sep=""))
  attr(zero, "J") <- J
  attr(zero, "wavelet") <- wf
  attr(zero, "boundary") <- boundary
  zero[[3*J+1]] <- x.wt[[3*J+1]]
  if(method == "modwt") {
    for(k in 1:(3*J))
      zero[[k]] <- matrix(0, m, n)
    x.mra[[3*J+1]] <- imodwt.2d(zero)
  } else {
    for(k in 1:J)
      zero[[3*(k-1)+1]] <- zero[[3*(k-1)+2]] <- zero[[3*k]] <-
        matrix(0, m/2^k, n/2^k)
    x.mra[[3*J+1]] <- idwt.2d(zero)
  }

## Details
  for(j in (3*J):1) {
    Jj <- ceiling(j/3)
    zero <- vector("list", 3*Jj+1)
    names(zero) <-
      c(matrix(rbind(paste("LH", 1:Jj, sep=""), paste("HL", 1:Jj, sep=""),
                     paste("HH", 1:Jj, sep="")), nrow=1),
        paste("LL", Jj, sep=""))
    attr(zero, "J") <- Jj
    attr(zero, "wavelet") <- wf
    attr(zero, "boundary") <- boundary
    zero[[j]] <- x.wt[[j]]
    if(method == "modwt") {
      for(k in names(zero)[-charmatch(names(zero)[j], names(zero))])
        zero[[k]] <- matrix(0, m, n)
      x.mra[[j]] <- imodwt.2d(zero)
    } else {
      for(k in 1:Jj)
        zero[[3*(k-1)+1]] <- zero[[3*(k-1)+2]] <- zero[[3*k]] <-
          matrix(0, m/2^k, n/2^k)
      zero[[3*Jj+1]] <- matrix(0, m/2^Jj, n/2^Jj)
      zero[[j]] <- x.wt[[j]]
      x.mra[[j]] <- idwt.2d(zero)
    }
  }

  names(x.mra) <-
    c(matrix(rbind(paste("LH", 1:J, sep=""), paste("HL", 1:J, sep=""),
                   paste("HH", 1:J, sep="")), nrow=1), paste("LL", Jj, sep=""))
   return(x.mra)
}

