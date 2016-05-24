###########################################################################
###########################################################################
###########################################################################

dwt.3d <- function(x, wf, J=4, boundary="periodic")
{
  nx <- dim(x)[1]
  storage.mode(nx) <- "integer"
  ny <- dim(x)[2]
  storage.mode(ny) <- "integer"
  nz <- dim(x)[3]
  storage.mode(nz) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  z <- array(0, dim=c(nx,ny,nz)/2)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 7*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("three_D_dwt", "cube"=as.double(x), "NX"=nx, "NY"=ny,
              "NZ"=nz, "filter.length"=L, "hpf"=h, "lpf"=g, "LLL"=z,
              "HLL"=z, "LHL"=z, "LLH"=z, "HHL"=z, "HLH"=z, "LHH"=z,
              "HHH"=z, PACKAGE="waveslim")[8:15]
    if(j < J) {
      index <- (7*(j-1)+1):(7*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      nx <- dim(x)[1]
      storage.mode(nx) <- "integer"
      ny <- dim(x)[2]
      storage.mode(ny) <- "integer"
      nz <- dim(x)[3]
      storage.mode(nz) <- "integer"
      z <- array(0, dim=c(nx,ny,nz)/2)
      storage.mode(z) <- "double"
    }
    else {
      index <- (7*(j-1)+1):(7*j+1)
      x.wt[index] <- out[c(2:8,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:8,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  class(x.wt) <- "dwt.3d"
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  return(x.wt)
}

###########################################################################
###########################################################################
###########################################################################

idwt.3d <- function(y)
{
  J <- attributes(y)$J
  LLL <- paste("LLL", J, sep="")

  wf <- attributes(y)$wavelet
  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  y.in <- y$LLL

  for(j in J:1) {
    HLL <- paste("HLL", j, sep="")
    LHL <- paste("LHL", j, sep="")
    LLH <- paste("LLH", j, sep="")
    HHL <- paste("HHL", j, sep="")
    HLH <- paste("HLH", j, sep="")
    LHH <- paste("LHH", j, sep="")
    HHH <- paste("HHH", j, sep="")
    
    nx <- dim(y.in)[1]
    storage.mode(nx) <- "integer"
    ny <- dim(y.in)[2]
    storage.mode(ny) <- "integer"
    nz <- dim(y.in)[3]
    storage.mode(nz) <- "integer"

    z <- array(0, dim=2*c(nx, ny, nz))
    storage.mode(z) <- "double"

    out <- .C("three_D_idwt", as.double(y.in), as.double(y[[HLL]]),
              as.double(y[[LHL]]), as.double(y[[LLH]]),
              as.double(y[[HHL]]), as.double(y[[HLH]]),
              as.double(y[[LHH]]), as.double(y[[HHH]]), 
              nx, ny, nz, L, h, g, "Y"=z, PACKAGE="waveslim")

    y.in <- out$Y
  }
  zapsmall(y.in)
}

###########################################################################
###########################################################################
###########################################################################

modwt.3d <- function(x, wf, J=4, boundary="periodic")
{
  nx <- dim(x)[1]
  storage.mode(nx) <- "integer"
  ny <- dim(x)[2]
  storage.mode(ny) <- "integer"
  nz <- dim(x)[3]
  storage.mode(nz) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  z <- array(0, dim=c(nx,ny,nz))
  storage.mode(z) <- "double"

  x.wt <- vector("list", 7*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("three_D_modwt", "cube"=as.double(x), "NX"=nx, "NY"=ny,
              "NZ"=nz, "J"=j, "filter.length"=L, "hpf"=h, "lpf"=g,
              "LLL"=z, "HLL"=z, "LHL"=z, "LLH"=z, "HHL"=z, "HLH"=z,
              "LHH"=z, "HHH"=z, PACKAGE="waveslim")[9:16]
    if(j < J) {
      index <- (7*(j-1)+1):(7*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      nx <- dim(x)[1]
      storage.mode(nx) <- "integer"
      ny <- dim(x)[2]
      storage.mode(ny) <- "integer"
      nz <- dim(x)[3]
      storage.mode(nz) <- "integer"
      z <- array(0, dim=c(nx,ny,nz))
      storage.mode(z) <- "double"
    }
    else {
      index <- (7*(j-1)+1):(7*j+1)
      x.wt[index] <- out[c(2:8,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:8,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  class(x.wt) <- "modwt.3d"
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  return(x.wt)
}

###########################################################################
###########################################################################
###########################################################################

imodwt.3d <- function(y)
{
  J <- attributes(y)$J
  LLL <- paste("LLL", J, sep="")

  wf <- attributes(y)$wavelet
  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  y.in <- y$LLL

  for(j in J:1) {
    HLL <- paste("HLL", j, sep="")
    LHL <- paste("LHL", j, sep="")
    LLH <- paste("LLH", j, sep="")
    HHL <- paste("HHL", j, sep="")
    HLH <- paste("HLH", j, sep="")
    LHH <- paste("LHH", j, sep="")
    HHH <- paste("HHH", j, sep="")
    
    nx <- dim(y.in)[1]
    storage.mode(nx) <- "integer"
    ny <- dim(y.in)[2]
    storage.mode(ny) <- "integer"
    nz <- dim(y.in)[3]
    storage.mode(nz) <- "integer"

    z <- array(0, dim=c(nx, ny, nz))
    storage.mode(z) <- "double"

    out <- .C("three_D_imodwt", as.double(y.in), as.double(y[[HLL]]),
              as.double(y[[LHL]]), as.double(y[[LLH]]),
              as.double(y[[HHL]]), as.double(y[[HLH]]),
              as.double(y[[LHH]]), as.double(y[[HHH]]), 
              nx, ny, nz, j, L, h, g, "Y"=z, PACKAGE="waveslim")

    y.in <- out$Y
  }
  zapsmall(y.in)
}

###########################################################################
###########################################################################
###########################################################################

mra.3d <-
  function(x, wf="la8", J=4, method="modwt", boundary="periodic")
{
  nx <- dim(x)[1]
  ny <- dim(x)[2]
  nz <- dim(x)[3]

  if(method == "modwt") {
    x.wt <- modwt.3d(x, wf, J, "periodic")
  } else {
    x.wt <- dwt.3d(x, wf, J, "periodic")
  }
  
  x.mra <- vector("list", 7*J+1)
  names(x.mra) <-
    c(matrix(rbind(paste("HLL", 1:J, sep=""), paste("LHL", 1:J, sep=""),
                   paste("LLH", 1:J, sep=""), paste("HHL", 1:J, sep=""),
                   paste("HLH", 1:J, sep=""), paste("LHH", 1:J, sep=""),
                   paste("HHH", 1:J, sep="")), nrow=1),
      paste("LLL", J, sep=""))

  ## Smooth
  zero <- vector("list", 7*J+1)
  names(zero) <- names(x.mra)
  attr(zero, "J") <- J
  attr(zero, "wavelet") <- wf
  attr(zero, "boundary") <- boundary
  zero[[7*J+1]] <- x.wt[[7*J+1]]
  if(method == "modwt") {
    class(x.wt) <- "modwt.3d"
    for(k in 1:(7*J))
      zero[[k]] <- array(0, dim=c(nx,ny,nz))
    x.mra[[7*J+1]] <- imodwt.3d(zero)
  } else {
    class(x.wt) <- "dwt.3d"
    for(k in 1:J)
      zero[[7*(k-1)+1]] <- zero[[7*(k-1)+2]] <- zero[[7*(k-1)+3]] <-
        zero[[7*(k-1)+4]] <- zero[[7*(k-1)+5]] <- zero[[7*(k-1)+6]] <-
          zero[[7*k]] <- array(0, dim=c(nx,ny,nz)/2^k)
    x.mra[[7*J+1]] <- idwt.3d(zero)
  }

  ## Details
  for(j in (7*J):1) {
    Jj <- ceiling(j/7)
    zero <- vector("list", 7*Jj+1)
    names(zero) <-
      c(matrix(rbind(paste("HLL", 1:Jj, sep=""), paste("LHL", 1:Jj, sep=""),
                     paste("LLH", 1:Jj, sep=""), paste("HHL", 1:Jj, sep=""),
                     paste("HLH", 1:Jj, sep=""), paste("LHH", 1:Jj, sep=""),
                     paste("HHH", 1:Jj, sep="")), nrow=1),
        paste("LLL", Jj, sep=""))
    attr(zero, "J") <- Jj
    attr(zero, "wavelet") <- wf
    attr(zero, "boundary") <- boundary
    zero[[j]] <- x.wt[[j]]
    if(method == "modwt") {
      for(k in names(zero)[-charmatch(names(zero)[j], names(zero))])
        zero[[k]] <- array(0, dim=c(nx,ny,nz))
      x.mra[[j]] <- imodwt.3d(zero)
    } else {
      for(k in 1:Jj)
        zero[[7*(k-1)+1]] <- zero[[7*(k-1)+2]] <- zero[[7*(k-1)+3]] <-
          zero[[7*(k-1)+4]] <- zero[[7*(k-1)+5]] <- zero[[7*(k-1)+6]] <-
            zero[[7*k]] <- array(0, dim=c(nx,ny,nz)/2^k)
      zero[[7*Jj+1]] <- array(0, dim=c(nx,ny,nz)/2^Jj)
      zero[[j]] <- x.wt[[j]]
      x.mra[[j]] <- idwt.3d(zero)
    }
  }
  return(x.mra)
}
