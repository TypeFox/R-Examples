###########################################################################
###########################################################################
###########################################################################

dwt.2d <- function(x, wf, J=4, boundary="periodic")
{
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  z <- matrix(0, m/2, n/2)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 3*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("two_D_dwt", "Image"=as.double(x), "Rows"=m, "Cols"=n, 
                "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z, "LH"=z,
                "HL"=z, "HH"=z, PACKAGE="waveslim")[7:10]
    if(j < J) {
      index <- (3*j-2):(3*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      m <- dim(x)[1]
      storage.mode(m) <- "integer"
      n <- dim(x)[2]
      storage.mode(n) <- "integer"
      z <- matrix(0, m/2, n/2)
      storage.mode(z) <- "double"
    }
    else {
      index <- (3*j):(3*(j+1)) - 2
      x.wt[index] <- out[c(2:4,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  attr(x.wt, "class") <- "dwt.2d"
  x.wt
}

###########################################################################
###########################################################################
###########################################################################

idwt.2d <- function(y)
{
  J <- attributes(y)$J

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  LL <- paste("LL", J, sep="")
  y.in <- y[[LL]]

  for(j in J:1) {
    LH <- paste("LH", j, sep="")
    HL <- paste("HL", j, sep="")
    HH <- paste("HH", j, sep="")
    
    m <- dim(y.in)[1]
    storage.mode(m) <- "integer"
    n <- dim(y.in)[2]
    storage.mode(n) <- "integer"
    x <- matrix(0, 2*m, 2*n)
    storage.mode(x) <- "double"

    out <- .C("two_D_idwt", as.double(y.in), as.double(y[[LH]]),
              as.double(y[[HL]]), as.double(y[[HH]]), m, n, L, h, g,
              "Y"=x, PACKAGE="waveslim")
    y.in <- out$Y
  }
  zapsmall(y.in)
}

###########################################################################
###########################################################################
###########################################################################

modwt.2d <- function(x, wf, J=4, boundary="periodic")
{
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  z <- matrix(0, m, n)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 3*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("two_D_modwt", "Image"=as.double(x), "Rows"=m, "Cols"=n,
              "Level"=j, "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z,
              "LH"=z, "HL"=z, "HH"=z, PACKAGE="waveslim")[8:11]
    if(j < J) {
      index <- (3*j-2):(3*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out$LL
    }
    else {
      index <- (3*j):(3*(j+1)) - 2
      x.wt[index] <- out[c(2:4,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  x.wt
}

###########################################################################
###########################################################################
###########################################################################

imodwt.2d <- function(y)
{
  J <- attributes(y)$J

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  LL <- paste("LL", J, sep="")
  y.in <- y[[LL]]

  for(j in J:1) {
    LH <- paste("LH", j, sep="")
    HL <- paste("HL", j, sep="")
    HH <- paste("HH", j, sep="")
    
    m <- dim(y.in)[1]
    storage.mode(m) <- "integer"
    n <- dim(y.in)[2]
    storage.mode(n) <- "integer"
    x <- matrix(0, m, n)
    storage.mode(x) <- "double"

    out <- .C("two_D_imodwt", as.double(y.in), as.double(y[[LH]]),
              as.double(y[[HL]]), as.double(y[[HH]]), m, n, j, L,
              h, g, "Y"=x, PACKAGE="waveslim")
    y.in <- out$Y
  }
  zapsmall(y.in)
}

###########################################################################
###########################################################################
###########################################################################

plot.dwt.2d <- function(x, cex.axis=1, plot=TRUE, ...)
{
  J <- attributes(x)$J
  X <- x[[paste("LL", J, sep="")]]
  for(j in J:1) {
    x.names <- sapply(c("LH","HL","HH"), paste, j, sep="")
    X <- rbind(cbind(X, x[[x.names[2]]]),
               cbind(x[[x.names[1]]], x[[x.names[3]]]))
  }
  M <- dim(X)[1]; N <- dim(X)[2]
  if(plot) {
    image(1:M, 1:N, X, col=rainbow(128), axes=FALSE, xlab="", ylab="", ...)
    x.label <- NULL
    lines(c(0,N,N,0,0) + 0.5, c(0,0,M,M,0) + 0.5)
    for(j in J:1) {
      lines(c(M/2^j,M/2^j) + 0.5, 2*c(0,N/2^j) + 0.5)
      lines(2*c(0,M/2^j) + 0.5, c(N/2^j,N/2^j) + 0.5)
    }
    at <- c((3*N+2)/2^(1:J+1),(N+2)/2^(J+1))
    labs <- c(paste("H",1:J,sep=""), paste("L",J,sep=""))
    axis(side=1, at=at, labels=labs, tick=FALSE, cex.axis=cex.axis)
    axis(side=2, at=at, labels=labs, tick=FALSE, cex.axis=cex.axis)
  }
  else
    return(X)
  invisible()
}

###########################################################################
###########################################################################
###########################################################################

denoise.dwt.2d <- function(x, wf = "la8", J = 4, method = "universal", 
                           H = 0.5, noise.dir = 3, rule = "hard")
{
  soft <- function(x, delta) sign(x) * pmax(abs(x) - delta, 0)
  hard <- function(x, delta) ifelse(abs(x) > delta, x, 0)

  n <- length(x)
  x.dwt <- dwt.2d(x, wf, J)
  if(noise.dir == 3)
    sigma.mad <- list(HH = mad(x.dwt$HH1), HL = mad(x.dwt$HL1), 
  		      LH = mad(x.dwt$LH1))
  else {
    noise <- x.dwt$jj
    sigma.mad <- list(HH = mad(noise), HL = mad(noise), LH = mad(noise))
  }
    thresh <- list(HH = rep(sqrt(2 * sigma.mad$HH^2 * log(n)), J), 
		   HL = rep(sqrt(2 * sigma.mad$HL^2 * log(n)), J),
		   LH = rep(sqrt(2 * sigma.mad$LH^2 * log(n)), J))

  if(method == "long-memory")
    thresh <- lapply(thresh, function(x,J,H) 2^(0:(J-1)*(H-1/2))*x, J=J, H=H)
  for(j in 1:J) {
    jj <- paste("HL", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$HL[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$HL[j])
    jj <- paste("LH", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$LH[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$LH[j])
    jj <- paste("HH", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$HH[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$HH[j])
  }
  idwt.2d(x.dwt)
}

###########################################################################
###########################################################################
###########################################################################

denoise.modwt.2d <- function(x, wf = "la8", J = 4, method = "universal", 
  H = 0.5, rule = "hard")
{
  soft <- function(x, delta) sign(x) * pmax(abs(x) - delta, 0)
  hard <- function(x, delta) ifelse(abs(x) > delta, x, 0)
  n <- length(x)
  x.modwt <- modwt.2d(x, wf, J)
  sigma.mad <- list(HH = sqrt(2) * mad(x.modwt$HH1),
		    HL = sqrt(2) * mad(x.modwt$HL1),
		    LH = sqrt(2) * mad(x.modwt$LH1))
    thresh <- list(HH = rep(sqrt(2 * sigma.mad$HH^2 * log(n))/2^(1:J), J), 
		   HL = rep(sqrt(2 * sigma.mad$HL^2 * log(n))/2^(1:J), J), 
		   LH = rep(sqrt(2 * sigma.mad$LH^2 * log(n))/2^(1:J), J))
  if(method == "long-memory")
    thresh <- lapply(thresh, function(x,J,H) 2^(0:(J-1)*(H-1/2))*x, J=J, H=H)
  for(j in 1:J) {
    jj <- paste("HL", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$HL[j])
    else 
      x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$HL[j])
    jj <- paste("LH", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$LH[j])
    else 
      x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$LH[j])
    jj <- paste("HH", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$HH[j])
    else 
     x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$HH[j])
  }
  imodwt.2d(x.modwt)
}

###########################################################################
###########################################################################
###########################################################################

dwpt.2d <- function(x, wf="la8", J=4, boundary="periodic")
{
  ## x <- xbox
  ## Define image dimensions (assign mode for C) and perform simple
  ## diagnostics.
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"
  if(log(m, 2) != trunc(log(m, 2)) | log(n, 2) != trunc(log(n, 2)))
    stop("One dimension is not a power of 2")
  if(2^J > m | 2^J > n)
    stop("Wavelet transform exceeds sample size in one dimension of DWPT")

  ## Extract wavelet and scaling filter coefficients, along with filter
  ## length, from the filter name provided.
  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  ## Create names for wavelet packet nodes (quad-tree structure).
  N <- sum(4^(1:J))
  level <- rep(1:J, 4^(1:J))
  x.wpt <- vector("list", N)
  c1 <- rep(1:J, 2^(1:J))
  c2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  cry <- paste("w", c1, ".", c2, sep="")
  x.wpt.names <- NULL
  for(j in 1:J) {
    xx <- matrix(cry[c1 == j], 2^j, 2^j)
    yy <- matrix(cry[c1 == j], 2^j, 2^j, byrow=TRUE)
    x.wpt.names <- c(x.wpt.names, as.matrix(paste(xx, "-", yy, sep="")))
  }
  names(x.wpt) <- x.wpt.names
  rm(j,xx,yy,c1,c2,cry)
  ## Define initial zero matrix to store wavelet sub-images.
  z <- matrix(0, m/2, n/2)
  storage.mode(z) <- "double"

  ## Implement the 2D DWPT in a nested loop structure.
  for(j in 1:J) {
    ## cat("j =", j, fill=TRUE)
    for(k in 0:(4^j/4-1)) {
      if(j > 1) {
        ## if j > 1, grab wavelet coefficient image and also its name.
        index <- min((1:N)[level == j-1]) + k
        parent <- x.wpt.names[index]
        ## cat("parent =", parent, fill=TRUE)
        x <- x.wpt[[parent]]
        tmp <- unlist(strsplit(parent, "\\-"))
      }
      else
        tmp <- c("w0.0", "w0.0")
      ## Deconstruct name into nodes for the x and y dimensions.
      node <- unlist(strsplit(tmp, "\\."))
      node <- as.integer(node[-c(1,3)])
      ## Preliminary assignments in order to keep wavelet coefficient
      ## sub-images in sequency order.
      if(node[1] %% 2 == 0) {
        Xlow <- paste("w", j, ".", 2 * node[1], sep="")
        Xhigh <- paste("w", j, ".", 2 * node[1] + 1, sep="")
      }
      else {
        Xlow <- paste("w", j, ".", 2 * node[1] + 1, sep="")
        Xhigh <- paste("w", j, ".", 2 * node[1], sep="")
      }
      if(node[2] %% 2 == 0) {
        Ylow <- paste("w", j, ".", 2 * node[2], sep="")
        Yhigh <- paste("w", j, ".", 2 * node[2] + 1, sep="")
      }
      else {
        Ylow <- paste("w", j, ".", 2 * node[2] + 1, sep="")
        Yhigh <- paste("w", j, ".", 2 * node[2], sep="")
      }
      ## Create names for the new wavelet coefficient images.
      LL <- paste(Xlow, "-", Ylow, sep="")
      LH <- paste(Xlow, "-", Yhigh, sep="")
      HL <- paste(Xhigh, "-", Ylow, sep="")
      HH <- paste(Xhigh, "-", Yhigh, sep="")
      ## cat(matrix(c(LH,LL,HH,HL), 2, 2), fill=TRUE)
      ## Perform the DWPT
      out <- .C("two_D_dwt", "Image"=as.double(x), "Rows"=m, "Cols"=n, 
                "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z, "LH"=z,
                "HL"=z, "HH"=z, PACKAGE="waveslim")[7:10]
      ## Pass wavelet coefficient images into the DWPT object.
      x.wpt[[LL]] <- out[["LL"]]
      x.wpt[[LH]] <- out[["LH"]]
      x.wpt[[HL]] <- out[["HL"]]
      x.wpt[[HH]] <- out[["HH"]]
    }
    ## Redefine zero matrix to its new (decimated) size.
    m <- dim(out[["LL"]])[1]
    storage.mode(m) <- "integer"
    n <- dim(out[["LL"]])[2]
    storage.mode(n) <- "integer"
    z <- matrix(0, m/2, n/2)
    storage.mode(z) <- "double"
  }
  attr(x.wpt, "J") <- J
  attr(x.wpt, "wavelet") <- wf
  attr(x.wpt, "boundary") <- boundary
  return(x.wpt)
}

###########################################################################
###########################################################################
###########################################################################

idwpt.2d <- function(y, y.basis)
{
  ## Error checking
  if(length(y) != length(y.basis))
    stop("DWPT object and basis selection must be the same length")
  ## Number of wavelet scales
  J <- attributes(y)$J
  ## Define wavelet/scaling filter coefficients and length
  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"
  ## Nested for loops
  names(y.basis) <- names(y)
  for(j in J:1) {
    for(nx in seq(0, 2^j - 1, by = 2)) {
      for(ny in seq(0, 2^j - 1, by = 2)) {
        ## Name the four wavelet coefficients sub-images
        LL <- paste("w", j, ".", nx, "-", "w", j, ".", ny, sep="")
        LH <- paste("w", j, ".", nx, "-", "w", j, ".", ny+1, sep="")
        HL <- paste("w", j, ".", nx+1, "-", "w", j, ".", ny, sep="")
        HH <- paste("w", j, ".", nx+1, "-", "w", j, ".", ny+1, sep="")
        if(any(y.basis[LL], y.basis[LH], y.basis[HL], y.basis[HH])) {
          m <- nrow(y[[LL]])
          storage.mode(m) <- "integer"
          n <- ncol(y[[LL]])
          storage.mode(n) <- "integer"
          XX <- matrix(0, 2*m, 2*n)
          storage.mode(XX) <- "double"
          ## parent indices to construct string
          pnx <- floor(nx / 2)
          pny <- floor(ny / 2)
          if((pnx %% 2 != 0) & (pny %% 2 != 0))
            ## Upper right-hand corner
            out <- .C("two_D_idwt", as.double(y[[HH]]),
                      as.double(y[[HL]]), as.double(y[[LH]]),
                      as.double(y[[LL]]), m, n, L, h, g, "Y"=XX,
                      PACKAGE="waveslim")$Y
          else {
            ## Upper left-hand corner
            if((pnx %% 2 == 0) & (pny %% 2 != 0))
              out <- .C("two_D_idwt", as.double(y[[LH]]),
                        as.double(y[[LL]]), as.double(y[[HH]]),
                        as.double(y[[HL]]), m, n, L, h, g, "Y"=XX,
                        PACKAGE="waveslim")$Y
            else {
              ## Lower right-hand corner
              if((pnx %% 2 != 0) & (pny %% 2 == 0))
                out <- .C("two_D_idwt", as.double(y[[HL]]),
                          as.double(y[[HH]]), as.double(y[[LL]]),
                          as.double(y[[LH]]), m, n, L, h, g, "Y"=XX,
                          PACKAGE="waveslim")$Y
              else {
                ## Lower left-hand corner
                if((pnx %% 2 == 0) & (pny %% 2 == 0))
                  out <- .C("two_D_idwt", as.double(y[[LL]]),
                            as.double(y[[LH]]), as.double(y[[HL]]),
                            as.double(y[[HH]]), m, n, L, h, g, "Y"=XX,
                            PACKAGE="waveslim")$Y
                else
                  stop("Ouch!")
              }
            }
          }
          if(j > 1) {
            pname <- paste("w", j-1, ".", pnx, "-", "w", j-1, ".", pny, sep="")
            y[[pname]] <- out
            y.basis[pname] <- 1
          }
        }
      }
    }
  }
  return(out)
}

