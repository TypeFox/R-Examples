
### Default grid sizes

default.gridsize <- function(d)
{
  if (d==1)      gridsize <- 401
  else if (d==2) gridsize <- rep(151,d)
  else if (d==3) gridsize <- rep(51, d)
  else if (d==4) gridsize <- rep(21, d)
  else gridsize <- NA
  
  return(gridsize)
}

default.bgridsize <- function(d)
{
  if (d==1)      gridsize <- 401
  else if (d==2) gridsize <- rep(151,d)
  else if (d==3) gridsize <- rep(31, d)
  else if (d==4) gridsize <- rep(11, d)
  else  gridsize <- NA
  
  return(gridsize)
}

default.bflag <- function(d, n)
{
  if (d==1) thr <- 1
  else if (d==2) thr <- 500
  else if (d>2) thr <- 1000
  
  return(n>thr)
}

########################################################################
## Linear binning
## Courtesy of M Wand 2005
## Extended by T Duong to 3- and 4-dim 2006
## Extended by Gramack & Gramacki to include unconstrained b/w 2015 
########################################################################

binning <- function(x, H, h, bgridsize, xmin, xmax, supp=3.7, w, gridtype="linear")
{
  x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)
  if (missing(w)) w <- rep(1,n)
  if (missing(h)) h <- rep(0,d)
  if (!missing(H)) h <- sqrt(diag(H))

  if (missing(bgridsize)) bgridsize <- default.gridsize(d)
  if (!(missing(xmin) & missing(xmax)))
  {
    range.x <- list()
    for (i in 1:d)
      range.x[[i]] <- c(xmin[i], xmax[i])
  }
  else
  {
    range.x <- list()
    for (i in 1:d)
      range.x[[i]] <- c(min(x[,i]) - supp*h[i], max(x[,i]) + supp*h[i])
  }
  if (!(missing(xmax))) {if (any(sweep(x, 2, FUN=">", xmax))) warning("Points in x greater than xmax don't contribute to binning grid counts.")}
  if (!(missing(xmin))) {if (any(sweep(x, 2, FUN="<", xmin))) warning("Points in x less than xmin don't contribute to binning grid counts.")}
      
  a <- unlist(lapply(range.x,min))
  b <- unlist(lapply(range.x,max))
  
  if (missing(gridtype)) gridtype <- rep("linear", d)
  gridtype.vec <- rep("", d)
  
  gpoints <- list()
  for (id in 1:d)
  {
    gridtype1 <- match.arg(gridtype[i], c("linear", "sqrt", "quantile", "log"))
    if (gridtype1=="linear")
      gpoints[[id]] <- seq(a[id],b[id],length=bgridsize[id])  
    else if (gridtype1=="log")
      gpoints[[id]] <- seq(exp(a[id]),exp(b[id]),length=bgridsize[id])
  }
  
  if (d==1) counts <- linbin.ks(x,gpoints[[1]], w=w) 
  if (d==2) counts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]], w=w)
  if (d==3) counts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]], w=w)
  if (d==4) counts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]], w=w)

  bin.counts <- list(counts=counts, eval.points=gpoints, w=w)
  if (d==1) bin.counts <- lapply(bin.counts, unlist)
  return(bin.counts)
}


########################################################################
## Linear binning
########################################################################

linbin.ks <- function(x, gpoints, w)
{
   n <- length(x)
   M <- length(gpoints)
   if (missing(w)) w <- rep(1, n)
   a <- gpoints[1]
   b <- gpoints[M]
   xi <- .C("massdist1d", x1=as.double(x[,1]), n=as.integer(n), a1=as.double(a), b1=as.double(b), M1=as.integer(M), weight=as.double(w), est=double(M), PACKAGE="ks")$est

   return(xi)
}

linbin2D.ks <- function(x, gpoints1, gpoints2, w)
{
   n <- nrow(x)
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   if (missing(w)) w <- rep(1, n)

   ## binning for interior points
   out <- .C("massdist2d", x1=as.double(x[,1]), x2=as.double(x[,2]), n=as.integer(n), a1=as.double(a1), a2=as.double(a2), b1=as.double(b1), b2=as.double(b2), M1=as.integer(M1), M2=as.integer(M2), weight=as.double(w), est=double(M1*M2), PACKAGE="ks")
   xi <- matrix(out$est, nrow=M1, ncol=M2)   
   return(xi)
}

linbin3D.ks <- function(x, gpoints1, gpoints2, gpoints3, w)
{
   n <- nrow(x)
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   if (missing(w)) w <- rep(1, n)

   ## binning for interior points
   out <- .C("massdist3d", x1=as.double(x[,1]), x2=as.double(x[,2]), x3=as.double(x[,3]), n=as.integer(n), a1=as.double(a1), a2=as.double(a2), a3=as.double(a3), b1=as.double(b1), b2=as.double(b2), b3=as.double(b3), M1=as.integer(M1), M2=as.integer(M2), M3=as.integer(M3), weight=as.double(w), est=double(M1*M2*M3), PACKAGE="ks")
   xi <- array(out$est, dim=c(M1,M2,M3))

   return(xi)
}

linbin4D.ks <- function(x, gpoints1, gpoints2, gpoints3, gpoints4, w)
{
   n <- nrow(x)
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3)
   M4 <- length(gpoints4)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   a4 <- gpoints4[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   b4 <- gpoints4[M4]
   if (missing(w)) w <- rep(1, n)

   ## binning for interior points
   out <- .C("massdist4d", x1=as.double(x[,1]), x2=as.double(x[,2]), x3=as.double(x[,3]), x4=as.double(x[,4]), n=as.integer(n), a1=as.double(a1), a2=as.double(a2), a3=as.double(a3), a4=as.double(a4), b1=as.double(b1), b2=as.double(b2), b3=as.double(b3), b4=as.double(b4), M1=as.integer(M1), M2=as.integer(M2), M3=as.integer(M3), M4=as.integer(M4), weight=as.double(w), est=double(M1*M2*M3*M4), PACKAGE="ks")
   xi <- array(out$est, dim=c(M1,M2,M3,M4))

   return(xi)
}

########################################################################
## Discrete convolution
########################################################################

symconv.1d <- function(keval, gcounts)
{
    M <- length(gcounts)
    N <- length(keval) 
    L <- (length(keval)+1)/2
  
    ## Smallest powers of 2 >= M + N
    P <- 2^(ceiling(log2(M + N)))
    
    ## Zero-padded version of keval an gcounts
    keval.zeropad <- rep(0, P)
    gcounts.zeropad <- rep(0, P)
    keval.zeropad[1:N] <- keval
    gcounts.zeropad[L:(L+M-1)] <- gcounts

    ## FFTs
    K <- fft(keval.zeropad)
    C <- fft(gcounts.zeropad)

    ## Invert element-wise product of FFTs and truncate and normalise it
    symconv.val <- Re(fft(K*C, inverse=TRUE)/P)[N:(N+M-1)]
    
    return(symconv.val)
}

symconv.nd <- function(keval, gcounts, d)
{
    M <- dim(gcounts)
    N <- dim(keval) 
    L <- (dim(keval)+1)/2
  
    ## Smallest powers of 2 >= M + N
    P <- 2^(ceiling(log2(M + N)))
    
    ## Zero-padded version of keval and gcounts
    keval.zeropad <- array(0, dim=P)
    gcounts.zeropad <- array(0, dim=P)
    if (d==2)
    {
        keval.zeropad[1:N[1], 1:N[2]] <- keval
        gcounts.zeropad[L[1]:(L[1]+M[1]-1), L[2]:(L[2]+M[2]-1)] <- gcounts
    }
    else if (d==3)
    {
        keval.zeropad[1:N[1], 1:N[2], 1:N[3]] <- keval
        gcounts.zeropad[L[1]:(L[1]+M[1]-1), L[2]:(L[2]+M[2]-1), L[3]:(L[3]+M[3]-1)] <- gcounts
    }
    else if (d==4)
    {
        keval.zeropad[1:N[1], 1:N[2], 1:N[3], 1:N[4]] <- keval
        gcounts.zeropad[L[1]:(L[1]+M[1]-1), L[2]:(L[2]+M[2]-1), L[3]:(L[3]+M[3]-1), L[4]:(L[4]+M[4]-1) ] <- gcounts
    }
    
    ## FFTs
    K <- fft(keval.zeropad)
    C <- fft(gcounts.zeropad)

    ## Invert element-wise product of FFTs and truncate and normalise it
    symconv.val <- Re(fft(K*C, inverse=TRUE)/prod(P))
    if (d==2) symconv.val <- symconv.val[N[1]:(N[1]+M[1]-1), N[2]:(N[2]+M[2]-1)]
    else if (d==3) symconv.val <- symconv.val[N[1]:(N[1]+M[1]-1), N[2]:(N[2]+M[2]-1), N[3]:(N[3]+M[3]-1)]
    else if (d==4) symconv.val <- symconv.val[N[1]:(N[1]+M[1]-1), N[2]:(N[2]+M[2]-1), N[3]:(N[3]+M[3]-1), N[4]:(N[4]+M[4]-1)]
    
    return(symconv.val)
}

