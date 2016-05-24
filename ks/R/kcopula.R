#############################################################################
## Kernel copula and copula density estimators
#############################################################################


## empirical pseudo-uniform transformation
## taken from pobs() function in copula package
pseudo.unif.empirical <- function (x, y) ##, na.last="keep", ties.method="average")
{
  if (missing(y)) y <- x
  if (is.vector(y)) y <- matrix(y, nrow=1)
  d <- ncol(x)
  u <- matrix(0, ncol=d, nrow=nrow(x))
  for (i in 1:d)
  {
    ecdf.fun <- ecdf(x=x[,i])
    u[,i] <- ecdf.fun(y[,i])
  }
  return(u)
  ##return(apply(x, 2, rank, na.last=na.last, ties.method=ties.method)/(nrow(x) + 1))
}

## kernel pseudo-uniform transformation
pseudo.unif.kernel <- function(x, y, hs, binned=TRUE)
{
  if (missing(y)) y <- x
  if (is.vector(y)) y <- matrix(y, nrow=1)
  d <- ncol(x)
  u <- list()
  ##u.eval.points <- list()
  for (i in 1:d)
  {
    u[[i]] <- kcde(x=x[,i], h=hs[i], eval.points=y[,i], binned=binned)
  }
  u2 <- sapply(u, getElement, "estimate")

  return(u2)
}

## Modified boundary beta kernel - first form (Chen, 1999)
dbeta.kernel <- function(x, eval.x, h)
{
  return (dbeta(x=eval.x, shape1=x/h^2+1, shape2=(1-x)/h^2+1))
}

## Modified boundary beta kernel - second form (Chen, 1999)
dbeta.kernel2 <- function(x, eval.x, h)
{
  rhox <- function(y, hy) {if (y==0) return (1) else return(2*hy^4 + 5/2 - sqrt(4*hy^8 + 6*hy^4 + 9/4 - y^2 -y/hy^2))}
  ind <- cut(x, c(0, 2*h^2, 1-2*h^2, 1), labels=FALSE, include.lowest=TRUE)
  dbk <- rep(0, length(eval.x))

  if (ind==1) {shape1 <- rhox(x,hy=h); shape2 <- (1-x)/h^2}
  else if (ind==2) {shape1 <- x/h^2;   shape2 <- (1-x)/h^2}
  else if (ind==3) {shape1 <- x/h^2;   shape2 <- rhox(1-x, hy=h)}

  return(dbeta(eval.x, shape1=shape1, shape2=shape2))
}

## Modified multivariate boundary beta product kernel

dmvbeta.prod.kernel2 <- function(x, eval.x, hs)
{
  d <- length(hs)
  db <- vector("list", d)
  for (i in 1:d) db[[i]] <- 0
  
  for (i in 1:d)
      db[[i]] <- dbeta.kernel2(x=x[i], eval.x=eval.x[[i]], h=hs[i])

  db <- expand.grid(db)
  db <- apply(db, 1, prod)
  return(db)
}

## Modified multivariate boundary beta spherically symmetric kernel

dmvbeta.symm.kernel2 <- function(x, eval.x, H)
{
  d <- ncol(H)
  eval.y <- sqrt(apply(eval.x^2, 1, sum))/sqrt(d)
  y <- sqrt(sum(x^2))/sqrt(d)
  return(dbeta.kernel2(x=y, eval.x=eval.y, h=sqrt(tr(H)))/d)
}


#############################################################################
## Kernel copula estimator
#############################################################################

kcopula <- function(x, H, hs, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, w, verbose=FALSE, marginal="kernel")
{
  d <- ncol(x)
  n <- nrow(x)

  if (missing(H)) H <- Hpi.kcde(x=x, binned=default.bflag(d=d,n=n))
  Fhat <- kcde(x=x, H=H, gridsize=gridsize, binned=binned, bgridsize=bgridsize, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, w=w, verbose=verbose, tail.flag="lower.tail")

  xlims <- sapply(Fhat$eval.points, range)
  xlims[1,] <- xlims[1,] - 0.1*abs(apply(xlims, 2, diff))
  xlims[2,] <- xlims[2,] + 0.1*abs(apply(xlims, 2, diff))

  ## generate pseudo-uniform values
  marginal1 <- match.arg(marginal, c("kernel", "empirical")) 
  if (missing(hs))
  {
    hs <- rep(0, d)
    for (i in 1:d) hs[i] <- hpi.kcde(x=x[,i], binned=TRUE)
  }
  if (marginal1=="kernel")
  {  
    ## kernel pseudo-uniform
    u <- list()
    u.eval.points <- list()
    for (i in 1:d)
    {
      u.eval.points[[i]] <- kcde(x=x[,i], h=hs[i], eval.points=Fhat$eval.points[[i]], xmin=xlims[1,i], xmax=xlims[2,i], binned=TRUE)
    }
    y <- pseudo.unif.kernel(x=x, y=x, hs=hs, binned=TRUE)
    ep <- lapply(u.eval.points, getElement, "estimate")
  }
  else if (marginal1=="empirical")
  {
    ## empirical pseudo-uniform
    y <- pseudo.unif.empirical(x=x, y=x)
    ep <- numeric()
    for (i in 1:d)
    {
      f <- ecdf(x=x[,i])
      ep <- c(ep, list(f(Fhat$eval.points[[i]])))
    }
  }
 
  Chat <- Fhat
  Chat$x <- y ## sapply(u, getElement, "estimate")
  Chat$x.orig <- x
  Chat$eval.points <- ep ##lapply(u.eval.points, getElement, "estimate")
  Chat$hs <- hs

  ## loess smoothing on a uniform grid
  if (d==2)
  {
    ## select smaller grid for d==2 for memory usage reasons
    subselect <- round(cbind(seq(1,length(Chat$eval.points[[1]]), length=101), seq(1,length(Chat$eval.points[[2]]), length=101)),0)
    eval.points.df <- data.frame(expand.grid(Chat$eval.points[[1]][subselect[,1]], Chat$eval.points[[2]][subselect[,2]]))
    names(eval.points.df) <- paste("x", 1:ncol(eval.points.df), sep="")
    eval.points.df <- data.frame(estimate=as.numeric(Chat$estimate[subselect[,1], subselect[,2]]), eval.points.df)
  }
  else if (d==3)
  {
    ## select smaller grid for d==3 for memory usage reasons
    subselect <- round(cbind(seq(1,length(Chat$eval.points[[1]]), length=21), seq(1,length(Chat$eval.points[[2]]), length=21), seq(1,length(Chat$eval.points[[3]]), length=21)),0)
    eval.points.df <- data.frame(expand.grid(Chat$eval.points[[1]][subselect[,1]], Chat$eval.points[[2]][subselect[,2]], Chat$eval.points[[3]][subselect[,3]]))
    names(eval.points.df) <- paste("x", 1:ncol(eval.points.df), sep="")
    eval.points.df <- data.frame(estimate=as.numeric(Chat$estimate[subselect[,1], subselect[,2], subselect[,3]]), eval.points.df)
  }
  
  if (d==2) Chat.loess <- loess(estimate ~ x1+x2, data=eval.points.df, span=0.1)
  else if (d==3)  Chat.loess <- loess(estimate ~ x1+x2+x3, data=eval.points.df, span=0.1)

  u.eval.points.regular <- list()
  for (i in 1:d)
    u.eval.points.regular[[i]] <- seq(0,1,length=length(Chat$eval.points[[i]]))
  u.eval.points.regular.df <- data.frame(expand.grid(u.eval.points.regular))
  names(u.eval.points.regular.df) <- paste("x", 1:ncol(u.eval.points.regular.df), sep="")
  
  Chat.smoothed <- Chat
  Chat.smoothed$eval.points <- u.eval.points.regular
  Chat.smoothed$estimate <- array(predict(Chat.loess, newdata=u.eval.points.regular.df), dim=dim(Chat$estimate))

  ## interpolate NA boundary values from loess smoothing
  gsdim <- dim(Chat$estimate)
  if (d==2)
  {
    Chat.smoothed$estimate[1,] <- 0
    Chat.smoothed$estimate[gsdim[1],] <- Chat.smoothed$estimate[gsdim[1]-1,]*1.001
    Chat.smoothed$estimate[,1] <-0
    Chat.smoothed$estimate[,gsdim[2]] <- Chat.smoothed$estimate[,gsdim[2]-1]*1.001
    Chat.smoothed$estimate[gsdim[1],gsdim[2]] <- 1
    Chat.smoothed$estimate[Chat.smoothed$estimate>1] <- 1
  }
  else if (d==3)
  {
    Chat.smoothed$estimate[,,1] <- 0
    for (k in 2:(gsdim[3]-1))
    {
       Chat.smoothed$estimate[,,k][1,] <- 0
       Chat.smoothed$estimate[,,k][gsdim[1],] <- Chat.smoothed$estimate[,,k][gsdim[1]-1,]*1.001
       Chat.smoothed$estimate[,,k][,1] <-0
       Chat.smoothed$estimate[,,k][,gsdim[2]] <- Chat.smoothed$estimate[,,k][,gsdim[2]-1]*1.001
     }
    Chat.smoothed$estimate[,,gsdim[3]] <- Chat.smoothed$estimate[,,gsdim[3]-1]*1.001
    Chat.smoothed$estimate[gsdim[1],gsdim[2],gsdim[3]] <- 1
    Chat.smoothed$estimate[Chat.smoothed$estimate>1] <- 1
  }
  
  Chat <- Chat.smoothed
  Chat$marginal <- marginal1
  class(Chat) <- "kcopula"
  
  return(Chat)
}


#############################################################################
## Kernel copula density estimator
#############################################################################

kcopula.de <- function(x, H, Hfun, hs, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, w, verbose=FALSE, compute.cont=FALSE, approx.cont=TRUE, boundary.supp, marginal="kernel", Hfun.pilot="dscalar")
{
  d <- ncol(x)
  n <- nrow(x)

  warning("kcopula.de is not as reliable as kdecop in the kdecopula package.")
  
  if (!missing(w))
    if (!(identical(all.equal(sum(w), n), TRUE)))
    {
      warning("Weights don't sum to sample size - they have been scaled accordingly\n")
      w <- w*n/sum(w)
    }
  if (missing(w)) w <- rep(1,n)
  if (missing(hs))
  {
    hs <- rep(0, d)
    for (i in 1:d) hs[i] <- hpi.kcde(x=x[,i], binned=TRUE)
  }

  ## generate pseudo-uniform values
  marginal1 <- match.arg(marginal, c("kernel", "empirical")) 
  if (marginal1=="kernel")
  {  
    ## kernel pseudo-uniform
    y <- pseudo.unif.kernel(x=x, y=x, hs=hs, binned=TRUE)
  }
  else if (marginal1=="empirical")
  {
    ## empirical pseudo-uniform
    y <- pseudo.unif.empirical(x=x, y=x)  
  }

  if (missing(gridsize)) gridsize <- default.gridsize(d)
  if (missing(bgridsize)) bgridsize <- default.gridsize(d)
  if (missing(Hfun)) Hfun <- Hpi

  if (missing(H))
  {  
    if (d==2)
      H <- do.call(Hfun, list(x=y, binned=default.bflag(d=d, n=nrow(y)), bgridsize=bgridsize, pilot=Hfun.pilot, verbose=verbose))
    else if (d==3)
      H <- do.call(Hfun, list(x=y, binned=default.bflag(d=d, n=nrow(y)), bgridsize=rep(31,d), pilot=Hfun.pilot, verbose=verbose)) 
  }
  
  #######################################################################
  ## Need to calibrate boundary.supp value for d=3
  #######################################################################
  if (missing(boundary.supp)) { if (d==2) boundary.supp <- 1; if (d==3) boundary.supp <- 0.5 }

  ## kernel copula density is boundary kernel estimator 
  if (d==2 | d==3) 
    chat <- kde.boundary(x=y, H=H, gridsize=gridsize, supp=supp, xmin=rep(0,d), xmax=rep(1,d), gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned, verbose=verbose, bgridsize=bgridsize)
 else
   stop("kcopula.de requires 2-d or 3-d data.")
   
  ## normalise KDE to integrate to 1 
  chat$estimate <- chat$estimate/sum(chat$estimate*apply(sapply(chat$eval.points, diff), 1, prod)[1])
  chat$names <- parse.name(x) 
  chat$x.orig <- x
  chat$hs <- hs

  ## compute prob contour levels
  if (compute.cont & missing(eval.points))
    chat$cont <- contourLevels(chat, cont=1:99, approx=approx.cont)

  chat$marginal <- marginal1
  class(chat) <- "kcopula.de"
  return(chat)
}



#############################################################################
## Boundary kernel estimator
#############################################################################

## indicator function for boundary region of [0,1] i.e. [0,h] + [1-h, h]

boundary.ind <- function(x, h, xmin, xmax, boundary.supp=1)
{
  if (is.vector(x)){ x <- matrix(x, ncol=1) }
  n <- nrow(x)
  d <- ncol(x)
  
  ## indicator for closeness to boundary of [0,1]^d
  bound.ind <- matrix(NA, nrow=n, ncol=d)
  for (j in 1:d) bound.ind[,j] <- (abs(x[,j]) <= boundary.supp*h[j]) | (abs(1-x[,j]) <= boundary.supp*h[j])
  bound.ind <- apply(bound.ind, 1, any)
  
  return(bound.ind)
}



## boundary kernel estimator using beta bounday kernels (2nd form)

kde.boundary <- function(x, H, h, hb, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, w, compute.cont=FALSE, approx.cont=TRUE, boundary.supp=1, verbose=FALSE)
{
  if (is.vector(x))
  {
    if (missing(H)) {d <- 1; n <- length(x)}
    else
    {
      if (is.vector(H)) { d <- 1; n <- length(x)}
      else {x <- matrix(x, nrow=1); d <- ncol(x); n <- nrow(x)}
    }
  }
  else {d <- ncol(x); n <- nrow(x)}

  if (!missing(w))
    if (!(identical(all.equal(sum(w), n), TRUE)))
    {
      warning("Weights don't sum to sample size - they have been scaled accordingly\n")
      w <- w*n/sum(w)
    }

  if (missing(w)) w <- rep(1,n)
  if (d==1) { if (missing(h)) h <- hpi(x=x, binned=default.bflag(d=d,n=n), bgridsize=bgridsize) }
  if (missing(H) & d>1)  H <- Hpi(x=x, binned=default.bflag(d=d,n=n), bgridsize=bgridsize)
      
  
  ## compute exact (non-binned) estimator
  if (missing(gridsize)) gridsize <- default.gridsize(d)

  ## 1-dimensional    
  if (d==1)
  {
    if (missing(eval.points))
    {
      fhat <- kde.boundary.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned)
    }
    else
      stop("Not yet implemented.")  ##fhat <- kde.points.1d(x=x, h=h, eval.points=eval.points, positive=positive, adj.positive=adj.positive, w=w)
  }
  ## multi-dimensional
  else
  {  
    if (is.data.frame(x)) x <- as.matrix(x)
    
    if (missing(eval.points))
    {
      if (d==2)
        fhat <- kde.boundary.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned, verbose=verbose)
      else if (d == 3)
        fhat <- kde.boundary.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned, verbose=verbose) 
      else 
        stop("Need to specify eval.points for more than 3 dimensions")
    }
    else
      stop("Not yet implemented.") ##fhat <- kde.points(x=x, H=H, eval.points=eval.points, w=w)     
  }

  fhat$binned <- binned
  fhat$names <- parse.name(x)  ## add variable names
  fhat$w <- w
  class(fhat) <- "kde"
  
  ## compute prob contour levels
  if (compute.cont & missing(eval.points))
    fhat$cont <- contourLevels(fhat, cont=1:99, approx=approx.cont)

  return(fhat)
 }



kde.boundary.grid.1d <- function(x, h, hb, gridsize, supp=3.7, xmin, xmax, gridtype, w, boundary.supp=0.5, binned=FALSE)
{
  if (missing(xmin)) xmin <- min(x)
  if (missing(xmax)) xmax <- max(x)
  if (missing(gridtype)) gridtype <- "linear"

  ## transform x into [0,1]
  x.star <- (x-xmin)/(xmax-xmin) 
  h.star <- h/(xmax-xmin)
  n <- length(x)

  gridtype1 <- match.arg(gridtype, c("linear", "sqrt")) 
  if (gridtype1=="linear")
    eval.x <- seq(0, 1, length=gridsize)
  else if (gridtype1=="sqrt")
  {
    eval.x.temp <- seq(0, 1, length=gridsize)
    eval.x <- sign(eval.x.temp) * eval.x.temp^2
  }
  gridtype.vec <- gridtype1

  ## indicator for closeness to boundary of [0,1]
  bound.ind <- boundary.ind(x=x.star, h=h.star, boundary.supp=boundary.supp) 

  n1 <- sum(!bound.ind)
  ## interior points - use normal kernel
  ## binned estimation only in the interior
  fhat.grid <- rep(0,length=gridsize)
  if (n1>0)
  {
    if (binned)
    {
      fhat.grid <- n1*kde(x=x.star[!bound.ind], h=h.star, xmin=0, xmax=1, binned=TRUE, bgridsize=gridsize)$estimate
    }
    else
    {
      fhat.grid <- n1*dnorm.mixt(x=eval.x, mus=x.star[!bound.ind], sigmas=rep(h.star, n1), props=w[!bound.ind]/n1)
    }
  }

  ## boundary points - use adjusted beta kernel
  hb.star <- 2*h.star
  for (i in 1:(n-n1))
    fhat.grid <- fhat.grid + dbeta.kernel2(x=x.star[bound.ind][i], eval.x=eval.x, h=hb.star)*w[bound.ind][i]
  fhat.grid <- fhat.grid/n

  ## backtransform
  eval.points <- (xmax-xmin)*eval.x + xmin 
  fhat.grid <- fhat.grid/(xmax-xmin)

  fhat <- list(x=x, eval.points=eval.points, estimate=fhat.grid, h=h, H=h^2, gridtype=gridtype.vec, gridded=TRUE)
  class(fhat) <- "kde"
  
  return(fhat)
}


kde.boundary.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, boundary.supp=1, binned=FALSE, verbose=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)

  if (missing(xmin)) xmin <- apply(x, 2, min) ##- h*supp
  if (missing(xmax)) xmax <- apply(x, 2, max) ## + h*supp
  if (missing(gridtype)) gridtype <- rep("linear", d)

  ## transform x into [0,1]^d
  x.star <- x
  for (j in 1:d) x.star[,j] <- (x[,j]-xmin[j])/(xmax[j]-xmin[j])
  H.star <- diag(1/(xmax-xmin)) %*% H %*% diag(1/(xmax-xmin))
  h.star <- sqrt(diag(H.star))
  
  ## initialise grid 
  if (is.null(gridx))
    gridx <- make.grid.ks(x.star, matrix.sqrt(H.star), tol=supp, gridsize=gridsize, xmin=rep(0,d), xmax=rep(1,d), gridtype=gridtype) 
  suppx <- make.supp(x.star, matrix.sqrt(H.star), tol=supp)

  if (is.null(grid.pts)) grid.pts <- find.gridpts(gridx, suppx)    
  fhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))

  ## indicator for closeness to boundary of [0,1]^d
  bound.ind <- boundary.ind(x=x.star, h=h.star, boundary.supp=boundary.supp)
  n1 <- sum(!bound.ind)

  if (verbose) pb <- txtProgressBar()
  if (binned)
  {
    ## interior points - use normal kernel
    fhat.grid <- n1*kde(x=x.star[!bound.ind,], H=H.star, xmin=rep(0,d), xmax=rep(1,d), binned=TRUE, bgridsize=gridsize)$estimate

    for (i in 1:n)
    {
      if (verbose) setTxtProgressBar(pb, i/n)
      if (bound.ind[i])
      {
        ## compute evaluation points 
        eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
        eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
        eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
        eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
        eval.x.len <- length(eval.x)

        ## convert bandwidth from normal kernel to beta kernel scale
        ## for boundary points 
        hb.star <- 2*h.star
        fhat <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star)
   
        ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
        for (j in 1:length(eval.y))
          fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
            fhat.grid[eval.x.ind, eval.y.ind[j]] + 
              w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
      }
    }
  }
  else
  {  
    for (i in 1:n)
    {
      if (verbose) setTxtProgressBar(pb, i/n)

      ## compute evaluation points 
      eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
      eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
      eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
      eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
      eval.x.len <- length(eval.x)
      
      ## interior points - use normal kernel
      if (!bound.ind[i])
      {
        eval.pts <- permute(list(eval.x, eval.y))
        fhat <- dmvnorm(eval.pts, x.star[i,], H.star)
      }
      else
      {
        ## convert bandwidth from normal kernel to beta kernel scale
        ## for boundary points 
        hb.star <- 2*h.star
        fhat <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star)
        ##fhat <- dmvbeta.symm.kernel2(x=x.star[i,], eval.x=expand.grid(eval.x, eval.y), H=2*H.star)
      }

      ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
      for (j in 1:length(eval.y))
        fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
          fhat.grid[eval.x.ind, eval.y.ind[j]] + 
          w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
    }
  }
  fhat.grid <- fhat.grid/n
  if (verbose) close(pb)
  
  ## back-transform
  gridx1 <- list((xmax[1]-xmin[1])*gridx[[1]] + xmin[1], (xmax[2]-xmin[2])*gridx[[2]] + xmin[2]) 
  fhat.grid <- fhat.grid/prod(xmax-xmin)
  
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, boundary=bound.ind)
  
  return(fhat.list)
}


kde.boundary.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, boundary.supp=0.5, verbose=FALSE, binned=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)

  if (missing(xmin)) xmin <- apply(x, 2, min) 
  if (missing(xmax)) xmax <- apply(x, 2, max) 
  if (missing(gridtype)) gridtype <- rep("linear", d)

  ## transform x into [0,1]^d
  x.star <- x
  for (j in 1:d) x.star[,j] <- (x[,j]-xmin[j])/(xmax[j]-xmin[j])
  H.star <- diag(1/(xmax-xmin)) %*% H %*% diag(1/(xmax-xmin))
  h.star <- sqrt(diag(H.star)) 

  ## initialise grid 
  if (is.null(gridx))
    gridx <- make.grid.ks(x.star, matrix.sqrt(H.star), tol=supp, gridsize=gridsize, xmin=rep(0,d), xmax=rep(1,d), gridtype=gridtype) 
  suppx <- make.supp(x.star, matrix.sqrt(H.star), tol=supp)
  if (is.null(grid.pts)) grid.pts <- find.gridpts(gridx, suppx)
  fhat.grid <- array(0, dim=c(length(gridx[[1]]), length(gridx[[2]]), length(gridx[[3]])))

  ## indicator for closeness to boundary of [0,1]^d
  bound.ind <- boundary.ind(x=x.star, h=h.star, boundary.supp=boundary.supp)
  n1 <- sum(!bound.ind)
  
  if (verbose) pb <- txtProgressBar() 
  if (binned)
  {
    ## interior points - use normal kernel
    fhat.grid <- n1*kde(x=x.star[!bound.ind,], H=H.star, xmin=rep(0,d), xmax=rep(1,d), binned=TRUE, bgridsize=gridsize)$estimate

    for (i in 1:n)
    {
      if (verbose) setTxtProgressBar(pb, i/n)
      
      if (bound.ind[i])
      {
        ## compute evaluation points
        eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
        eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
        eval.z <- gridx[[3]][grid.pts$xmin[i,3]:grid.pts$xmax[i,3]]
        eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
        eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
        eval.z.ind <- c(grid.pts$xmin[i,3]:grid.pts$xmax[i,3])
        eval.x.len <- length(eval.x)
        eval.pts <- permute(list(eval.x, eval.y))
        
        ## convert bandwidth from normal kernel to beta kernel scale
        ## for boundary points 
        hb.star <- 2*h.star
        fhat.xy <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star[1:2])

        ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
        for (k in 1:length(eval.z))
        {
          fhat <- w[i]*cbind(fhat.xy, dbeta.kernel2(x=x.star[i,3], eval.x=eval.z[k], h=hb.star[3]))
          for (j in 1:length(eval.y))
            fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
              fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
        }
      }
    }
  }
  else
  {  
    for (i in 1:n)
    {
      if (verbose) setTxtProgressBar(pb, i/n)

      ## compute evaluation points
      eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
      eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
      eval.z <- gridx[[3]][grid.pts$xmin[i,3]:grid.pts$xmax[i,3]]
      eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
      eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
      eval.z.ind <- c(grid.pts$xmin[i,3]:grid.pts$xmax[i,3])
      eval.x.len <- length(eval.x)
      eval.pts <- permute(list(eval.x, eval.y))
      
      ## interior points - use normal kernel
      if (!bound.ind[i])
      {
        eval.pts <- permute(list(eval.x, eval.y))
        ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
        
        for (k in 1:length(eval.z))
        {
          fhat <- w[i]*dmvnorm(cbind(eval.pts, eval.z[k]), x[i,], H)
          for (j in 1:length(eval.y))
            fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
              fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
        }
      }
      else
      {
        ## convert bandwidth from normal kernel to beta kernel scale
        hb.star <- 2*h.star
        fhat.xy <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star[1:2])
        
        for (k in 1:length(eval.z))
        {
          fhat <- w[i]*cbind(fhat.xy, dbeta.kernel2(x=x.star[i,3], eval.x=eval.z[k], h=hb.star[3]))
          for (j in 1:length(eval.y))
            fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
              fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
        }
        
      }
    }
  }
  fhat.grid <- fhat.grid/n
  if (verbose) close(pb)
  
  ## back-transform
  gridx1 <- list((xmax[1]-xmin[1])*gridx[[1]] + xmin[1], (xmax[2]-xmin[2])*gridx[[2]] + xmin[2], (xmax[3]-xmin[3])*gridx[[3]] + xmin[3])
  fhat.grid <- fhat.grid/prod(xmax-xmin)
  
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, boundary=bound.ind)
  return(fhat.list)
}



plot.kcopula <- function(x, ...)
{
  plot.kcde(x, ...)
}

plot.kcopula.de <- function(x, ...)
{
  plot.kde(x, ...)
}



#############################################################################
## predict methods
#############################################################################

predict.kcopula <- function(object, ..., x, u)
{
  if (missing(u)) {if (object$marginal=="kernel") u <- pseudo.unif.kernel(x=object$x.orig, y=x, hs=object$hs)}
  return(predict.kde(object, ..., x=u))
}

predict.kcopula.de <- function(object, ..., x, u)
{
  if (missing(u)) {if (object$marginal=="kernel") u <- pseudo.unif.kernel(x=object$x.orig, y=x, hs=object$hs)}
 
  return(predict.kde(object, ..., x=u))
}


contourLevels.kcopula.de <- function(x, ...)
{
  x1 <- x; class(x1) <- "kde"  
  return(contourLevels(x=x1, ...))
}




    
