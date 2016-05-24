###############################################################################
### Multivariate kernel density derivative estimate 
###############################################################################

kdde <- function(x, H, h, deriv.order=0, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w, deriv.vec=TRUE, verbose=FALSE)
{
  r <- deriv.order
  
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

  if (missing(h) & d==1) h <- hpi(x=x, nstage=2, binned=TRUE, bgridsize=bgridsize, deriv.order=r)
  if (missing(H) & d>1)
  {
    H <- Hpi(x=x, nstage=2, binned=default.bflag(d=d, n=n), bgridsize=bgridsize, deriv.order=r, verbose=verbose)
  }
  
  ## compute binned estimator
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
    
    ##if (d>1) { if (!identical(diag(diag(H)), H)) warning("binned estimation for non-diagonal bandwidth matrix H can be inaccurate") }
    
    if (positive & is.vector(x))
    {
      y <- log(x)
      fhat <- kdde.binned(x=y, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w)
      fhat$estimate <- fhat$estimate/exp(fhat$eval.points)
      fhat$eval.points <- exp(fhat$eval.points)
      fhat$x <- x
    }
    else
      fhat <- kdde.binned(x=x, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, deriv.vec=deriv.vec, verbose=verbose)
  }
  else
  {
    ## compute exact (non-binned) estimator
    if (missing(gridsize)) gridsize <- default.gridsize(d)
    
    ## 1-dimensional    
    if (d==1)
    {
      if (!missing(H) & !missing(h))
        stop("Both H and h are both specified")
      
      if (missing(h))
        h <- sqrt(H)

      if (missing(eval.points))
        fhat <- kdde.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r)
      else
        fhat <- kdde.points.1d(x=x, h=h, eval.points=eval.points, w=w, deriv.order=r)
    }
    ## multi-dimensional
    else
    {   
      if (is.data.frame(x)) x <- as.matrix(x)
      
      if (missing(eval.points))
      {
        if (d==2) 
          fhat <- kdde.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r, deriv.vec=deriv.vec, verbose=verbose)
        else if (d==3)
          ##stop("Exact kdde not yet implemented for 3 dimensions")
          fhat <- kdde.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r, deriv.vec=deriv.vec, verbose=verbose) 
        else 
          stop("Need to specify eval.points for more than 3 dimensions")
      }
      else
        fhat <- kdde.points(x=x, H=H, eval.points=eval.points, w=w, deriv.order=r, deriv.vec=deriv.vec)
    }

  }
  fhat$binned <- binned
  fhat$names <- parse.name(x)
  class(fhat) <- "kdde"
  return(fhat)
 }



###############################################################################
## Multivariate binned kernel density derivative estimate
###############################################################################

kdde.binned <- function(x, H, h, deriv.order, bgridsize, xmin, xmax, bin.par, w, deriv.vec=TRUE, deriv.index, verbose=FALSE)
{
  r <- deriv.order
  if (length(r)>1) stop("deriv.order should be a non-negative integer")

  ## linear binning
  if (missing(bin.par))
  {
    if (is.vector(x)) {d <- 1; n <- length(w)}
    else {d <- ncol(x); n <- nrow(x)}

    if (missing(w)) w <- rep(1,n)
    
    if (d==1)
      if (missing(H)) { H <- as.matrix(h^2)} 
      else {h <- sqrt(H); H <- as.matrix(H)}

    if (d==1) Hd <- H else Hd <- diag(diag(H))
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
    
    bin.par <- binning(x=x, H=Hd, h=h, bgridsize, xmin, xmax, supp=3.7+max(r), w=w)
  }
  else
  {
    if (!is.list(bin.par$eval.points)) { d <- 1; bgridsize <- length(bin.par$eval.points)}
    else  { d <- length(bin.par$eval.points); bgridsize <- sapply(bin.par$eval.points, length)} 

    w <- bin.par$w
    if (d==1)
      if (missing(H)) H <- as.matrix(h^2)
      else {h <- sqrt(H); H <- as.matrix(H)}
  }

 
  if (d==1)
  {
    fhat <- kdde.binned.1d(h=h, deriv.order=r, bin.par=bin.par)
    eval.points <- fhat$eval.points
    est <- fhat$estimate
  }
  else
  {
    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, only.index=TRUE, deriv.vec=deriv.vec)
    fhat.grid <- kdde.binned.nd(H=H, deriv.order=r, bin.par=bin.par, verbose=verbose, deriv.vec=deriv.vec)
  }

  if (missing(x)) x <- NULL
  
  if (d==1)
  {
    if (r==0) fhat <- list(x=x, eval.points=unlist(eval.points), estimate=est, h=h, H=h^2, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w)
    else
      fhat <- list(x=x, eval.points=unlist(eval.points), estimate=est, h=h, H=h^2, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=r)
  }
  else
  {
    if (r==0)
      fhat <- list(x=x, eval.points=fhat.grid$eval.points, estimate=fhat.grid$estimate[[1]], H=H, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w)
    else
      fhat <- list(x=x, eval.points=fhat.grid$eval.points, estimate=fhat.grid$estimate, H=H, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=ind.mat)
  }
  
  class(fhat) <- "kdde"
  
  return(fhat)
}

kdde.binned.1d <- function(h, deriv.order, bin.par)
{
  r <- deriv.order
  n <- sum(bin.par$counts)
  a <- min(bin.par$eval.points)
  b <- max(bin.par$eval.points)
  M <- length(bin.par$eval.points)
  L <- min(ceiling((4+r)*h*(M-1)/(b-a)), M-1)
  delta <- (b-a)/(M-1)
  N <- 2*L-1
  grid1 <- seq(-(L-1), L-1)
  
  keval <- dnorm.deriv(x=delta*grid1, mu=0, sigma=h, deriv.order=r)/n
  est <- symconv.1d(keval, bin.par$counts)
  ##keval <- dnorm.deriv(x=(b-a)*(0:L)/(M-1), mu=0, sigma=h, deriv.order=r)/n
  ##est <- symconv.ks.old(keval, bin.par$counts, skewflag=(-1)^r)

  return(list(eval.points=bin.par$eval.points, estimate=est))
}

kdde.binned.nd <- function(H, deriv.order, bin.par, verbose=FALSE, deriv.vec=TRUE)
{
  d <- ncol(H)
  r <- deriv.order
  n <- sum(bin.par$counts)
  a <- sapply(bin.par$eval.points, min)
  b <- sapply(bin.par$eval.points, max)
  M <- sapply(bin.par$eval.points, length)
  L <- pmin(ceiling((4+r)*max(sqrt(abs(diag(H))))*(M-1)/(b-a)), M-1)
  delta <- (b-a)/(M-1)
  N <- 2*L-1 

  ##supp <- 3.7 
  ##L <- numeric(2)
  ##delta <- c(0, 0)
  #eig <- eigen(H)
  #for (i in 1:2) {
  #    delta[i] <- (b[i] - a[i]) / (M[i] - 1)
  #    if (all(H[!diag(nrow(H))] == 0)) L[i] <- min(abs(floor(supp*sqrt(eig$values[i])/delta[i])), M[i]-1)
  #    else  L[i] <- min(abs(floor(supp*sqrt(max(eig$values))*eig$vectors[i,1] / delta[i])), M[i] - 1)
  #}
  if (min(L)==0) warning("Binning grid too coarse for current (small) bandwidth: consider increasing gridsize")

  if(d==2)
  {
      grid1 <- seq(-(L[1]-1), L[1]-1)
      grid2 <- seq(-(L[2]-1), L[2]-1)
      xgrid <- expand.grid(delta[1]*grid1, delta[2]*grid2)
  }
  else if(d==3)
  {
     grid1 <- seq(-(L[1]-1), L[1]-1)
     grid2 <- seq(-(L[2]-1), L[2]-1)
     grid3 <- seq(-(L[3]-1), L[3]-1)
     xgrid <- expand.grid(delta[1]*grid1, delta[2]*grid2, delta[3]*grid3)
  }
  else if (d==4)
  {
     grid1 <- seq(-(L[1]-1), L[1]-1)
     grid2 <- seq(-(L[2]-1), L[2]-1)
     grid3 <- seq(-(L[3]-1), L[3]-1)
     grid4 <- seq(-(L[4]-1), L[4]-1)
     xgrid <- expand.grid(delta[1]*grid1, delta[2]*grid2, delta[3]*grid3, delta[4]*grid4)
  }
  
  ##if (d==2) xgrid <- expand.grid((b[1]-a[1])*(0:L[1])/M[1], (b[2]-a[2])*(0:L[2])/M[2])
  ##if (d==3) xgrid <- expand.grid((b[1]-a[1])*(0:L[1])/M[1], (b[2]-a[2])*(0:L[2])/M[2], (b[3]-a[3])*(0:L[3])/M[3])
  ##if (d==4) xgrid <- expand.grid((b[1]-a[1])*(0:L[1])/M[1], (b[2]-a[2])*(0:L[2])/M[2], (b[3]-a[3])*(0:L[3])/M[3], (b[4]-a[4])*(0:L[4])/M[4])
  
  deriv.index <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, only.index=TRUE, deriv.vec=TRUE) 
  deriv.index.minimal <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, only.index=TRUE, deriv.vec=FALSE)

  keval <- dmvnorm.deriv(x=xgrid, mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, deriv.vec=FALSE)
  keval <- keval$deriv/n
  if (r==0) keval <- as.matrix(keval, ncol=1)
  est <- list()
  if (verbose) pb <- txtProgressBar() 

  ## loop over only unique partial derivative indices
  nderiv <- nrow(deriv.index.minimal)
  if (!(is.null(nderiv)))
    for (s in 1:nderiv)
    {
      if (deriv.vec) deriv.rep.index <- which.mat(deriv.index.minimal[s,], deriv.index)
      else deriv.rep.index <- s
      kevals <- array(keval[,s], dim=N)
      if (r==0) sf <- rep(1,d)
      else sf <- (-1)^deriv.index.minimal[s,]
      est.temp <- symconv.nd(kevals, bin.par$counts, d=d)
      for (s2 in 1:length(deriv.rep.index)) est[[deriv.rep.index[s2]]] <- est.temp ##zapsmall(est.temp)
      if (verbose) setTxtProgressBar(pb, s/nderiv)
    }
  else
  {
    for (s in 1:ncol(keval))
    {  
      kevals <- array(keval[,s], dim=N)
      if (r==0) sf <- rep(1,d)
      else sf <- (-1)^deriv.index[s,]
      est[[s]] <- symconv.nd(kevals, bin.par$counts, d=d)
      ##est[[s]] <- zapsmall(est[[s]])
      if (verbose) setTxtProgressBar(pb, s/ncol(keval))
    }
  }
  if (verbose) close(pb)
 
  return(list(eval.points=bin.par$eval.points, estimate=est, deriv.order=r))
}


#############################################################################
#### Univariate kernel density derivative estimate on a grid
#############################################################################

kdde.grid.1d <- function(x, h, gridsize, supp=3.7, positive=FALSE, adj.positive, xmin, xmax, gridtype, w, deriv.order=0)
{
  r <- deriv.order
  if (r==0)
    fhatr <- kde(x=x, h=h, gridsize=gridsize, supp=supp, positive=positive, adj.positive=adj.positive, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w)
  else
  {  
    if (missing(xmin)) xmin <- min(x) - h*supp
    if (missing(xmax)) xmax <- max(x) + h*supp
    if (missing(gridtype)) gridtype <- "linear"
  
    y <- x
    gridtype1 <- match.arg(gridtype, c("linear", "sqrt"))
    if (gridtype1=="linear")
     gridy <- seq(xmin, xmax, length=gridsize)
    else if (gridtype1=="sqrt")
    {
      gridy.temp <- seq(sign(xmin)*sqrt(abs(xmin)), sign(xmax)*sqrt(abs(xmax)), length=gridsize)
      gridy <- sign(gridy.temp) * gridy.temp^2
    }
    gridtype.vec <- gridtype1
    
    n <- length(y)
    est <- dnorm.deriv.mixt(x=gridy, mus=y, sigmas=rep(h, n), props=w/n, deriv.order=r)
    fhatr <- list(x=y, eval.points=gridy, estimate=est, h=h, H=h^2, gridtype=gridtype.vec, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=deriv.order)
      
    class(fhatr) <- "kde"
  }
  
  return(fhatr)
}



##############################################################################
## Bivariate kernel density derivative estimate on a grid
## Computes all mixed partial derivatives for a given deriv.order
##############################################################################

kdde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, deriv.order=0, deriv.vec=TRUE, verbose=FALSE)
{
  d <- 2
  r <- deriv.order
  if (r==0)
    fhatr <- kde(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose)
  else
  {  
    ## initialise grid 
    n <- nrow(x)
    if (is.null(gridx))
      gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
    
    suppx <- make.supp(x, matrix.sqrt(H), tol=supp)
    
    if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    

    nderiv <- d^r
    fhat.grid <- list()
    for (k in 1:nderiv)
      fhat.grid[[k]] <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))
    if (verbose) pb <- txtProgressBar()
    for (i in 1:n)
    {
      ## compute evaluation points 
      eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
      eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
      eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
      eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
      eval.x.len <- length(eval.x)
      eval.pts <- permute(list(eval.x, eval.y))

      ## Create list of matrices for different partial derivatives
      fhat <- dmvnorm.deriv(x=eval.pts, mu=x[i,], Sigma=H, deriv.order=r)
      
      ## place vector of density estimate values `fhat' onto grid 'fhat.grid'
      for (k in 1:nderiv)
        for (j in 1:length(eval.y))
          fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] <- fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] + w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len),k]
      if (verbose) setTxtProgressBar(pb, i/n)
    }
    if (verbose) close(pb)
    
    for (k in 1:nderiv) fhat.grid[[k]] <- fhat.grid[[k]]/n
    gridx1 <- list(gridx[[1]], gridx[[2]]) 

    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, only.index=TRUE)

    if (!deriv.vec)
    {
      fhat.grid.vech <- list()
      deriv.ind <- unique(ind.mat)

      for (i in 1:nrow(deriv.ind))
      {
        which.deriv <- which.mat(deriv.ind[i,], ind.mat)[1]
        fhat.grid.vech[[i]] <- fhat.grid[[which.deriv]]
      }
      ind.mat <- deriv.ind
      fhat.grid <- fhat.grid.vech
    }

    fhatr <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=deriv.order, deriv.ind=ind.mat)
  }

  return(fhatr)
}

kdde.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, deriv.order=0, deriv.vec=TRUE, verbose=FALSE)
{
  d <- 3
  r <- deriv.order
  if (r==0)
    fhatr <- kde(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose)
  else
  {
     
    ## initialise grid 
    n <- nrow(x)
    if (is.null(gridx))
      gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
    
    suppx <- make.supp(x, matrix.sqrt(H), tol=supp)
    
    if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    

    nderiv <- d^r
    fhat.grid <- list()
    for (k in 1:nderiv)
      fhat.grid[[k]] <- array(0, dim=gridsize)
    if (verbose) pb <- txtProgressBar()
    for (i in 1:n)
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

      ## Create list of matrices for different partial derivatives
      ##fhat <- dmvnorm.deriv(x=eval.pts, mu=x[i,], Sigma=H, deriv.order=r)

      ## place vector of density estimate values `fhat' onto grid 'fhat.grid'
      for (ell in 1:nderiv)
          for (k in 1:length(eval.z))
          {
              fhat <- w[i]*dmvnorm.deriv(cbind(eval.pts, eval.z[k]), x[i,], H, deriv.order=r)
              for (j in 1:length(eval.y))
                  fhat.grid[[ell]][eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
                      fhat.grid[[ell]][eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                          fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
          }
      
      if (verbose) setTxtProgressBar(pb, i/n)
    }
    
    if (verbose) close(pb)
    
    for (k in 1:nderiv) fhat.grid[[k]] <- fhat.grid[[k]]/n
    ##gridx1 <- list(gridx[[1]], gridx[[2]]) 
    gridx1 <- list(gridx[[1]], gridx[[2]], gridx[[3]])
    
    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, only.index=TRUE)

    if (!deriv.vec)
    {
      fhat.grid.vech <- list()
      deriv.ind <- unique(ind.mat)

      for (i in 1:nrow(deriv.ind))
      {
        which.deriv <- which.mat(deriv.ind[i,], ind.mat)[1]
        fhat.grid.vech[[i]] <- fhat.grid[[which.deriv]]
      }
      ind.mat <- deriv.ind
      fhat.grid <- fhat.grid.vech
    }

    
    fhatr <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=deriv.order, deriv.ind=ind.mat)
  }

  return(fhatr)
}


#############################################################################
## Multivariate kernel density estimate using normal kernels,
## evaluated at each sample point
#############################################################################

kdde.points.1d <- function(x, h, eval.points, w, deriv.order=0) 
{
  r <- deriv.order
  n <- length(x)
  fhat <- dnorm.deriv.mixt(x=eval.points, mus=x, sigmas=rep(h,n), props=w/n, deriv.order=r)
  
  return(list(x=x, eval.points=eval.points, estimate=fhat, h=h, H=h^2, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=r))
}


kdde.points <- function(x, H, eval.points, w, deriv.order=0, deriv.vec=TRUE) 
{
  n <- nrow(x)
  ##Hs <- numeric(0)
  ##for (i in 1:n)
  ##  Hs <- rbind(Hs, H)
  Hs <- replicate(n, H, simplify=FALSE) 
  Hs <- do.call(rbind, Hs)
  r <- deriv.order
  fhat <- dmvnorm.deriv.mixt(x=eval.points, mus=x, Sigmas=Hs, props=w/n, deriv.order=r, deriv.vec=deriv.vec, add.index=TRUE)
  
  return(list(x=x, eval.points=eval.points, estimate=fhat$deriv, H=H, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=fhat$deriv.ind))
}

#############################################################################
## Plot method for KDDE 
#############################################################################

plot.kdde <- function(x, ...)
{
  fhat <- x
  opr <- options()$preferRaster; if (!is.null(opr)) if (!opr) options("preferRaster"=TRUE)
  if (is.null(fhat$deriv.order))
  {
    class(fhat) <- "kde"
    plot(fhat, ...)
  }
  else
  {  
    if (is.vector(fhat$x))
    {
      plotkdde.1d(fhat, ...)
      invisible()
    }
    else
    {
      d <- ncol(fhat$x)
      if (d==2) 
      {
        plotret <- plotkdde.2d(fhat, ...)
        invisible(plotret)
      }
      else if (d==3)
      {
        ##stop("plot.kdde not yet implemented for 3-d data")
        plotkdde.3d(fhat, ...)
        invisible()
      }
      else 
        stop ("Plot function only available for 1, 2 or 3-d data")
    }
  }
  if (!is.null(opr)) options("preferRaster"=opr)
}

plotkdde.1d <- function(fhat, ylab="Density derivative function", ...)
{
  class(fhat) <- "kde"
  plot(fhat, ylab=ylab, ...)
}


plotkdde.2d <- function(fhat, which.deriv.ind=1, cont=c(25,50,75), abs.cont, display="slice", zlab="Density derivative function",...)
{
  disp1 <- match.arg(display, c("persp", "slice", "image", "filled.contour", "filled.contour2")) 
  
  if (disp1=="slice" | disp1=="filled.contour" | disp1=="filled.contour2")
  {
    if (missing(abs.cont))
    {
      abs.cont <- as.matrix(contourLevels(fhat, approx=TRUE, cont=cont), ncol=length(cont))
      abs.cont <- c(abs.cont[1,], rev(abs.cont[2,]))
    }
  } 
 
  fhat.temp <- fhat 
  fhat.temp$deriv.ind <- fhat.temp$deriv.ind[which.deriv.ind,]
  fhat.temp$estimate <- fhat.temp$estimate[[which.deriv.ind]]
  fhat <- fhat.temp
  class(fhat) <- "kde"

  plot(fhat, display=display, abs.cont=abs.cont, zlab=zlab,  ...) 
}



plotkdde.3d <- function(fhat, which.deriv.ind=1, cont=c(25,50,75), abs.cont, colors, col.fun=topo.colors, ...)
{
  if (missing(abs.cont))
  {
      abs.cont <- as.matrix(contourLevels(fhat, approx=TRUE, cont=cont), ncol=length(cont))
      abs.cont <- c(abs.cont[1,], rev(abs.cont[2,]))
  }
  
  fhat.temp <- fhat 
  fhat.temp$deriv.ind <- fhat.temp$deriv.ind[which.deriv.ind,]
  fhat.temp$estimate <- fhat.temp$estimate[[which.deriv.ind]]
  fhat <- fhat.temp
  class(fhat) <- "kde"

  if (missing(colors))
  {
      colors <- rev(col.fun(length(abs.cont)))
      nc <- length(colors)
      colors[(nc/2+1):nc] <- rev(colors[(nc/2+1):nc])
  }
  plot(fhat, abs.cont=abs.cont, colors=colors, ...) 
}


  
#############################################################################
## ContourLevels method for KDDE 
#############################################################################

contourLevels.kdde <- function(x, prob, cont, nlevels=5, approx=TRUE, which.deriv.ind=1, ...)
{ 
  fhat <- x
  if (is.vector(fhat$x))
  {
    d <- 1; n <- length(fhat$x)
    if (!is.null(fhat$deriv.order))  
    {
      fhat.temp <- fhat 
      fhat.temp$deriv.ind <-fhat.temp$deriv.ind[which.deriv.ind]
      fhat <- fhat.temp
    }    
  }
  else
  {
    d <- ncol(fhat$x); n <-nrow(fhat$x)
    if (!is.matrix(fhat$x)) fhat$x <- as.matrix(fhat$x)

    if (!is.null(fhat$deriv.order))  
    {
      fhat.temp <- fhat 
      fhat.temp$estimate <- fhat.temp$estimate[[which.deriv.ind]]
      fhat.temp$deriv.ind <-fhat.temp$deriv.ind[which.deriv.ind,]
      fhat <- fhat.temp
    }
  } 

  if (is.null(x$w)) w <- rep(1, n)
  else w <- x$w

  if (is.null(fhat$gridded))
  {
    if (d==1) fhat$gridded <- fhat$binned
    else fhat$gridded <- is.list(fhat$eval.points)
  }


  if (missing(prob) & missing(cont))
    hts <- pretty(fhat$estimate, n=nlevels) 
  else
  {
    if (approx & fhat$gridded)
      dobs <- predict.kde(fhat, x=fhat$x)
    else
        dobs <- kdde(x=fhat$x, H=fhat$H, eval.points=fhat$x, w=w, deriv.order=fhat$deriv.order)$estimate[,which.deriv.ind] 
    
    if (is.null(fhat$deriv.order))
    {
       if (!missing(prob) & missing(cont)) hts <- quantile(dobs[dobs>=0], prob=prob)
       if (missing(prob) & !missing(cont)) hts <- quantile(dobs[dobs>=0], prob=(100-cont)/100)
    }
    else
    {   
      if (!missing(prob) & missing(cont)) hts <- rbind(quantile(dobs[dobs<0], prob=prob), quantile(dobs[dobs>=0], prob=prob))
      if (missing(prob) & !missing(cont)) hts <- rbind(quantile(dobs[dobs<0], prob=(100-cont)/100), quantile(dobs[dobs>=0], prob=(100-cont)/100))
    }
  }
  
  return(hts)
}

#############################################################################
## predict method for KDDE 
#############################################################################

predict.kdde <- function(object, ..., x)
{
  fhat <- object
  if (is.vector(fhat$H)) d <- 1 else d <- ncol(fhat$H)
  if (d==1) n <- length(x)
  else
  {
    x <- as.matrix(x) 
    n <- nrow(x)
  }
  
  if (!is.null(fhat$deriv.ind))
  {
    if (is.vector(fhat$deriv.ind)) pk.mat <- predict.kde(fhat, x=x, ...)
    else
    {
      nd <- nrow(fhat$deriv.ind)
      pk.mat <- matrix(0, ncol=nd, nrow=n)
      for (i in 1:nd)
        {
          fhat.temp <- fhat
          fhat.temp$estimate <- fhat$estimate[[i]]
          pk.mat[,i] <- predict.kde(fhat.temp, x=x, ...)
        }
    }
  }
  else
  {
    pk.mat <- predict.kde(fhat, x=x, ...)
  }
  return(drop(pk.mat))
}


