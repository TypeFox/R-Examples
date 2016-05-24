
##==============================================================================
## Creation of a one-dimensional finite difference grid
##==============================================================================

setup.grid.1D <- function(x.up=0,	x.down=NULL, L=NULL,
  N=NULL, dx.1 =NULL, p.dx.1 = rep(1,length(L)),
  max.dx.1 = L, dx.N =NULL, p.dx.N = rep(1,length(L)),
  max.dx.N = L) {

## Check on the input
  if (is.null(L)){
    if (is.null(N) && is.null(dx.1))
      stop ("either the length (L) and the end point (x.down) or the number of boxes,(N) and size of boxes (dx.1 should be specified")
    if (is.null(N) && ! is.null(dx.1)) 
      if (is.null(x.down))
       stop("either the length (L) or the end points (x.up, x.down) and N should be specified")  
    if (! is.null(x.down))
      L = x.down-x.up
    else
      L = N * dx.1
  }
  
  for (i in 1:length(L)) {
  if (is.null(x.down[i]) && is.null(L[i]))
    stop (paste("Error in setup.grid.1D! Grid zone",i,": either the length (L) or the end point (x.down) should be specified"))
  }

## If the interval lengths are given, create the end points of each zone 

  if (is.null(x.down[1])) {
    if (any(L<0)) stop (paste("Error in setup.grid.1D! L[",which(L<0),"] < 0",sep=""))
    x.down <- x.up + cumsum(L)
  }
  
## If the interval lengths are given, create the end points of each zone 

  if (is.null(L[1])) {
    L <- diff(c(x.up,x.down))
    if (any(L<0)) stop (paste("Error in setup.grid.1D! L[",which(L<0),"] < 0",sep=""))
  }

## Check wether all info is available in each zone 

  for (i in 1:length(L)) {
    if (is.null(N[i]) && is.null(dx.1[i]) && is.null(dx.N[i]))
    stop (paste("Error in setup.grid.1D! Both N[",i,"], dx.1[",i,"] and dx.N[",i,"] are NULL"))
  }

## Calculation of the grid cell sizes

  dx <- vector()

## Power law function that controls the increase in grid size. The 
## root of this function determines the power factor that corresponds to a 
## desired number of cells N 

  p1 <- 1.00001
  p2 <- 3
  
  f.root <- function(p,dx,N,L) dx*(p^(N)-1)/(p-1) - L

## Loop over all grid zones
  for (i in 1:length(L)) {

## Option 1: only N[i] is specified, all grid cells are equal
    if (!(is.null(N[i])) && (is.null(dx.1[i])) && (is.null(dx.N[i])))  { 
      if (N[i] < 1) stop (paste("Error in setup.grid.1D! N[",i,"] < 1",sep=""))
      A.dx  <- rep(L[i]/N[i],N[i])
    }
    
## Option 2: dx.1[i] specified, N[i] = NULL and dx.N[i] = NULL 

    if ((is.null(N[i])) && !(is.null(dx.1[i])) && (is.null(dx.N[i])))  { 
      if (dx.1[i]<=0)
        stop (paste("Error in setup.grid.1D! dx.1[",i,"] <= 0",sep=""))
      if (dx.1[i] > L[i])
        stop (paste("Error in setup.grid.1D! dx.1[",i,"] > L[",i,"]",sep=""))
  
     # use gradual increase of grid cell size at upstream interface 
      A.dx <- vector()
      pos <- vector()
      A.dx[1] <- dx.1[i]
      pos[1] <- dx.1[i]
      j <- 1
      while (pos[j] < L[i]) {
        j <- j + 1
        A.dx[j] <- min(max.dx.1[i],A.dx[j-1]*p.dx.1[i])
        pos[j] <- pos[j-1] + A.dx[j]
      }
     # rescaling to fit the interval length 
      A.dx <- A.dx*(L[i]/pos[length(pos)])
    }

## Option 3: dx.N[i] specified, N[i] = NULL and dx.1[i] = NULL 

    if ((is.null(N[i])) && (is.null(dx.1[i])) && !(is.null(dx.N[i])))  { 
      if (dx.N[i]<=0)
        stop (paste("Error in setup.grid.1D! dx.N[",i,"] <= 0",sep=""))
      if (dx.N[i] > L[i])
        stop (paste("Error in setup.grid.1D! dx.N[",i,"] > L[",i,"]",sep=""))
  
     # use gradual increase of grid cell size at downstream interface 
      A.dx <- vector()
      pos <- vector()
      A.dx[1] <- dx.N[i]
      pos[1] <- dx.N[i]
      j <- 1
      while (pos[j] < L[i]) {
        j <- j + 1
        A.dx[j] <- min(max.dx.N[i],A.dx[j-1]*p.dx.N[i])
        pos[j] <- pos[j-1] + A.dx[j]
      }
     # rescaling to fit the interval length 
      A.dx <- A.dx*(L[i]/pos[length(pos)])
     # reversing the A.dx vector
      A.dx <- A.dx[length(A.dx):1]
    }

## Option 4: N[i] and dx.1[i] are specified, dx.N[i] = NULL 

    if (!(is.null(N[i])) && !(is.null(dx.1[i])) && (is.null(dx.N[i])))  { 

      if (N[i] < 1) stop (paste("Error in setup.grid.1D! N[",i,"] < 1",sep=""))
      if (dx.1[i]<=0)
        stop (paste("Error in setup.grid.1D! dx.1[",i,"] <= 0",sep=""))
      if (dx.1[i] > L[i])
        stop (paste("Error in setup.grid.1D! dx.1[",i,"] > L[",i,"]",sep=""))
  
      # estimate power in power law 
      if (f.root(p=p1,dx.1[i],N[i],L[i])*f.root(p=p2,dx.1[i],N[i],L[i]) >= 0)
      {
      A.dx <- rep(L[i]/N[i],times=N[i])
      } else {
      p.estim <- uniroot(f=f.root,dx=dx.1[i],N=N[i],L=L[i],lower=p1,upper=p2)$root
      # use gradual increase of grid cell size at upper interface 
      A.dx <- vector()
      pos <- vector()
      A.dx[1] <- dx.1[i]
      pos[1] <- dx.1[i]
      j <- 1
      while (j < N[i]) {
        j <- j + 1
        A.dx[j] <- min(max.dx.1[i],A.dx[j-1]*p.estim)
        pos[j] <- pos[j-1] + A.dx[j]
      }
     # rescaling to fit the interval length 
      A.dx <- A.dx*(L[i]/pos[length(pos)])
      }
   }

## Option 5: N[i] and dx.N[i] are specified, dx.1[i] = NULL 

    if (!(is.null(N[i])) && (is.null(dx.1[i])) && !(is.null(dx.N[i])))  { 

      if (N[i] < 1) stop (paste("Error in setup.grid.1D! N[",i,"] < 1",sep=""))
      if (dx.N[i]<=0)
        stop (paste("Error in setup.grid.1D! dx.N[",i,"] <= 0",sep=""))
      if (dx.N[i] > L[i])
        stop (paste("Error in setup.grid.1D! dx.N[",i,"] > L[",i,"]",sep=""))
  
      if (f.root(p=p1,dx.N[i],N[i],L[i])*f.root(p=p2,dx.N[i],N[i],L[i]) >= 0)
      {
      A.dx <- rep(L[i]/N[i],times=N[i])
      } else {
     # estimate power in power law 
      p.estim <- uniroot(f=f.root,dx=dx.N[i],N=N[i],L=L[i],lower=p1,upper=p2)$root
      # use gradual increase of grid cell size at upper interface 
      A.dx <- vector()
      pos <- vector()
      A.dx[1] <- dx.N[i]
      pos[1] <- dx.N[i]
      j <- 1
      while (j < N[i]) {
        j <- j + 1
        A.dx[j] <- min(max.dx.N[i],A.dx[j-1]*p.estim)
        pos[j] <- pos[j-1] + A.dx[j]
      }
     # rescaling to fit the interval length 
      A.dx <- A.dx*(L[i]/pos[length(pos)])
      }
     # reversing the A.dx vector
      A.dx <- A.dx[length(A.dx):1]
    }

## Option 6: dx.1[i] and dx.N[i] are specified, N[i] can be NULL 

    if (!(is.null(dx.1[i])) && !(is.null(dx.N[i])))  { 

      if (dx.1[i]<=0)
        stop (paste("Error in setup.grid.1D! dx.1[",i,"] <= 0",sep=""))
      if (dx.1[i] > L[i]/2)
        stop (paste("Error in setup.grid.1D! dx.1[",i,"] > L[",i,"]",sep=""))
      if (dx.N[i]<=0)
        stop (paste("Error in setup.grid.1D! dx.N[",i,"] <= 0",sep=""))
      if (dx.N[i] > L[i]/2)
        stop (paste("Error in setup.grid.1D! dx.N[",i,"] > L[",i,"]",sep=""))

      L.A <- 0.5*L[i]
      L.B <- 0.5*L[i]

      # estimate power in power law 
      if (!is.null(N[i])) {
      
      N.A <- as.integer(N[i]/2) 
      N.B <- N[i] - N.A
      L.A <- (N.A/N[i])*L[i]
      L.B <- L[i] - L.A
      
      if (f.root(p=p1,dx.1[i],N.A,L.A)*f.root(p=p2,dx.1[i],N.A,L.A) >= 0)
      {
        p.dx.1[i] <- 1
        } else {
        p.dx.1[i] <- uniroot(f=f.root,dx=dx.1[i],N=N.A,L=L.A,lower=p1,upper=p2)$root
      }
      if (f.root(p=p1,dx.N[i],N.B,L.B)*f.root(p=p2,dx.N[i],N.B,L.B) >= 0)
      {
        p.dx.N[i] <- 1
        } else {
        p.dx.N[i] <- uniroot(f=f.root,dx=dx.N[i],N=N.B,L=L.B,lower=p1,upper=p2)$root
      }
      }
      # use gradual increase of grid cell size at upsteam interface
      A.dx  <- vector()
      pos <- vector()
      A.dx[1] <- dx.1[i]
      pos[1] <- dx.1[i]
      j <- 1
      while ((is.null(N[i])&&(pos[j] < L.A))||(!is.null(N[i])&&(j < N.A))) {
         j <- j + 1
         A.dx[j] <- min(max.dx.1[i],A.dx[j-1]*p.dx.1[i])
         pos[j] <- pos[j-1] + A.dx[j]
      }
      A.dx <- A.dx*(L.A/pos[length(pos)])

      # use gradual increase of grid cell size at downsteam interface
      B.dx  <- vector()
      pos <- vector()
      B.dx[1] <- dx.N[i]
      pos[1] <- dx.N[i]
      j <- 1
      while ((is.null(N[i])&&(pos[j] < L.B))||(!is.null(N[i])&&(j < N.B))) {
        j <- j + 1
        B.dx[j] <- min(max.dx.N[i],B.dx[j-1]*p.dx.N[i])
        pos[j] <- pos[j-1] + B.dx[j]
      }
      B.dx <- B.dx*(L.B/pos[length(pos)])
      # assemble the distance vector 
      A.dx <- c(A.dx,B.dx[length(B.dx):1])
    }

  ## calculation of present zone has finished, add distances to existing array 
  
    dx <- c(dx,A.dx)
  }

  ## positions and distances

  x.int  <- x.up + diffinv(dx)
  x.mid  <- colMeans(rbind(x.int[-1],x.int[-(length(dx)+1)]))
  dx.aux <- c(dx[1]/2,diff(x.mid),dx[length(dx)]/2)

  ## Packaging of results

  Res <- list(x.up = x.up,
	  			    x.down = x.down[length(x.down)],
		  		    x.mid = x.mid,   # position of centre of the grid cells, vector of length N
              x.int  = x.int,  # position of the grid cell interfaces , vector of length N+1
              dx = dx  ,       # thickness of the grid cells , vector length N
              dx.aux = dx.aux, # auxiliary vector with distances between centre of adjacent cells, first and last: half of cell size, vector of length N+1
              N = length(dx))  # total number of grid cells
  class(Res) <- "grid.1D"
  return(Res)
}

##==============================================================================
## S3 method: Plotting of a one-dimensional finite difference grid
##==============================================================================

plot.grid.1D <- function(x,...) {
  mf <- par(mfrow=c(2,1))
  on.exit(par(mf))
    
  plot(x$x.int,main="position of cells",ylab="x",xlab="index",...)
  points(0.5+(1:x$N),x$x.mid,pch=16)
  legend("topleft",c("x.int","x.mid"),pch=c(1,16))
  plot(x$dx.aux,main="box thickness",ylab="dx",xlab="index",...)
  points(0.5+(1:x$N),x$dx,pch=16)
  legend("topleft",c("dx.aux","dx.mid"),pch=c(1,16))
}
