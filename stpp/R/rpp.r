.make.grid <- function(nx,ny,poly)
  {
    if (missing(poly)) poly <- matrix(c(0,0,1,0,1,1,0,1),4,2,T)
    
    if ((nx < 2) || (ny < 2)) stop("the grid must be at least of size 2x2")
    
    xrang <- range(poly[, 1], na.rm = TRUE)
    yrang <- range(poly[, 2], na.rm = TRUE)
    xmin <- xrang[1]
    xmax <- xrang[2]
    ymin <- yrang[1]
    ymax <- yrang[2]
    
    xinc <- (xmax-xmin)/nx
    yinc <- (ymax-ymin)/ny
    
    xc <- xmin-xinc/2
    yc <- ymin-yinc/2
    xgrid <- rep(0,nx)
    ygrid <- rep(0,ny)
    xgrid[1] <- xc + xinc
    ygrid[1] <- yc + yinc

    for (i in 2:nx)
      {
        xgrid[i] <- xgrid[i-1]+xinc
      }
    for (i in 2:ny)
      {
        ygrid[i] <- ygrid[i-1]+yinc
      }

    yy <- matrix(xgrid,nx,ny)
    xx <- t(yy)
    yy <- matrix(ygrid,nx,ny)

    X <- as.vector(xx)
    Y <- as.vector(yy)

    poly <- rbind(poly,poly[1,])
    pts <- inpip(pts=cbind(X,Y),poly)

    X[pts] <- TRUE
    X[X!=TRUE] <- FALSE
    mask <- matrix(X,ncol=ny,nrow=nx,byrow=TRUE)

    invisible(return(list(x=xgrid,y=ygrid,X=xx,Y=yy,pts=pts,xinc=xinc,yinc=yinc,mask=matrix(as.logical(mask),nx,ny))))
  }

.rhpp <- function(lambda, s.region, t.region, npoints=NULL, replace=TRUE, discrete.time=FALSE)
  {
    if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
    if (missing(t.region)) t.region <- c(0,1)
    
    s.area <- areapl(s.region)
    t.region <- sort(t.region)
    t.area <- t.region[2]-t.region[1]

    if (missing(lambda) & !(is.null(npoints)))
      {
        if (t.area==0) lambda <- npoints/(s.area)
        else lambda <- npoints/(s.area * t.area)
      }
    
    pattern <- list()
    index.t <- list()
   if (is.numeric(lambda))
    {
      if (is.null(npoints)==TRUE)
        {
          if (t.area==0)
            { npoints <- round(rpois(n=1,lambda=lambda * s.area),0) }
          else
            { npoints <- round(rpois(n=1,lambda=lambda * s.area * t.area),0) }
        }
      xy <- matrix(csr(poly=s.region,npoints=npoints),ncol=2)
      x <- xy[,1]
      y <- xy[,2]
      npoints <- length(x)
      
      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          if ((length(vect)<npoints) & (replace==FALSE))
            stop("when replace=FALSE and discrete.time=TRUE, the length of seq(t.region[1],t.region[2],by=1) must be greater than the number of points")
          names(vect) <- 1:length(vect) 
          M <- sample(vect,npoints,replace=replace)
          times <- M
          names(times) <- NULL
          samp <- as.numeric(names(M))
        }
      else
        {
          times <- runif(npoints,min=t.region[1],max=t.region[2])
          samp <- sample(1:npoints,npoints,replace=replace)
          times <- times[samp]
        }
      times <- sort(times)
      index.times <- sort(samp)
      pattern.interm <- list(x=x,y=y,t=times,index.t=index.times)      
      pattern <- cbind(x=pattern.interm$x,y=pattern.interm$y,t=pattern.interm$t)
      index.t <- pattern.interm$index.t
    }
   else
     stop("lambda must be numeric")

    invisible(return(list(pts=pattern,index.t=index.t)))
  }



.ripp <- function(lambda, s.region, t.region, npoints=NULL, replace=TRUE, discrete.time=FALSE, nx=100, ny=100, nt=100, lmax=NULL, Lambda=NULL, ...)
{
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  s.area <- areapl(s.region)
  t.region <- sort(t.region)
  t.area <- t.region[2]-t.region[1]

  if (missing(lambda) & !(is.null(npoints)))
    {
      if (t.area==0) lambda <- npoints/(s.area)
      else lambda <- npoints/(s.area * t.area)
    }

  lambdamax <- lmax
  pattern <- list()
  index.t <- list()

  if (is.function(lambda))
    {
      s.grid <- .make.grid(nx,ny,s.region)
      s.grid$mask <- matrix(as.logical(s.grid$mask),nx,ny)
      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          if (nt>length(vect))
            {
              nt <- length(vect)
              warning("nt used is less than the one given in argument")
              t.grid <- list(times=vect,tinc=1)
            }
          else
            {
              vect <- round(seq(floor(t.region[1]),ceiling(t.region[2]),length=nt))
              t.grid <- list(times=vect,tinc=round(t.area/(nt-1)))
            }
        }
      else
        t.grid <- list(times=seq(t.region[1],t.region[2],length=nt),tinc=(t.area/(nt-1)))

      if (is.null(Lambda))
        {
	    Lambda <- array(NaN,dim=c(nx,ny,nt)) 
          for(it in 1:nt)
            {
              L <- lambda(as.vector(s.grid$X),as.vector(s.grid$Y),t.grid$times[it],...)
              M <- matrix(L,ncol=ny,nrow=nx,byrow=TRUE)
              M[!(s.grid$mask)] <- NaN
              Lambda[,,it] <- M
            }
        }
      
      if (is.null(npoints)==TRUE)
        {
          if (t.area==0)
            { en <- sum(Lambda,na.rm=TRUE)*s.grid$xinc*s.grid$yinc }
          else
            {
              en <- sum(Lambda,na.rm=TRUE)*s.grid$xinc*s.grid$yinc*t.grid$tinc 
              npoints <- round(rpois(n=1,lambda=en),0)
            }
        }

      if (is.null(lambdamax))
        lambdamax <- max(Lambda,na.rm=TRUE)
      npts <- round(lambdamax/(s.area*t.area),0)
      if (npts==0) stop("there is no data to thin")
      
      xy <- matrix(csr(poly=s.region,npoints=npts),ncol=2)
      x <- xy[,1]
      y <- xy[,2]
      
      if ((replace==FALSE) & (nt < max(npts,npoints))) stop("when replace=FALSE, nt must be greater than the number of points used for thinning")
      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          times.init <- sample(vect,nt,replace=replace)
        }
      else
        times.init <- runif(nt,min=t.region[1],max=t.region[2])
      
      samp <- sample(1:nt,npts,replace=replace)
      times <- times.init[samp]
      prob <-  lambda(x,y,times,...)/lambdamax
      u <- runif(npts)
      retain <- u <= prob
      if (sum(retain==FALSE)==length(retain))
        {
          lambdas <- matrix(0,nrow=nx,ncol=ny)
          for(ix in 1:nx){for(iy in 1:ny){
            lambdas[ix,iy] <- median(Lambda[ix,iy,],na.rm=TRUE)}}
          lambdamax <- max(lambdas,na.rm=TRUE)
          prob <-  lambda(x,y,times,...)/lambdamax
          retain <- u <= prob
          if (sum(retain==F)==length(retain)) stop ("no point was retained at the first iteration, please check your parameters")
        }
      x <- x[retain]
      y <- y[retain]
      samp <- samp[retain]
      samp.remain <- (1:nt)[-samp]
      times <- times[retain]
      
      neffec <- length(x)
      while(neffec < npoints)
        {
          xy <- as.matrix(csr(poly=s.region,npoints=npoints-neffec))
          if(dim(xy)[2]==1){wx <- xy[1]; wy <- xy[2]}
          else{wx <- xy[,1]; wy <- xy[,2]}
          if(replace==FALSE)
            { wsamp <- sample(samp.remain,npoints-neffec,replace=replace) }
          else{ wsamp <- sample(1:nt,npoints-neffec,replace=replace) }
          wtimes <- times.init[wsamp]
#              lambdamax <- maxlambda[wsamp]
          prob <-  lambda(wx,wy,wtimes,...)/lambdamax
          u <- runif(npoints-neffec)
          retain <- u <= prob
          x <- c(x,wx[retain])
          y <- c(y,wy[retain])
          times <- c(times,wtimes[retain])
          samp <- c(samp,wsamp[retain])
          samp.remain <- (1:nt)[-samp]
          neffec <- length(x)
        }
      times <- sort(times)
      index.times <- sort(samp)
      pattern.interm <- list(x=x,y=y,t=times,index.t=index.times)
          pattern <- cbind(x=pattern.interm$x,y=pattern.interm$y,t=pattern.interm$t)
      index.t <- pattern.interm$index.t
    }

  if (is.character(lambda))
    {
      if (is.null(Lambda))
        stop("Lambda must be specified")
      nx <- dim(Lambda)[1]
      ny <- dim(Lambda)[2]
      nt <- dim(Lambda)[3]
      
      s.grid <- .make.grid(nx,ny,s.region)
      s.grid$mask <- matrix(as.logical(s.grid$mask),nx,ny)

      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          if (nt>length(vect))
            {
              nt <- length(vect)
              warning("nt used is less than the one given in argument")
              t.grid <- list(times=vect,tinc=1)
            }
          else
            {
              vect <- round(seq(floor(t.region[1]),ceiling(t.region[2]),length=nt))
              t.grid <- list(times=vect,tinc=round(t.area/(nt-1)))
            }
        }
      else
        t.grid <- list(times=seq(t.region[1],t.region[2],length=nt),tinc=(t.area/(nt-1)))
          
      if (is.null(npoints))
        {
	    en <- sum(Lambda,na.rm=TRUE)*s.grid$xinc*s.grid$yinc*t.grid$tinc 
          npoints <- round(rpois(n=1,lambda=en),0)
        }
      if (is.null(lambdamax))
        lambdamax <- max(Lambda,na.rm=TRUE)
#      npts <- round(lambdamax/(s.area*t.area),0)
      npts <- npoints
      if (npts==0) stop("there is no data to thin")

      if ((replace==FALSE) & (nt < max(npts,npoints))) stop("when replace=FALSE, nt must be greater than the number of points used for thinning")
      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          times.init <- sample(vect,nt,replace=replace)
        }
      else
        times.init <- runif(nt,min=t.region[1],max=t.region[2])

      samp <- sample(1:nt,npts,replace=replace)
      times <- times.init[samp]

      retain.eq.F <- FALSE
      while(retain.eq.F==FALSE)
        {
          xy <- matrix(csr(poly=s.region,npoints=npts),ncol=2)
          x <- xy[,1]
          y <- xy[,2]
    
          prob <- NULL
          for(nx in 1:length(x))
            {
              nix <- findInterval(vec=s.grid$x,x=x[nx])
              niy <- findInterval(vec=s.grid$y,x=y[nx])
              nit <- findInterval(vec=t.grid$times,x=times[nx])
		  if (nix==0 | niy==0 | nit==0) 
			prob=c(prob,NA)
		  else	
                  prob <- c(prob,Lambda[nix,niy,nit]/lambdamax)
            }
          
          M <- which(is.na(prob))
          if (length(M)!=0)
            {
              x <- x[-M]
              y <- y[-M]
              times <- times[-M]
              prob <- prob[-M]
              npts <- length(x)
            }
          
          u <- runif(npts)
          retain <- u <= prob
          if (sum(retain==F)==length(retain)) retain.eq.F <- FALSE
          else retain.eq.F <- TRUE
        }
#      if (sum(retain==FALSE)==length(retain)) stop ("no point was retained at the first iteration, please check your parameters")
      

      if (sum(retain==FALSE)==length(retain))
        {
          lambdas <- matrix(0,nrow=nx,ncol=ny)
          for(ix in 1:nx){for(iy in 1:ny){
            lambdas[ix,iy] <- median(Lambda[ix,iy,],na.rm=TRUE)}}
          lambdamax <- max(lambdas,na.rm=TRUE)
          prob <-  lambda(x,y,times,...)/lambdamax
          retain <- u <= prob
          if (sum(retain==FALSE)==length(retain)) stop ("no point was retained at the first iteration, please check your parameters")
        }

      x <- x[retain]
      y <- y[retain]
      samp <- samp[retain]
      samp.remain <- (1:nt)[-samp]
      times <- times[retain]

      neffec <- length(x)
      if (neffec > npoints)
        {
          retain <- 1:npoints
          x <- x[retain]
          y <- y[retain]
          samp <- samp[retain]
          samp.remain <- (1:nt)[-samp]
          times <- times[retain]
        }

      while(neffec < npoints)
        {
          xy <- as.matrix(csr(poly=s.region,npoints=npoints-neffec))
          if(dim(xy)[2]==1){wx <- xy[1]; wy <- xy[2]}
          else{wx <- xy[,1]; wy <- xy[,2]}
          if(replace==FALSE)
            { wsamp <- sample(samp.remain,npoints-neffec,replace=replace) }
          else{ wsamp <- sample(1:nt,npoints-neffec,replace=replace) }
          wtimes <- times.init[wsamp]
#              lambdamax <- maxlambda[wsamp]

          prob <- NULL
          for(nx in 1:length(wx))
            {
              nix <- findInterval(vec=s.grid$x,x=wx[nx])
              niy <- findInterval(vec=s.grid$y,x=wy[nx])
              nit <- findInterval(vec=t.grid$times,x=wtimes[nx])
		  if (nix==0 | niy==0 | nit==0) 
			prob=c(prob,NA)
		  else	
                  prob <- c(prob,Lambda[nix,niy,nit]/lambdamax)
            }
          M <- which(is.na(prob))
          if (length(M)!=0)
            {
              wx <- wx[-M]
              wy <- wy[-M]
              wtimes <- wtimes[-M]
              prob <- prob[-M]
            }
          if (neffec > 0)
            {
              u <- runif(length(prob))
              retain <- u <= prob
              x <- c(x,wx[retain])
              y <- c(y,wy[retain])
              times <- c(times,wtimes[retain])
              samp <- c(samp,wsamp[retain])
              samp.remain <- (1:nt)[-samp]
              neffec <- length(x)
            }
        }
      times <- sort(times)
      index.times <- sort(samp)
      pattern.interm <- list(x=x,y=y,t=times,index.t=index.times)
      pattern <- cbind(x=pattern.interm$x,y=pattern.interm$y,t=pattern.interm$t)
      index.t <- pattern.interm$index.t
    }
      
  invisible(return(list(pts=pattern,index.t=index.t)))
}
  

rpp <- function(lambda, s.region, t.region, npoints=NULL, nsim=1, replace=TRUE, discrete.time=FALSE, nx=100, ny=100, nt=100, lmax=NULL, ...)
{

  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  s.area <- areapl(s.region)
  t.region <- sort(t.region)
  t.area <- t.region[2]-t.region[1]

  if (missing(lambda) & !(is.null(npoints)))
    {
      if (t.area==0) lambda <- npoints/(s.area)
      else lambda <- npoints/(s.area * t.area)
    }

  lambdamax <- lmax
  pattern <- list()
  index.t <- list()
  ni <- 1

  #
  # Homogeneous Poisson Process
  #

  if (is.numeric(lambda) & length(lambda)==1)
    {
      while(ni<=nsim)
        {
          hpp <- .rhpp(lambda=lambda, s.region=s.region, t.region=t.region, npoints=npoints, replace=replace, discrete.time=discrete.time)
          if (nsim==1)
            {
              pattern <- as.3dpoints(hpp$pts)
              index.t <- hpp$index.t
            }
          else
            {
              pattern[[ni]] <- as.3dpoints(hpp$pts)
              index.t[[ni]] <- hpp$index.t
            }
          ni <- ni+1
        }
      Lambda <- NULL
    }
    
  #
  # Inhomogeneous Poisson Process
  #

  else if (is.function(lambda))
    {
      s.grid <- .make.grid(nx,ny,s.region)
      s.grid$mask <- matrix(as.logical(s.grid$mask),nx,ny)
      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          if (nt>length(vect))
            {
              nt <- length(vect)
              warning("nt used is less than the one given in argument")
              t.grid <- list(times=vect,tinc=1)
            }
          else
            {
              vect <- round(seq(floor(t.region[1]),ceiling(t.region[2]),length=nt))
              t.grid <- list(times=vect,tinc=round(t.area/(nt-1)))
            }
        }
      else
        t.grid <- list(times=seq(t.region[1],t.region[2],length=nt),tinc=(t.area/(nt-1)))

      Lambda <- array(NaN,dim=c(nx,ny,nt)) 
      for(it in 1:nt)
        {
          L <- lambda(as.vector(s.grid$X),as.vector(s.grid$Y),t.grid$times[it],...)
          M <- matrix(L,ncol=ny,nrow=nx,byrow=TRUE)
          M[!(s.grid$mask)] <- NaN
          Lambda[,,it] <- M
        }
          
      while(ni<=nsim)
        {
          ipp <- .ripp(lambda=lambda, s.region=s.region, t.region=t.region, npoints=npoints, replace=replace, discrete.time=discrete.time, nx=nx, ny=ny, nt=nt, lmax=lmax, Lambda=Lambda, ...)
          
          if (nsim==1)
            {
              pattern <- as.3dpoints(ipp$pts)
              index.t <- ipp$index.t
            }
          else
            {
              pattern[[ni]] <- as.3dpoints(ipp$pts)
              index.t[[ni]] <- ipp$index.t
            }
          ni <- ni+1
        }
    }
     
 else if (is.array(lambda))
    {
      if (length(dim(lambda))!=3) stop ("lambda must be a 3D-array")
      Lambda = lambda
      lambda = "a"

      while(ni<=nsim)
        {
          ipp <- .ripp(lambda=lambda, s.region=s.region, t.region=t.region, npoints=npoints, replace=replace, discrete.time=discrete.time, nx=nx, ny=ny, nt=nt, lmax=lmax, Lambda=Lambda, ...)
          
          if (nsim==1)
            {
              pattern <- as.3dpoints(ipp$pts)
              index.t <- ipp$index.t
            }
          else
            {
              pattern[[ni]] <- as.3dpoints(ipp$pts)
              index.t[[ni]] <- ipp$index.t
            }
          ni <- ni+1
        }
    }
    else stop("lambda must be either a single positive value or a function or a 3D-array")
  
  invisible(return(list(xyt=pattern,index.t=index.t,s.region=s.region,t.region=t.region,lambda=lambda,Lambda=Lambda)))
}



