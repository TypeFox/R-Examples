.set.cov <- function(separable,model,param,sigma2)
  {
    mods <- 0
    modt <- 0
    mod <- 0    

    models <- c("exponential","cauchy","stable","wave","gneiting","cesare","matern","none")

    for(i in 1:length(model))
      {
        M <- which(models==model[i])
        if (length(M)==0) stop("the model is not implemented")
      }

    for (i in 1:length(unique(model)))
      {
        if (((isTRUE(separable)) & ((model[i]==models[5]) | (model[i]==models[6]))) | ((!(isTRUE(separable))) & ((model[i]==models[1]) | (model[i]==models[2]) | (model[i]==models[3]) | (model[i]==models[4]) | (model[i]==models[7])))) stop("'stcov' does not match with 'model'")
      }

    if (isTRUE(separable))
      {
        if ((length(model)!=1) & (length(model)!=2))
          stop("for separable covariance functions, 'model' must be of length 1 or 2")
        if (length(model)==1)
          {
            if (model=="none")
              {
                mods <- 0
                modt <- 0
              }
            if (model=="exponential")
              {
                mods <- 1
                modt <- 1
              }
            if (model=="stable")
              {
                mods <- 2
                if ((param[1] >1) | (param[1]<0)) stop("Stable model parameter must lie in (0,1]")
                modt <- 2
                if ((param[2] >1) | (param[2]<0)) stop("Stable model parameter must lie in (0,1]")
              }
            if (model=="cauchy")
              {
                mods <- 3
                if (param[1]<=0) stop("Cauchy model parameter must be strictly positive")
                modt <- 3
                if (param[2]<=0) stop("Cauchy model parameter must be strictly positive")
              }
            if (model=="wave")
              {
                mods <- 4
                modt <- 4
                }
            if (model=="matern") 
			{
                  mods <- 7
			if (param[2]<=0 | param[1]<=0) stop("Matern model parameters must be strictly positive")
                  modt <- 7
			if (param[3]<=0 | param[4]<=0) stop("Matern model parameters must be strictly positive")
                  }
          }
            if (length(model)==2)
              {
                if (model[1]=="none")
                    mods <- 0
                if (model[2]=="none")
                  modt <- 0
                if (model[1]=="exponential")
                    mods <- 1
                if (model[2]=="exponential")
                  modt <- 1
                if (model[1]=="stable")
                    {
                      mods <- 2
                      if ((param[1] >1) | (param[1]<0)) stop("Stable model parameter must lie in (0,1]")
                    }
                if (model[2]=="stable")
                  {
                    modt <- 2
                    if ((param[2] >1) | (param[2]<0)) stop("Stable model parameter must lie in (0,1]")
                  }
                if (model[1]=="cauchy")
                  {
                    mods <- 3
                    if (param[1]<=0) stop("Cauchy model parmaeter must be strictly positive")
                  }
                if (model[2]=="cauchy")
                  {
                    modt <- 3
                    if (param[2]<=0) stop("Cauchy model parameter must be strictly positive")
                  }
                if (model[1]=="wave")
                    mods <- 4
                if (model[2]=="wave")
                  modt <- 4
                if (model[1]=="matern")
			{
                  mods <- 7
			if (param[3]<=0 | param[1]<=0) stop("Matern model parameters must be strictly positive")
                  }
		    if (model[2]=="matern")
			{
                  modt <- 7
			if (param[2]<=0 | param[4]<=0) stop("Matern model parameters must be strictly positive")
                  }
              }
      }
    if (!(isTRUE(separable)))
      {
        if (length(model)!=1)
          stop("for non-separable covariance functions, 'model' must be of length 1")
        if (model=="gneiting")
          {
            mod <- 5
            if (param[6]<1) stop("for Gneiting's covariance function, the sixth parameter must be greater than 1")
            if ((param[3]<=0) | (param[3]>1)) stop("for Gneiting's covariance function, the third parameter must lie in (0,1]")
            if ((param[4]<=0) | (param[4]>1)) stop("for Gneiting's covariance function, the fourth parameter must lie in (0,1]")
            if ((param[5]!=1) & (param[5]!=2) & (param[5]!=3)) stop("for Gneiting's covariance function, the fifth parameter must be 1, 2 or 3")
            if ((param[2]!=1) & (param[2]!=2) & (param[2]!=3)) stop("for Gneiting's covariance function, the second parameter must be 1, 2 or 3")
            if ((param[2]==1) & ((param[1]<0) | (param[1]>2))) stop("for Gneiting's covariance function, if the second parameter equals 1, the first parameter must lie in [0,2]") 
            if ((param[2]==2) & (param[1]<=0)) stop("for Gneiting's covariance function, if the second parameter equals 2, the first parameter must be strictly positive")            
          }
        if (model=="cesare")
          {
            mod <- 6
            if (((param[1]>2) | (param[1]<1)) | ((param[2]>2) | (param[2]<1))) stop("for De Cesare's model, the first and second parameters must lie in [1,2]")
            if (param[3]<3/2) stop("for De Cesare's model, the third parameter must be greater than 3/2")
          }
      }

    return(model=c(mods,modt,mod))
  }

.matern = function (d, scale = 1, alpha = 1, nu = 0.5) 
{
    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    d <- d * alpha
    d[d == 0] <- 1e-10
    k <- 1/((2^(nu - 1)) * gamma(nu))
    res <- scale * k * (d^nu) * besselK(d, nu)
    return(res)
}


.covst <- function(dist,times,separable=TRUE,model,param=c(1,1,1,1,1,1),sigma2=1,scale=c(1,1),plot=TRUE,nlevels=10)
{

  nt <- length(times)
  np <- length(dist)

  model <- .set.cov(separable,model,param,sigma2)
  
  gs <- array(0, dim = c(np,nt))
  storage.mode(gs) <- "double"

  gs <- .Fortran("covst",
                 (gs),
                 as.double(dist),
                 as.integer(np),
                 as.double(times),
                 as.integer(nt),
                 as.integer(model),
                 as.double(param),
                 as.double(sigma2),
                 as.double(scale))[[1]]


  if (plot==TRUE)
    {
      image(dist,times,gs,col=grey((1000:1)/1000),xlab="h",ylab="t",cex.axis=1.5,cex.lab=2,font=2)
      contour(dist,times,gs,add=T,col=4,labcex=1.5,nlevels=nlevels)
    }
  
  return(gs)

}


.gauss3D <- function(nx=100,ny=100,nt=100,xlim,ylim,tlim,separable=TRUE,model="exponential",param=c(1,1,1,1,1,1),scale=c(1,1),var.grf=1,mean.grf=0,exact=TRUE)
{

  N <- c(nx,ny,nt)

  mod <- .set.cov(separable,model,param,var.grf)
  
  g <- floor(log(2*(N-1))/log(2))+1
  M <- 2^g

  count <- 1
  changG <- TRUE
  while(changG==TRUE)
    {
      L <- rep(-9999,M[1]*M[2]*M[3])  
      
      storage.mode(L) <- "double"
      
      res <- .Fortran("circ",
                      (L),
                      as.integer(M),
                      as.integer(N),
                      as.double(xlim),
                      as.double(ylim),
                      as.double(tlim),
                      as.integer(mod), 
                      as.double(param),
                      as.double(var.grf),
                      as.double(scale))[[1]]
 
      L <- array(res,dim=M)
      FTL <- fft(L)

      if (isTRUE(exact))
        {      
          if (min(Re(FTL))<0)
            {
              g <- g+1
              M <- 2^g
              changG <- TRUE
              count <- count+1
            }
          else
            changG <- FALSE
        }
      else
        {
          FTL[Re(FTL)<0] <- 0
          changG <- FALSE
        }

    }
  print(count)
  
  X <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
  X <- fft(X,inverse=TRUE)
  A <- sqrt(FTL)*X
  G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
  G <- G+mean.grf

## Remarks
# --------
#
#  The right way is:
#    X1 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
#    X2 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
#    X <- fft(X1+1i*X2,inverse=TRUE)
#    A <- sqrt(FTL)*X
#    G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
#  but it doesn't change the results and it is slightly faster. 
#  
# Taking the first elements of FTL and then computing G (changing M by N)
# is faster than taking the first elements of G, but it does not provide
# the same results
  
  return(G)
}



rlgcp <- function(s.region, t.region, replace=TRUE, npoints=NULL, nsim=1, nx=100, ny=100, nt=100,separable=TRUE,model="exponential",param=c(1,1,1,1,1,1),scale=c(1,1),var.grf=1,mean.grf=0,lmax=NULL,discrete.time=FALSE,exact=FALSE)
{
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)
  
  t.region <- sort(t.region)
  s.area <- areapl(s.region)
  t.area <- t.region[2]-t.region[1]
  tau <- c(start=t.region[1],end=t.region[2],step=(t.region[2]-t.region[1])/(nt-1))
  bpoly <- bbox(s.region)

  lambdamax <- lmax
  pattern <- list()
  index.t <- list()
  Lambdafin <- list()
  ni <- 1

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

  while(ni<=nsim)
    {
      S <- .gauss3D(nx=nx,ny=ny,nt=nt,xlim=range(s.region[,1]),ylim=range(s.region[,2]),tlim=range(t.region),separable=separable,model=model,param=param,scale=scale,var.grf=var.grf,mean.grf=mean.grf,exact=exact)

      Lambda <- exp(S)

      mut <- rep(0,nt)
      for (it in 1:nt)
        {
          Lambda[,,it][s.grid$mask==FALSE] <- NaN
          mut[it] <- sum(Lambda[,,it],na.rm=TRUE)
        }
      
      if (is.null(npoints))
        {
          en <- sum(Lambda,na.rm=TRUE)*s.grid$xinc*s.grid$yinc*t.grid$tinc
          npoints <- round(rpois(n=1,lambda=en),0)
        }

      if (is.null(lambdamax))
        lambdamax <- max(Lambda,na.rm=TRUE)
  
      npts <- round(lambdamax/(s.area*t.area),0)
      if (npts==0) stop("there is no data to thin")
  
      if ((replace==FALSE) & (nt < max(npts,npoints))) stop("when replace=FALSE, nt must be greater than the number of points used for thinning")

      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          times.init <- sample(vect,nt,replace=replace)
        }
      else
        times.init <- runif(nt,min=t.region[1],max=t.region[2])
      
      samp <- sample(1:nt,npts,replace=replace,prob=mut/max(mut,na.rm=TRUE))
      times <- times.init[samp]

      retain.eq.F <- FALSE
      while(retain.eq.F==FALSE)
        {
          xy <- matrix(csr(poly=s.region,npoints=npts),ncol=2)
          x <- xy[,1]
          y <- xy[,2]
      
          prob <- NULL
          for(ix in 1:length(x))
            {
              nix <- findInterval(vec=s.grid$x,x=x[ix])
              niy <- findInterval(vec=s.grid$y,x=y[ix])
              nit <- findInterval(vec=t.grid$times,x=times[ix])
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
          if (sum(retain==FALSE)==length(retain)) retain.eq.F <- FALSE
          else retain.eq.F <- TRUE
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
          if (dim(xy)[2]==1) {wx <- xy[1]; wy <- xy[2]}
          else {wx <- xy[,1]; wy <- xy[,2]}
          
          if (isTRUE(replace))
            wsamp <- sample(1:nt,npoints-neffec,replace=replace,prob=mut/max(mut,na.rm=TRUE))
          else
            wsamp <- sample(samp.remain,npoints-neffec,replace=replace,prob=mut[samp.remain]/max(mut[samp.remain],na.rm=TRUE))
          
          wtimes <- times.init[wsamp]
          prob <- NULL
          for(ix in 1:length(wx))
            {
              nix <- findInterval(vec=s.grid$x,x=wx[ix])
              niy <- findInterval(vec=s.grid$y,x=wy[ix])
              nit <- findInterval(vec=t.grid$times,x=wtimes[ix])
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
      pattern.interm <- cbind(x=x,y=y,t=times)

      if (nsim==1)
        {
          pattern <- as.3dpoints(pattern.interm)
          index.t <- index.times
          Lambdafin <- Lambda
        }
      else
        {
          pattern[[ni]] <- as.3dpoints(pattern.interm)
          index.t[[ni]] <- index.times
          Lambdafin[[ni]] <- Lambda
        }
      ni <- ni+1
    }

  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region,Lambda=Lambdafin,index.t=index.t)))
}









