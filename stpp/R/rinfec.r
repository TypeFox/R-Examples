rinfec <- function(npoints, s.region, t.region, nsim=1, alpha, beta, gamma, 
s.distr="exponential", t.distr="uniform", maxrad, delta, h="step", g="min", 
recent=1, lambda=NULL, lmax=NULL, nx=100, ny=100, nt=1000, 
t0, inhibition=FALSE, ...)
{
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)
  if (missing(t0)) t0 <- t.region[1]

  if (missing(maxrad)) maxrad <- c(0.05,0.05)
  maxrads <- maxrad[1]
  maxradt <- maxrad[2]

  if (missing(delta)) delta <- maxrads

  if (s.distr=="poisson")
    {
      if (is.matrix(lambda))
        {
          nx <- dim(lambda)[1]
          ny <- dim(lambda)[2]
          Lambda <- lambda
        }
      s.grid <- .make.grid(nx,ny,s.region)
      s.grid$mask <- matrix(as.logical(s.grid$mask),nx,ny)
      if (is.function(lambda))
        {
          Lambda <- lambda(as.vector(s.grid$X),as.vector(s.grid$Y),...)
          Lambda <- matrix(Lambda,ncol=ny,byrow=T)
        }
      if (is.null(lambda))
        {
          Lambda <- matrix(npoints/areapl(s.region),ncol=ny,nrow=nx)
        }
      Lambda[s.grid$mask==FALSE] <- NaN
    }
  
  if (h=="step")
    {
      hk <- function(t,t0,alpha,beta,gamma)
        {
          res <- rep(0,length(t))
          mask <- (t <= t0+alpha+gamma) & (t >= t0+alpha)
          res[mask] <- beta
          return(res)
        }
    }

  if (h=="gaussian")
    {
      hk <- function(t,t0,alpha,beta,gamma)
        {
          t0 <- alpha+t0
          d <- gamma/8
          top <- beta*sqrt(2*pi)*d
#          res <- top*(1/(sqrt(2*pi)*d))*exp(-(((t-(t0+gamma/2))^2)/(2*d^2)))
          res <- top*dnorm(t,mean=t0+gamma/2,sd=d)
          return(res)
        }

      d <- gamma/8
      top <- beta*sqrt(2*pi)*d
      ug <- top*dnorm(0,mean=alpha+gamma/2,sd=d)
    }
  
  if (g=="min")
    {
      pk <- function(mu,recent)
        {
          if (recent=="all")
            res <- min(mu)
          else
            {
              if (is.numeric(recent))
                {
                  if(recent<=length(ITK))
                    res <- min(mu[(length(ITK)-recent+1):length(ITK)])
                  else
                    res <- min(mu)
                }
              else
                stop("'recent' must be numeric")
            }
          return(res)
        }
    }
  
  if (g=="max")
    {
      pk <- function(mu,recent)
        {
          if (recent=="all")
            res <- max(mu)
          else
            {
              if (is.numeric(recent))
                {
                  if(recent<=length(ITK))
                    res <- max(mu[(length(ITK)-recent+1):length(ITK)])
                  else
                    res <- max(mu)
                }
              else
                stop("'recent' must be numeric")
            }
          return(res)
        }
    }
    if (g=="prod")
    {
      pk <- function(mu,recent)
        {
          if (recent=="all")
            res <- prod(mu)
          else
            {
              if (is.numeric(recent))
                {
                  if(recent<=length(ITK))
                    res <- prod(mu[(length(ITK)-recent+1):length(ITK)])
                  else
                    res <- prod(mu)
                }
              else
                stop("'recent' must be numeric")
            }
          return(res)
        }
    }
  
  t.grid <- list(times=seq(t.region[1],t.region[2],length=nt),tinc=((t.region[2]-t.region[1])/(nt-1)))
  
  pattern <- list()
  ni <- 1
  while(ni<=nsim)
    {
      xy <- csr(npoints=1,poly=s.region)
      npts <- 1
      pattern.interm <- cbind(x=xy[1],y=xy[2],t=t0)
      Mu <- rep(0,nt)
      ti <- t0
      ITK <- findInterval(vec=t.grid$times,x=ti)
      
      continue <- FALSE
      while(continue==FALSE)
        {
          while (npts<npoints)
            {
              ht <- hk(t=t.grid$times,t0=ti[npts],alpha=alpha,beta=beta,gamma=gamma)
#              ht <- ht/integrate(hk,lower=t.region[1],upper=t.region[2],subdivisions=nt,t0=t0,alpha=alpha,beta=beta,gamma=gamma)$value
              ht <- ht/(sum(ht)*t.grid$tinc)
           
              mu <- Mu+ht
              if (sum(mu)==0) 
                mu2 <- mu
              else
                mu2 <- mu/max(mu)
              
              if (t.distr=="uniform")
                tk <- ti[npts] + runif(1,min=0,max=maxradt)

              if (t.distr=="exponential")
                  tk <- ti[npts] + rexp(1,rate=1/maxradt)
#                  tk <- ti[npts]+rexp(1,rate=1/(sum(mu)*t.grid$tinc))

              if (tk>=t.region[2])
                {
			N = npts
                  npts <- npoints
                  continue <- TRUE
			M = paste("only",N,"points have been generated over",npoints,". Process stopped when getting a time t greater than max(t.region): try a larger t.region or check your parameters",sep=" ")
			warning(M)
                }
              else
                {
                  itk <- findInterval(vec=t.grid$times,x=tk)

                 if ((mu[itk]==0) || ((h=="gaussian") && (mu[itk]<ug)))
                   {
			  N = npts
                    npts <- npoints
                    continue <- TRUE
			  M = paste("only",N,"points have been generated over",npoints,"(process stopped when getting mu(t)=0",sep=" ")
			  warning(M)
                   }
                 else
                   {
                     Mu <- mu
                     if (inhibition==FALSE)
                       {
                         cont <- FALSE
                         while(cont==FALSE)
                           {
	                       uk <- runif(1)
                             xpar <- pattern.interm[npts,1]
                             ypar <- pattern.interm[npts,2]

                             out <- TRUE
                             while(out==TRUE)
                               {
                                 if (s.distr=="uniform")
                                   {
                                     xp <- xpar+runif(1,min=-maxrads,max=maxrads)
                                     yp <- ypar+runif(1,min=-maxrads,max=maxrads)
                                   }
                                 if (s.distr=="gaussian")
                                   {
                                     xp <- xpar+rnorm(1,mean=0,sd=maxrads/2)
                                     yp <- ypar+rnorm(1,mean=0,sd=maxrads/2)
                                   }
                                 if (s.distr=="exponential")
                                   {
                                     xp <- xpar+sample(c(-1,1),1)*rexp(1,rate=1/maxrads)
                                     yp <- ypar+sample(c(-1,1),1)*rexp(1,rate=1/maxrads)
                                   }
                                 if (s.distr=="poisson")
                                   {
                                     if (is.null(lmax))
                                       lmax <- max(Lambda,na.rm=TRUE)
                                     retain.eq.F <- FALSE
                                     while(retain.eq.F==FALSE)
                                       {
                                         xy <- csr(poly=s.region,npoints=1)
                                         xp <- xy[1]
                                         yp <- xy[2]
                                         up <- runif(1)
                                         nix <- findInterval(vec=s.grid$x,x=xp)
                                         niy <- findInterval(vec=s.grid$y,x=yp)
						     if (nix!=0 & niy!=0)
							{
                                          Up <- Lambda[nix,niy]/lmax
                                          retain <- up <= Up
                                          if ((retain==FALSE) || is.na(retain))
                                            retain.eq.F <- FALSE
                                          else
                                            retain.eq.F <- TRUE
							}
							else 
							  retain.eq.F <- FALSE
                                       }
                                   }
                                 M <- inout(pts=rbind(c(xp,yp),c(xp,yp)),poly=s.region)
                                 if ((sqrt((xp - pattern.interm[npts,1])^2 + (yp - pattern.interm[npts,2])^2) < maxrads)) M <- c(M,M)
                                 if (sum(M)==4)
                                   out <- FALSE
                               }
                             
                             if (all(sqrt((xp - pattern.interm[,1])^2 + (yp - pattern.interm[,2])^2) < delta))
                               umax <- 1
                             else
					umax <- pk(mu=mu2[ITK],recent=recent)
#                               umax <- mu2[itk]

                             if (uk < umax)
                               {
                                 npts <- npts+1
                                 ITK <- c(ITK,itk)
                                 ti <- c(ti,tk)
                                 pattern.interm <- rbind(pattern.interm,c(xp,yp,tk))
                                 cont <- TRUE
                               }
                           }                          
                       }
                     else
                       {
				 uk <- runif(1)
                         xy <- csr(npoints=1,poly=s.region)
                         xp <- xy[1]
                         yp <- xy[2]
                         
                         if (all((sqrt((xp - pattern.interm[,1])^2 + (yp - pattern.interm[,2])^2)) > delta))
                           umax <- 1
                         else
                           umax <- pk(mu=mu2[ITK],recent=recent)
                         if (uk < umax)
                           {
                             npts <- npts+1
                             ITK <- c(ITK,itk)
                             ti <- c(ti,tk)
                             pattern.interm <- rbind(pattern.interm,c(xp,yp,tk))
                           }
                       }
                   }
               }
            }
          continue <- TRUE
        }
      dimnames(pattern.interm) <- list(NULL,c("x","y","t"))
      if (nsim==1)
        {
          pattern <- as.3dpoints(pattern.interm)
          ni <-  ni+1
        }
      else
        {
          pattern[[ni]] <- as.3dpoints(pattern.interm)
          ni <- ni+1
        }
    }
  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region)))
}

 
