.spatial.inhibition <- function(npoints,h,theta,delta,p,recent="all",s.region,inhibition=TRUE,...)
  {
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)

  if (!(is.function(h)))
    {
      models <- c("step","gaussian")
      if (sum(h==models)==0)
        {
          message <- paste("this model is not implemented, please choose among: ",paste(models,"",sep=" ",collapse="and "))
          stop(message)
        }
      if (h=="step")
        {
          hk <- function(d,theta,delta)
            {
              res <- rep(1,length(d))
		  if (inhibition==TRUE) res[d<=delta] <- theta
		  else res[d>=delta] <- theta
              return(res)
            }
        }
      if (h=="gaussian")
        {
          hk <- function(d,theta,delta)
            {
              if (inhibition==TRUE) 
			{
			res=NULL
			for(i in 1:length(d))
				{	
				if (d[i]<=delta) res=c(res,0)
				if (d[i]>(delta+theta/2)) res=c(res,1)
				if (d[i]>delta & d[i]<=(delta+theta/2)) res=c(res,exp(-((d[i]-delta-theta/2)^2)/(2*(theta/8)^2)))
				}
			}
		  else
			{
			res=NULL
			for(i in 1:length(d))
				{	
				if (d[i]<delta) res=c(res,1)
				else res=c(res,exp(-((d[i]-delta)^2)/(2*(theta/8)^2)))
				}
			}
	   	  return(res)
		}
	  }
	}
 else
   	{   
          hk <- function(d,theta,delta)
            {
            res <- h(d,theta,delta)
            return(res)
		}
	}

	pk <- function(d,h,recent,theta,delta)
        {
          if (recent=="all")
		{
		 if (p=="min") res <- min(h(d=d,theta=theta,delta=delta))
		 if (p=="max") res <- max(h(d=d,theta=theta,delta=delta))
		 if (p=="prod") res <- prod(h(d=d,theta=theta,delta=delta))
		}
          else
            {
              if (is.numeric(recent))
                {
                  if(recent<=length(d))
				{
				if (p=="min") res <- min(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				if (p=="max") res <- max(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				if (p=="prod") res <- prod(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				}
                  else
				{
				if (p=="min") res <- min(h(d=d,theta=theta,delta=delta))
				if (p=="max") res <- max(h(d=d,theta=theta,delta=delta))
				if (p=="prod") res <- prod(h(d=d,theta=theta,delta=delta))
				}
                }
              else
                stop("'recent' must be numeric")
            }
          return(res)
        }

  xy <- csr(npoints=1,poly=s.region)
  npts <- 1
  pattern.interm <- cbind(x=xy[1],y=xy[2])
  if (inhibition==TRUE)
    {
      while(npts < npoints)
        {
          prob <- runif(1)
          xy <- csr(npoints=1,poly=s.region)
          if (all((sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2)) > delta))
            umax <- 1
            else			
             umax <- pk(d=sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2),hk,recent,theta,delta)             
		
          if (prob<umax)
            {
              pattern.interm <- rbind(pattern.interm,c(xy[1],xy[2]))
              npts <- npts+1
            }
        }
    }
  else
    {
      while(npts < npoints)
        {
          prob <- runif(1)
          continue <- FALSE
          while(continue==FALSE)
            {
              xy <- csr(npoints=1,poly=s.region)
              if (all(sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2) < delta))
                umax <- 1            
              else		
                umax <- pk(d=sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2),hk,recent,theta,delta)             	
                
              if (prob < umax)
                {
                  pattern.interm <- rbind(pattern.interm,c(xy[1],xy[2]))
                  npts <- npts+1
                  continue <- TRUE
                }
            }
        }
    }
  invisible(return(list(pts=pattern.interm,s.region=s.region)))
}


.temporal.inhibition <- function(npoints,h,theta,delta,p,recent="all",t.region,discrete.time=FALSE,replace=FALSE,inhibition=TRUE)
  {  
  if (missing(t.region)) t.region <- c(0,1)

  if (!(is.function(h)))
	{
      models <- c("step","gaussian")
      if (sum(h==models)==0)
       	{
          	message <- paste("this model is not implemented, please choose among: ",paste(models,"",sep=" ",collapse="and "))
          	stop(message)
        	}
      if (h=="step")
      	{
          	hk <- function(d,theta,delta)
            	{
              	res <- rep(1,length(d))
		  	if (inhibition==TRUE) res[d<=delta] <- theta
		  	else res[d>=delta] <- theta
              	return(res)
            	}
        	}
      if (h=="gaussian")
      	{
          	hk <- function(d,theta,delta)
            	{
              	if (inhibition==TRUE) 
				{
				res=NULL
				for(i in 1:length(d))
					{	
					if (d[i]<=delta) res=c(res,0)
					if (d[i]>(delta+theta/2)) res=c(res,1)
					if (d[i]>delta & d[i]<=(delta+theta/2)) res=c(res,exp(-((d[i]-delta-theta/2)^2)/(2*(theta/8)^2)))
					}
				}
		  	else
				{
				res=NULL
				for(i in 1:length(d))
					{	
					if (d[i]<delta) res=c(res,1)
					else res=c(res,exp(-((d[i]-delta)^2)/(2*(theta/8)^2)))
					}
				}
	   	  	return(res)
			}
	  	}
	}
  else
	{
       hk <- function(d,theta,delta)
            {
            res <- h(d,theta,delta)
            return(res)
		}
    	}

  pk <- function(d,h,recent,theta,delta)
 	{
      if (recent=="all")
		{
		if (p=="min") res <- min(h(d=d,theta=theta,delta=delta))
		if (p=="max") res <- max(h(d=d,theta=theta,delta=delta))
		if (p=="prod") res <- prod(h(d=d,theta=theta,delta=delta))
		}
      else
            {
            if (is.numeric(recent))
			{
                  if(recent<=length(d))
				{
				if (p=="min") res <- min(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				if (p=="max") res <- max(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				if (p=="prod") res <- prod(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				}
                  else
				{
				if (p=="min") res <- min(h(d=d,theta=theta,delta=delta))
				if (p=="max") res <- max(h(d=d,theta=theta,delta=delta))
				if (p=="prod") res <- prod(h(d=d,theta=theta,delta=delta))
				}
                	}
		 else stop("'recent' must be numeric")
            }
	return(res)
	}
 
  if (discrete.time==FALSE)
    ti <- runif(1,min=t.region[1],max=t.region[1]+delta)
  else
    ti <- sample(floor(t.region[1]):ceiling(t.region[1]+delta),1)
  times <-  ti
  npts <- 1
  if (inhibition==TRUE)
    {
      while(npts < npoints)
        {
          if (discrete.time==FALSE)
            ti <- runif(1,min=t.region[1],max=t.region[2])
          else
            ti <- sample(floor(t.region[1]):ceiling(t.region[2]),1)

          prob <- runif(1)
          
          if (all(abs(ti - times) > delta))
            umax <- 1
          else
            umax <- pk(d=abs(ti - times),hk,recent,theta,delta)
          if (prob<umax)
            {
              times <- c(times,ti)
              npts <- npts+1
            }
        }
    }
  else
    {
      while(npts < npoints)
        {
          prob <- runif(1)
          
          continue <- FALSE
          while(continue==FALSE)
            {
              if (discrete.time==FALSE)
                ti <- runif(1,min=t.region[1],max=t.region[2])
              else
                ti <- sample(floor(t.region[1]):ceiling(t.region[2]),1)
              
              if (abs(ti - times[npts]) < delta)
                umax <- 1            
              else
                umax <- pk(d=abs(ti - times),hk,recent,theta,delta)
              if (prob < umax)
                {
                  times <- c(times,ti)
                  npts <- npts+1
                  continue <- TRUE
                }
            }
        }
    }

  samp <- sample(1:npoints,npoints,replace=replace)
  times <- sort(times[samp])
  
  invisible(return(list(times=times,t.region=t.region)))
}



rinter <- function(npoints,s.region,t.region,hs="step",gs="min",thetas=0,deltas,ht="step",gt="min",thetat=1,deltat,recent="all",nsim=1,discrete.time=FALSE,replace=FALSE,inhibition=TRUE,...)
  {
  
  if (is.null(npoints)) stop("please specify the number of points to generate")
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  if (!(is.function(hs)))
    {
      if ((inhibition==TRUE) & (thetas==0) & (hs=="step") & (npoints * pi * deltas^2/4 > areapl(s.region)))
        stop(paste("s.region is too small to fit", npoints, "points", "at minimum distance", deltas))
    
      if ((inhibition==TRUE) & (thetas==0) & (hs=="step") & ((max(t.region)-min(t.region))/deltat<npoints))
        stop(paste("t.region is too small to fit", npoints, "points", "at minimum time interval", deltat))
    }

  pattern <- list()
  ni <- 1
  while(ni<=nsim)
    {
      pattern.interm <- .spatial.inhibition(npoints,h=hs,theta=thetas,delta=deltas,p=gs,recent=recent,s.region=s.region,inhibition=inhibition,...)$pts
      times.interm <- .temporal.inhibition(npoints,h=ht,theta=thetat,delta=deltat,p=gt,recent=recent,t.region=t.region,inhibition=inhibition,discrete.time=discrete.time,replace=replace,...)$times
      
      if (nsim==1)
        {
          pattern <- as.3dpoints(cbind(pattern.interm,times.interm))
          ni <-  ni+1
        }
      else
        {
          pattern[[ni]] <- as.3dpoints(cbind(pattern.interm,times.interm))
          ni <- ni+1
        }
    }
  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region)))
}




