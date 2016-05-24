partUtil <- function(obs,query,ncl,th,xcoord='x',ycoord='y',tcoord='t',fac=1.0)
{
  ## partition the observed and query data in 3D
  ## obs: observed data frame
  ##query: query data frame
  ##ncl: number of clusters
  ##th: threshold
  ## coordinate names;
  ## Value:
  ## query: query points divided into regions;
  ## obs: observed points divided into the same regions;
  ## fac: tuning parameters for load balancing;
  
  ## calculate temporal domains with limited overlap
  batchTime <- function(range,th,cutoff=0.3){
    ## range: the range of temporal domain
    ## th   : the threshold of temporal domain
    ## cutoff: the maximum overlap allowed
    
    max <- floor(range/(2*th))
    nbin<-seq(2,max)
    overlap <- 2*th+pmax(nbin-2,0)*2*th
    pct <- overlap/range
    if (min(pct)>=0.3){
      index <- 1
    }
    else{
      index <- max(which(pct<0.3))
    }
    nbin[index]
  }
  
  timeRange <- (range(query[,tcoord],obs[,tcoord],na.rm=T))+c(-1e-10,.001)
  lonr <- range(query[,xcoord],obs[,xcoord])+c(-1e-10,.001)
  latr <- range(query[,ycoord],obs[,ycoord])+c(-1e-10,.001)
  batch.Time <- batchTime(diff(timeRange),th[2])
  batch.Space <- ceiling(ncl*fac/batch.Time)
  batch.X <- ceiling(sqrt(batch.Space))
  batch.Y <- ceiling(sqrt(batch.Space))
  
  ## devide the query points first
  tstampb <- seq(timeRange[1],timeRange[2],len=batch.Time+1)
  tsid <- findInterval(query[,tcoord],tstampb)
  lonb <- seq(lonr[1],lonr[2],len=batch.X+1)
  latb <- seq(latr[1],latr[2],len=batch.Y+1)
  lonid <- findInterval(query[,xcoord],lonb)
  latid <- findInterval(query[,ycoord],latb)
  
  query.id <- expand.grid(x=unique(lonid),y=unique(latid),t=unique(tsid))
  
  ## obtain the observation intervals
  obsInt <- function(bins,th)
  {
    nbatch <- length(bins)-1
    tstampo <- cbind(bins[-(nbatch+1)],bins[-1])
    tstampo[-1,1]<- tstampo[-1,1]-th
    tstampo[-(nbatch),2]<- tstampo[-(nbatch),2]+th
    tstampo
  }
  
  tstampo<- obsInt(tstampb,th[2])
  latbo <- obsInt(latb,th[1])
  lonbo <- obsInt(lonb,th[1])
  
  findInt <- function(x)
  {
    ## find the observations
    ## x: interval ID in longitude
    ## y: interval ID in latitude
    ## t: interval ID in time domain
    id1 <- which(obs[,tcoord]<tstampo[x[3],2] & obs[,tcoord]>tstampo[x[3],1])
    if(length(id1)>0){
      tmp <- obs[id1,]
      id2 <-  which(tmp[,xcoord]<lonbo[x[1],2] & tmp[,xcoord]>lonbo[x[1],1])
      if(length(id2)>0){
        tmp2 <- tmp[id2,]
        id3 <- which(tmp2[,ycoord]<latbo[x[2],2] & tmp2[,ycoord]>latbo[x[2],1])
        if(length(id3)>0){r <- id1[id2[id3]]}
        else{r <- NULL}
      }
      else{
        r <- NULL
      }
    }
    else{
      r <- NULL
    }
    r
  }
  
  tolist <- function(x)
  {
    if(class(x)=="list"){
      r <- x
    }
    else if(is.null(dim(x))){
      r <- split(x,seq(along=x))
    }
    else if(length(dim(x))==2){
      r <- split(x, rep(1:ncol(x), each = nrow(x)))
    }
    else{
      stop(x," is not a matrix or vector")
    }
    r
  }
  
  query.list <- apply(query.id,1,function(x) which(lonid==x[1]&latid==x[2]&tsid==x[3]))
  query.list <- tolist(query.list)
  obs.list <- apply(query.id,1,findInt)
  obs.list <- tolist(obs.list)
  
  ## remove empty query
  nQuery <- sapply(query.list,length)
  
  list(query=query.list[nQuery>0],obs=obs.list[nQuery>0])
}

partSpUtil <- function(obs,query,ncl,th,xcoord='x',ycoord='y',fac=1.2){
  ## partition the observed and query data in 2D
  ## obs: observed data frame
  ##query: query data frame
  ##ncl: number of clusters
  ##th: threshold
  ## coordinate names;
  ## Value:
  ## query: query points divided into regions;
  ## obs: observed points divided into the same regions;
  ## fac: tuning parameters for load balancing;
  
  ## extend another dummy dimension t
  batch.Space <- ceiling(ncl*fac)
  batch.X <- ceiling(sqrt(batch.Space))
  batch.Y <- ceiling(sqrt(batch.Space))
  
  ## devide the query points first
  lonr <- range(query[,xcoord])+c(-1e-10,.001)
  latr <- range(query[,ycoord])+c(-1e-10,.001)
  lonb <- seq(lonr[1],lonr[2],len=batch.X+1)
  latb <- seq(latr[1],latr[2],len=batch.Y+1)
  lonid <- findInterval(query[,xcoord],lonb)
  latid <- findInterval(query[,ycoord],latb)
  
  query.id <- expand.grid(x=unique(lonid),y=unique(latid))
  
  ## obtain the observation intervals
  obsInt <- function(bins,th)
  {
    nbatch <- length(bins)-1
    tstampo <- cbind(bins[-(nbatch+1)],bins[-1])
    tstampo[-1,1]<- tstampo[-1,1]-th
    tstampo[-(nbatch),2]<- tstampo[-(nbatch),2]+th
    tstampo
  }
  
  latbo <- obsInt(latb,th)
  lonbo <- obsInt(lonb,th)
  
  findInt <- function(x)
  {
    ## find the observations
    ## x: interval ID in longitude
    ## y: interval ID in latitude
    id1 <- rep(1:nrow(obs))
    if(length(id1)>0){
      tmp <- obs[id1,]
      id2 <-  which(tmp[,xcoord]<lonbo[x[1],2] & tmp[,xcoord]>lonbo[x[1],1])
      if(length(id2)>0){
        tmp2 <- tmp[id2,]
        id3 <- which(tmp2[,ycoord]<latbo[x[2],2] & tmp2[,ycoord]>latbo[x[2],1])
        if(length(id3)>0){r <- id1[id2[id3]]}
        else{r <- NULL}
      }
      else{
        r <- NULL
      }
    }
    else{
      r <- NULL
    }
    r
  }

  tolist <- function(x)
  {
    if(class(x)=="list"){
      r <- x
    }
    else if(is.null(dim(x))){
      r <- split(x,seq(along=x))
    }
    else if(length(dim(x))==2){
      r <- split(x, rep(1:ncol(x), each = nrow(x)))
    }
    else{
      stop(x," is not a matrix or vector")
    }
    r
  }
  
  query.list <- apply(query.id,1,function(x) which(lonid==x[1]&latid==x[2]))
  query.list <- tolist(query.list)
  obs.list <- apply(query.id,1,findInt)
  obs.list <- tolist(obs.list)
  
  ## remove empty query
  nQuery <- sapply(query.list,length)
  
  list(query=query.list[nQuery>0],obs=obs.list[nQuery>0])
}