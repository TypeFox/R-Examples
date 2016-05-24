STIKhat <- function(xyt, s.region, t.region, dist, times, lambda, correction="isotropic", infectious = FALSE) 
{

  if (infectious==TRUE && any(correction!= "isotropic")) stop("For infectious processes, the STIK-function is currently only implemented with the isotropic edge correction method.")	

  correc=c("none","isotropic","border","modified.border","translate")
  id <- match(correction, correc, nomatch = NA)
  if (any(nbg <- is.na(id))) {
        mess <- paste("unrecognised correction method:", paste(dQuote(correction[nbg]), 
            collapse = ", "))
        stop(mess, call. = FALSE)
    }
  id=unique(id)	
  correc2=rep(0,5)
  correc2[id]=1	

  if (missing(s.region)) s.region <- sbox(xyt[,1:2],xfrac=0.01,yfrac=0.01)
  if (missing(t.region)) 
	{
      xr = range(xyt[,3],na.rm=TRUE)
      xw = diff(xr)
      t.region <- c(xr[1]-0.01*xw,xr[2]+0.01*xw)
	}    
  bsupt <- max(t.region)
  binft <- min(t.region)
  bdry=owin(poly=list(x=s.region[,1],y=s.region[,2]))

  if (missing(dist))
  {
   rect=as.rectangle(bdry)
   maxd=min(diff(rect$xrange),diff(rect$yrange))/4
   dist=make.even.breaks(maxd,npos=15)$r
  } 
  if (missing(times))
  {
  maxt=(bsupt-binft)/4
  times=make.even.breaks(maxt,npos=15)$r
   }

   dist <- sort(dist)
   if(dist[1]==0) dist=dist[-1]
   times <- sort(times)
   if(times[1]==0) times=times[-1]

  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[, 1]
  ptsy <- pts[, 2]
  ptst <- xytimes
  npt <- length(ptsx)
  ndist <- length(dist)
  ntimes <- length(times)

  area <- areapl(s.region)*(bsupt-binft)

  np <- length(s.region[, 1])
  polyx <- c(s.region[, 1], s.region[1, 1])
  polyy <- c(s.region[, 2], s.region[1, 2])
  hkhat <- array(0, dim = c(ndist,ntimes,5))

  if(missing(lambda))
    {
      misl <- 1
      lambda <- rep(npt/area,npt)
    }
  else misl <- 0
  if (length(lambda)==1) lambda <- rep(lambda,npt)

  if (infectious==TRUE){infd <- 1} else infd <- 0
  storage.mode(hkhat) <- "double"

  wbi=array(0,dim=c(npt,ndist,ntimes))
  wbimod=array(0,dim=c(npt,ndist,ntimes))
  wt = array(0,dim=c(npt,npt)) 

  pppxy = ppp(x=ptsx,y=ptsy,window=bdry)

#  correction=="border" and "modified border"

   if(any(correction=="border")|any(correction=="modified.border"))
   {
   bi=bdist.points(pppxy)
   bj=.bdist.times(xytimes,t.region)

   for(i in 1:ndist) 
	{ 
   for(j in 1:ntimes)
      {
	  wbi[,i,j] = (bi>dist[i])*(bj>times[j])/sum((bi>dist[i])*(bj>times[j])/lambda)
        wbimod[,i,j] = (bi>dist[i])*(bj>times[j])/(eroded.areas(bdry,dist[i])*.eroded.areat(t.region,times[j]))
  	} }
	wbi[is.na(wbi)]=0
   }

# correction=="translate"
 	
  if(any(correction=="translate"))
  {
  wtt = .overlap.tint(xytimes,t.region)
  wts = edge.Trans(pppxy)
  wt = wtt*wts
  wt=1/wt
  }
	

  klist <- .Fortran("stikfunction", as.double(ptsx),
                    as.double(ptsy), as.double(ptst), 
                    as.integer(npt), as.double(polyx),
                    as.double(polyy), as.integer(np),
                    as.double(dist), as.integer(ndist),
                    as.double(times), as.integer(ntimes),
                    as.double(bsupt), as.double(binft),
                    as.double(lambda), as.integer(infd),
                    (hkhat), as.double(wbi),
			  as.double(wbimod), as.double(wt),
			  as.integer(correc2))
  hkhat <- klist[[16]]


  hkhat[,,c(1,2,5)]=hkhat[,,c(1,2,5)]/area
    	
  Khpp <- matrix(0,ncol=length(times),nrow=length(dist))
  for(i in 1:length(dist))Khpp[i,]<-pi*(dist[i]^2)*times
  
  if (infectious==TRUE) {Ktheo <- Khpp}
  else Ktheo <- 2*Khpp

  if(length(id)==1) Khat=as.array(hkhat[,,id])
  else
  {
  Khat=list()
  for(i in 1:length(id)) Khat[[i]]=hkhat[,,id[i]]	
  names(Khat)=correc[id]
  }
  correction=correc[id]

  invisible(return(list(Khat=Khat,Ktheo=Ktheo,dist=dist,times=times,correction=correction,infectious=infectious)))  
}


.overlap.tint=function(times,t.region)
{
 if (missing(t.region)) t.region=range(times)
 ntimes=length(times)

 a = diff(range(t.region))
 
 wt=matrix(a,ncol=ntimes,nrow=ntimes)
 for(i in 1:ntimes)
 { for(j in 1:ntimes)
 {
 if (i!=j)
 {
 b = a-abs(times[i]-times[j])
 wt[i,j]=a/b
 }}}
 invisible(return(wt))
}

.bdist.times=function(times, t.region)
{
 if (missing(t.region)) t.region=range(times)
 ntimes=length(times)
 a=min(t.region)
 b=max(t.region)

 bj=NULL
 for(j in 1:ntimes)
 bj=c(bj,min(c(abs(times[j]-a),abs(times[j]-b))))
 
 invisible(return(bj))
}

.eroded.areat=function(t.region,dist)
{
 a = diff(range(t.region))
 b = a-dist
}



