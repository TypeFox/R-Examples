PCFhat <- function(xyt, s.region, t.region, dist, times, lambda, ks="box", hs, kt="box", ht, correction="isotropic") 
{
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
   dist=make.even.breaks(maxd,bstep=maxd/512)$r
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
  pcfhat <- array(0, dim = c(ndist,ntimes,5))

  frac=1	
  if (missing(hs))
	{
       d=dist(pts)
       if (ks=="gaussian") hs = dpik(d,kernel="normal",range.x=c(min(d),max(d)/frac))
	 else hs = dpik(d,kernel=ks,range.x=c(min(d),max(d)/frac))
      } 
   if (missing(ht))
	{
	 d=dist(ptst)
       if (kt=="gaussian") ht = dpik(d,kernel="normal",range.x=c(min(d),max(d)/frac))
	 else ht = dpik(d,kernel=kt,range.x=c(min(d),max(d)/frac))
	}

   kernel=c(ks=ks,hs=hs,kt=kt,ht=ht)

      if (ks=="box") ks=1 	
	else if (ks=="epanech") ks=2
	else if (ks=="gaussian") ks=3
	else if (ks=="biweight") ks=4

      if (kt=="box") kt=1 	
	else if (kt=="epanech") kt=2
	else if (kt=="gaussian") kt=3
	else if (kt=="biweight") kt=4

  if(missing(lambda))
    {
      misl <- 1
      lambda <- rep(npt/area,npt)
    }
  else misl <- 0
  if (length(lambda)==1) lambda <- rep(lambda,npt)

  storage.mode(pcfhat) <- "double"
 
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
 

  klist <- .Fortran("pcffunction", as.double(ptsx),
                    as.double(ptsy), as.double(ptst), 
                    as.integer(npt), as.double(polyx),
                    as.double(polyy), as.integer(np),
                    as.double(dist), as.integer(ndist),
                    as.double(times), as.integer(ntimes),
                    as.double(bsupt), as.double(binft),
      		  as.double(lambda), as.integer(ks), as.integer(kt),
                    as.double(hs), as.double(ht), 
			  (pcfhat), as.double(wbi),
			  as.double(wbimod), as.double(wt),
			  as.integer(correc2))
  pcfhat <- klist[[19]]

  pcfhat[,,c(1,2,5)]=pcfhat[,,c(1,2,5)]/area

  pcfhat <- pcfhat/(4*pi*dist)

  if(length(id)==1) PCFhat=as.array(pcfhat[,,id])
  else
  {
  PCFhat=list()
  for(i in 1:length(id)) PCFhat[[i]]=pcfhat[,,id[i]]	
  names(PCFhat)=correc[id]
  }
  correction=correc[id]

  invisible(return(list(pcf=PCFhat,dist=dist,times=times,kernel=kernel,correction=correction)))  
}
