dvariogram <-
function(obs, vth=NULL,vlen = NULL)
{
  z <- matrix(obs[,4],ncol=1)
  G <- (dist(z))^2/2
  tstamp <- matrix(obs[,3],ncol=1)
  coords <- matrix(obs[,1:2],ncol=2)
  dx <- diff(range(coords[,1]))
  dy <- diff(range(coords[,2]))
  dmax <- sqrt(dx^2 + dy^2)
  tmax <- diff(range(tstamp))
  if(is.null(vth)) vth <- c(dmax,tmax)*.75
  if(is.null(vlen)) vlen <- c(15,floor(tmax))
  dmat <- dist(coords)
  tmat <- dist(tstamp)
  dbins <- seq(-1e-5,vth[1]+0.01,len=vlen[1]+1)
  tbins <- seq(-0.01,vth[2]+0.01,len=vlen[2]+1)
  dii <- findInterval(dmat,dbins) ## all.inside =F ignore > maxdist
  tii <- findInterval(tmat,tbins)
  tmp <- aggregate(x=as.vector(G),by=list(dbins=(dii),tbins=(tii)),FUN=mean)
  ddist <- dbins[-1]
  tdist <- tbins[-1]
  gamma.long <- data.frame(s=ddist[tmp[,1]],t=tdist[tmp[,2]],gamma=tmp[,3])
  vmm <- matrix(NA,ncol=length(dbins),nrow=length(tbins))
  colnames(vmm) <- paste('(',round(dbins,2),
	',',c(round(dbins[-1],2),999),')',sep='')
  row.names(vmm) <- paste('(',round(tbins,2),
	',',c(round(tbins[-1],2),999),')',sep='')
  vmm[as.matrix(tmp[,2:1])] <- tmp[,3]
  nmm <- matrix(0,ncol=length(dbins),nrow=length(tbins))
  tmp <- aggregate(x=as.vector(G), 
	by=list(dbins=as.vector(dii),tbins=as.vector(tii)), 
	FUN=length)
  gamma.long$n <- tmp[,3]
  nmm[as.matrix(tmp[,2:1])] <- tmp[,3]
  colnames(nmm) <- paste('(',round(dbins,2),
	',',c(round(dbins[-1],2),999),')',sep='')
  row.names(nmm) <- paste('(',round(tbins,2),
	',',c(round(tbins[-1],2),999),')',sep='')
  vmm <- vmm[-nrow(vmm),-ncol(vmm)]
  nmm <- nmm[-nrow(nmm),-ncol(nmm)]
  ii <- with(gamma.long, which(is.na(s) | is.na(t)))
  gamma.s0 <- data.frame(x=ddist,n=nmm[1,],gamma=vmm[1,])
  gamma.0t <- data.frame(x=tdist,n=nmm[,1],gamma=vmm[,1])
  list(gamma=vmm,gamma.s0=gamma.s0,gamma.0t=gamma.0t,gamma.long=gamma.long[-ii,],np = nmm, ddist=ddist,tdist=tdist)
}
