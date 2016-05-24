dfitvariogram <-
function(vout,nbr)
{
  ## vout : output from dvariogram
  ## nbr  : orginal observation points
  ## value: obj use to calculate fitted variogram for space time kriging
  fs0 <- working.fitvariog1(vout$gamma.s0)
  vfs0 <- get(paste('v',fs0$model,sep=''))
  f0t <- working.fitvariog1(vout$gamma.0t)
  vf0t <- get(paste('v',f0t$model,sep=''))
  chk <- all(fs0$ret,f0t$ret)
  if (!chk) {
	cat('[space or time variogram fit un-success]\n')
	return(list(scoef=NULL,tcoef=NULL,smodel='null',
		tmodel='null',k=0,ret=F))
  }
  ## sill0 <- with(vout$gamma.long,max(gamma))
  sill0 <- var(nbr[,4])
  sills0 <- with(fs0,sum(coef[c(1,3)]))
  sill0t <- with(f0t,sum(coef[c(1,3)]))
  out <- dadjustsills(sill0,sills0,sill0t)
  sill0 <- out[1]; sills0 <- out[2]; sill0t <- out[3]
  stn.ratio.s0 <- max(1e-3,fs0$coef[1]/sum(fs0$coef[c(1,3)]))
  stn.ratio.0t <- max(1e-3,f0t$coef[1]/sum(f0t$coef[c(1,3)]))
  fs0$coef[1] <- sills0 * stn.ratio.s0 
  f0t$coef[1] <- sill0t * stn.ratio.0t
  fs0$coef[3] <- sills0 * (1 - stn.ratio.s0)
  f0t$coef[3] <- sill0t * (1- stn.ratio.0t)
  k1 <- (sills0 + sill0t - sill0)/ (sills0*sill0t)
  k2 <- (sill0 - sill0t) / sills0
  k3 <- (sill0 - sills0) / sill0t
  list(scoef=fs0$coef,tcoef=f0t$coef,smodel=fs0$model,
	tmodel=f0t$model,k=k1,ks=c(k1,k2,k3),sill0=sill0,ret=k1>0)
}
