work.blk.kriging <-
function(query,obs,fout,future=T,verbose=T)
{
  ## query : query points (matrix)
  ## obs : neighbor points with a given thresholds
  ## fout : output fitted variogram
  ## future : whether include future observations to prediction
  ## value : ordinary block kriging esimtates
  ## product sum model
  
  tstamp <- matrix(obs[,3],ncol=1)
  coords <- matrix(obs[,1:2],ncol=2)
  z <- obs[,4]  
  
  ## calculate semi-variance between neighbors
  vfs0 <- get(paste('v',fout$smodel,sep=''))
  vf0t <- get(paste('v',fout$tmodel,sep=''))

  Dmat <- dist(coords)
  Tmat <- dist(tstamp)
  
  Gd0 <- with(fout,as.vector(vfs0(Dmat,scoef[1],scoef[2],scoef[3])))
  G0t <- with(fout,as.vector(vf0t(Tmat,tcoef[1],tcoef[2],tcoef[3])))
  tmp <- Gd0 + G0t - fout$k * Gd0 * G0t

  ## no spatial temporal structure, add small nuggets
  chk <- with(fout,(scoef[3]==0) & (tcoef[3]==0))
  if(chk){
    ## zero nuggets add ad-hoc small nuggets effect
    tmp <- tmp + 1e-5
  }  
  
  fitGamma <- tritomat(tmp,length(z))
  ## call it nxn where n is sample size;

  ## calculate semivariance between query and neighbors
  s0 <- query[,4] ## geographic coordinate ID
  t0 <- query[,3] ## time stamp ID
  
  # extract coordinates;
  loc0 <- matrix(query[!duplicated(s0),1:2],ncol=2) ## call it  s x1
  s.map <- match(s0,unique(s0))
  tstamp0 <- query[!duplicated(t0),3]  ## call it  t x1
  t.map <- match(t0,unique(t0))
  
  b0 <- query[,5] ## block ID

  dmat <- as.matrix(rdist(loc0,coords))  ## s x n
  tmat0 <- outer(tstamp0,obs[,3],FUN="-")
  tmat <- abs(tmat0) ## t x n
  
  gd0 <- with(fout,vfs0(dmat,scoef[1],scoef[2],scoef[3]))  ## s x n
  g0t <- with(fout,vf0t(tmat,tcoef[1],tcoef[2],tcoef[3]))  ## t x n
  
  # loc1 <- matrix(query[,1:2],ncol=2)
  # tstamp1 <- query[,3]
  # dmat <- as.matrix(rdist(loc1,coords))  
  # tmat <- abs(outer(tstamp1,obs[,3],FUN="-"))
  # gd1 <- with(fout,vfs0(dmat,scoef[1],scoef[2],scoef[3]))  
  # g1t <- with(fout,vf0t(tmat,tcoef[1],tcoef[2],tcoef[3])) 
  #fitgamma <- gd1 + g1t - fout$k * gd1 * g1t
  #dim(fitgamma) <- dim(dmat)
  # tmp2 <- rowsum(fitgamma,b0) 

  
  
  ## aggregation;
  tmp1 <- rowsum(gd0[s.map,]+g0t[t.map,]-fout$k*gd0[s.map,]*g0t[t.map,],b0) 
  # call it bxn where m is grid size;
  r <- data.frame(table(b0))
  fitgamma0 <- t(tmp1/r$Freq)
  
  ## call it nxb where b is block size
  
  ## calculate block kriging
  G <- rbind(cbind(fitGamma,1),c(rep(1,ncol(fitGamma)),0)) ## (n+1) x (n+1)
  g <- rbind(fitgamma0,1) # (n+1)xb
  lambda <- solve(G,g)  # (n+1) x b
  tmp <- lambda[-nrow(G),] # n x b
  m <- lambda[nrow(G),] # vector b
  zhat <- crossprod(tmp,z) # b X 1
  sigmasq <- m + apply(tmp * fitgamma0,2,sum)
  sigmasq <- pmax(0,sigmasq)  # ordinary kriging may be negative
  
  r$krig <- as.vector(zhat)
  r$se <- sqrt(sigmasq)
  r
  
}
