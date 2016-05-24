## version 2
## include local Kriging arguments to avoid over smoothing
require(Matrix)
work.blk.kriging2 <-
function(query,obs,fout,future=T,verbose=T,nmin=0,nmax=Inf,th=c(Inf,Inf))
{
  ## query : query points (matrix)
  ## obs : neighbor points with a given thresholds
  ## fout : output fitted variogram
  ## future : whether include future observations to prediction
  ## verbose : whether output detailed debug messages
  ## nmin : the number of nearest neighbors within spatial temporal thresholds, if less than nmin, a missing value would be generated
  ## nmax : for local Kriging, the number of nearest neighbors that should be used for Kriging
  ## th: for local Kriging, only observations within this thresholds from prediction locations are used
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

  work.blk.kriging.nbr <- function(obs,query,nmax,th,fout){
    ## find a list of observations to each query point within thresholds 
    ## at most nmax neighbors
    ## arguments:
    ## obs : the observed x y and t coordinates
    ## query: the query x y and t coordinates
    ## nmax : maximum number of neighbors to return
    ## th : space and time thresholds
    ## fout: fitted product sum variogram object
    ## Value:
    ## a sparse matrix object of class Matrix
    ## rows indicate query points, column the observations 
    ## and non-zero entries denote the semivariance between the
    ## query and the observed points
    ## Details
    ## first use FNN and nmax to find spatial neighbors, expand nmax to allow for temporal constraints
    ## then apply the space and time constraints according to th
    ## then calculate the semivariance and order the results according to nmax, 
    ## if more than nmax retaining the nmax values only
    return(0)
  }
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
