work.kriging.vec <-
function(query,obs,fout,subset=T,nmin=2,nmax=10,future=T,verbose=T)
{
  ## query : query points (matrix)
  ## obs : neighbor points with a given thresholds
  ## fout : output fitted variogram
  ## subset : subset neighbors using estimated sills
  ## nmin   : for local kriging: the minimal number of neighbors
  ## nmax   : the number of neighbors used for prediction
  ## future : whether include future observations to prediction
  ## value : ordinary kriging esimtates
  ## product sum model
  loc0 <- matrix(query[,1:2],ncol=2)
  tstamp <- matrix(obs[,3],ncol=1)
  coords <- matrix(obs[,1:2],ncol=2)
  z <- obs[,4]
  
  ## return
  out <- matrix(NA,nrow(query),2)
  
  
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


  ## calculate semivariance between query and neighbors
  dmat <- as.matrix(rdist(loc0,coords))
  tmat0 <- outer(query[,3],obs[,3],FUN="-")
  tmat <- abs(tmat0)
  
  gd0 <- with(fout,vfs0(as.vector(dmat),scoef[1],scoef[2],scoef[3]))
  g0t <- with(fout,vf0t(as.vector(tmat),tcoef[1],tcoef[2],tcoef[3]))
  fitgamma <- gd0 + g0t - fout$k * gd0 * g0t
  dim(fitgamma) <- dim(dmat)
  
  sel <- matrix(TRUE,nrow(query),nrow(obs)) ## initialize neighbor selections
  
  ## use sills as maximum distances
  if (subset){
    ii <- (dmat <= fout$scoef[2]) | (tmat <= fout$tcoef[2])
    sel[!ii] <- FALSE
  }
  
  ## remove future observations
  if (!future){
    sel[tmat0<0] <- FALSE
  }
  
  ## select nmax closest neighbors
  reduce <- function(x,nmax)
  {
    ii <- which(!is.na(x))
    oo <- order(x[ii])
    ii2 <- ii[oo]
    jj <- ii2[-seq(1,nmax)]
    x[jj] <- NA
    x
  }
  if(!is.null(nmax)){
    nnbr <- rowSums(sel)
    tmp <- ifelse(sel,fitgamma,NA)
    chk <- sum(nnbr>nmin)
    if(chk==0){
      return(out)
    }
    else if(chk==1){
      tmp2 <- reduce(tmp[nnbr>nmin,],nmax=nmax)
    }
    else{
      tmp2 <- t(apply(tmp[nnbr>nmin,],1,reduce,nmax=nmax))
    }
    tmp[nnbr>nmin,] <- tmp2
    sel <- !is.na(tmp)      
  }

  ## calculate kriging
  nnbr <- rowSums(sel)
  chk <- sum(nnbr>nmin)
  tmp <- ifelse(sel[nnbr>nmin,],fitgamma[nnbr>nmin,],NA)
  krig <- function(x){
    nid <- which(!is.na(x))
    l.gamma <- x[nid]
    work.kriging(fitGamma[nid,nid],l.gamma,z[nid])   
  }
  if(chk==0){
    return(out)
  }
  else if(chk==1){
    tmp2 <-krig(tmp)
  }
  else{
    tmp2 <- t(apply(tmp,1,krig))
  }  

  out[nnbr>nmin,] <- tmp2
  out
}
