dsubsample <-
function(obs,nbin=NULL,Large=2000)
{
  ## obs : coords and timestamps of observed points
  ## nbin : number of lon,lat,day to bin data
  ## Large : upper limit of neighbor points
  ## value : subsample data if large, otherwise return same
  ##       : according to number of points in each strata
  ##       : whether time space kriging is appropriate
  tmp <- if( is.null(nrow(obs))) matrix(obs,1,4) else obs
  colnames(tmp) <- c('x','y','t','z')
  ## if(nrow(tmp) <= Large) return(tmp)
  lonr <- range(tmp[,1])+c(-1e-5,.001)
  latr <- range(tmp[,2])+c(-1e-5,.001)
  tstampr <- range(tmp[,3])+c(-.001,.001)
  ##if(is.null(nbin)) nbin <- c(10,10,length(unique(tmp[,3])))
  if(is.null(nbin)) nbin <- c(10,10,floor(diff(tstampr))+1)
  lonb <- seq(lonr[1],lonr[2],len=nbin[1]+1)
  latb <- seq(latr[1],latr[2],len=nbin[2]+1)
  tstampb <- seq(tstampr[1],tstampr[2],len=nbin[3]+1)
  lonid <- findInterval(tmp[,1],lonb)
  latid <- findInterval(tmp[,2],latb)
  tsid <- findInterval(tmp[,3],tstampb)
  uid <- lonid*100000+latid*1000 + tsid
  nt <- length(unique(tsid))
  ns <- length(unique(lonid*1000+latid))
  if(nrow(tmp) <= Large)
  {
	 out <- list(nbr=tmp,ns=ns,nt=nt)
	 return(out)
  }
  tmp2 <- aggregate(tmp[,4],by=list(x=lonid,y=latid,t=tsid),length)
  S <- nrow(tmp2)
  Ns <- tmp2[,4]
  nn <- min(Large,nrow(tmp))
  if( nn > S){
	nns <- dsample.pps(Ns-1, nn-S) + 1
  }else{
	nns <- dsample.pps(Ns,nn)
  }
  ids <- dsample.strata(uid,Ns,nns)
  list(nbr = tmp[ids,], ns =ns,nt=nt)
}
