tsbk <-
function(query,obs,xcoord='x',ycoord='y',tcoord='t',zcoord='z',bcoord='block',gcoord='g',
		vth=NULL,vlen=NULL,llim=c(3,3),verbose=T,Large=2000,future=T)
{

  # check input  
  l.query <- check_input(query,xcoord,ycoord,tcoord,zcoord)
  l.query <- check_na(l.query[,c(xcoord,ycoord,tcoord)],'query')
  l.obs <- check_input(obs,xcoord,ycoord,tcoord,zcoord)
  l.obs <- check_na(l.obs,'observed')
  
  l.query.coord <- na.omit(as.matrix(query[,c(gcoord,bcoord)]))
  block.lst <- table(l.query.coord[,2])
  n.block <- length(block.lst)
  
  ## return values
  fit <- data.frame(block.lst,rep(mean(l.obs[,4]),n.block),sd(l.obs[,4]),4)
  names(fit) <- c(bcoord,'Freq','krig','se','flag')
  if( nrow(l.obs)<=5 ){
      if(verbose) cat('k= ',nrow(l.obs),'\n')
      return(fit)
  }
  ssout <- dsubsample(l.obs,Large=Large)
  nbr <- ssout$nbr
  if(verbose)
  {
    with(ssout,cat('k= ',nrow(nbr),'ns=',ns,'nt=',nt,'\n'))
  }
  if( (ssout$ns > llim[1]) && (ssout$nt > llim[2]) )
  {
    vout <- dvariogram(nbr,vth,vlen)
    vout <- dsmooth.variogram(vout)
    fout <- dfitvariogram(vout,nbr)
    if(fout$ret){
      kout <- work.blk.kriging(query=cbind(l.query,l.query.coord),
                              obs=nbr,fout=fout,future=future,
                              verbose=verbose)
      fit[,1:4] <- kout
      fit[,5] <- 0
    }else{
      if (verbose) cat("variogram not fitted\n")
      fit[,5] <- 4
      #fit <- c(mean(nbr[,4]),sd(nbr[,4]),4) ## variogram not fit
    }
  }
  else if(ssout$nt <= llim[2]){
   if (verbose) cat('insufficient time points.\n')
   fit[,5] <- 1  
   #fit <- c(mean(nbr[,4]),sd(nbr[,4]),1)
  }
  else if(ssout$ns <= llim[1]){
   if (verbose) cat('insufficient space points.\n')
   fit[,5] <- 2
   #fit <- c(mean(nbr[,4]),sd(nbr[,4]),2)
  }else{
    if (verbose) cat('insufficient space & time points.\n')
    fit[,5] <- 3
    #fit <- c(mean(nbr[,4]),sd(nbr[,4]),3)
  }
  fit
}
