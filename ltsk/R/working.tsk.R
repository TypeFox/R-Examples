working.tsk <-
function(q0,obs,subset,nmin,nmax,vth,vlen,llim,verbose,Large,future)
{
 	## working function for tsk
	fit <- cbind(rep(mean(obs[,4]),nrow(q0)),sd(obs[,4]),4)
  vout <- NULL
  fout <- NULL
	if( nrow(obs)<=5 ){
		if(verbose) cat('k= ',nrow(obs),'\n')
    r <- list(krig=fit,variog=NULL,fitvariog=NULL)
		return(r)
	}
	## chkres <- chknb(obs[ii,],llim)
 	## alternative check based upon discussion with Jin Aug 04
	ssout <- dsubsample(obs,Large=Large)
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
		  kout <- work.kriging.vec(query=q0,obs=nbr,fout=fout,subset=subset,nmin=nmin,nmax=nmax,future=future,verbose=verbose)
      fit <- cbind(kout,0)
    }else{
      if (verbose) cat("variogram not fitted\n")
      fit[,3] <- 4
      #fit <- c(mean(nbr[,4]),sd(nbr[,4]),4) ## variogram not fit
    }
	}else if(ssout$nt <= llim[2]){
	  if (verbose) cat('insufficient time points.\n')
    fit[,3] <- 1  
		#fit <- c(mean(nbr[,4]),sd(nbr[,4]),1)
	}else if(ssout$ns <= llim[1]){
		if (verbose) cat('insufficient space points.\n')
    fit[,3] <- 2
		#fit <- c(mean(nbr[,4]),sd(nbr[,4]),2)
	}else{
		if (verbose) cat('insufficient space & time points.\n')
    fit[,3] <- 3
		#fit <- c(mean(nbr[,4]),sd(nbr[,4]),3)
	}
	list(krig=fit,variog=vout,fitvariog=fout)
}
