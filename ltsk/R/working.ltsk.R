working.ltsk <-
function(q0,obs,th,vth,vlen,llim,verbose,Large,future)
{
 	## working function for ltsk
	fit <- rep(0,3)
	ii <- dnb(q0,obs,th,future=future)
	if( length(ii)<=5 ){
		if(verbose) cat(q0,'k= ',length(ii),'\n')
		return(c(0,0,4))
	}
	## chkres <- chknb(obs[ii,],llim)
 	## alternative check based upon discussion with Jin Aug 04
	ssout <- dsubsample(obs[ii,],Large=Large)
	nbr <- ssout$nbr
	if(verbose)
	{
		with(ssout,cat(q0,'k= ',nrow(nbr),'ns=',ns,'nt=',nt,'\n'))
	}
	if( (ssout$ns > llim[1]) && (ssout$nt > llim[2]) )
	{
		vout <- dvariogram(nbr,vth,vlen)
		vout <- dsmooth.variogram(vout)
		fout <- dfitvariogram(vout,nbr)
    if(fout$ret){
		  gout <- cal.gamma(q0,nbr,fout)
		  fit[1:2] <- with(gout,work.kriging(Gamma,gamma,dat[,4]))
		  fit[3] <- 0 ## sucess
    }else{
      fit <- c(mean(nbr[,4]),sd(nbr[,4]),4) ## variogram not fit
    }
	}else if(ssout$nt <= llim[2]){
	    if (verbose) cat('insufficient time points.\n')
		fit <- c(mean(nbr[,4]),sd(nbr[,4]),1)
	}else if(ssout$ns <= llim[1]){
		if (verbose) cat('insufficient space points.\n')		
		fit <- c(mean(nbr[,4]),sd(nbr[,4]),2)
	}else{
		if (verbose) cat('insufficient space & time points.\n')
		fit <- c(mean(nbr[,4]),sd(nbr[,4]),3)
	}
	fit
}
