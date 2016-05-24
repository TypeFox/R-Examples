working.cltsk <-
function(q0,obs,th,bins,vth,vlen,llim,verbose,Large,future)
{
  
  ## Calculate variogram using the largest bins
  ii <- dnb(q0,obs,th,future=future) 
  if( length(ii)<=5 ){
    if(verbose) cat(q0,'k= ',length(ii),'\n')
    r <- (list(fout=NULL,nbr=NULL,ret=3))
  }
  else{
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
        r <- list(fout=fout,nbr=nbr,ret=0)
      }
      else{
        r <- (list(fout=NULL,nbr=nbr,ret=4))
      }
    }else if(ssout$nt <= llim[2]){
      if (verbose) cat('insufficient time points.\n')
      r <- (list(fout=NULL,nbr=nbr,ret=1))
    }else if(ssout$ns <= llim[1]){
      if (verbose) cat('insufficient space points.\n')    
      r <- (list(fout=NULL,nbr=nbr,ret=2))
    }else{
      if (verbose) cat('insufficient space & time points.\n')
      r <- (list(fout=NULL,nbr=nbr,ret=3))
    }    
  }
  
  ## prepare output
  output <- rep(NA,nrow(bins))
  ## working function for Kriging
  working.cltsk.calkriging <- function(q0, nbr, binsth,fout,verbose,future)
  {
    jj <- dnb(q0,nbr,binsth,future=future)
    if(verbose)
    {
      cat(binsth,'k=',length(jj),'\n')
    }
    output <- NA
    if( length(jj) >=3)
    {
      gout <- cal.gamma(q0,nbr[jj,],fout)
      output <- with(gout,work.kriging(Gamma,gamma,dat[,4])[1])
    }
    output
  }
  if( !is.null(r$fout))
  {
    output <- apply(bins,1, working.cltsk.calkriging,q0=q0,nbr=r$nbr,fout=r$fout,verbose=verbose,future=future)
  }
  ## exit here
  output
}
