`plt.MTM0` <-
  function(frange,prange, plxy, M, freqs, amp , a,  dof=dof, Fv=Fv, COL=2)
  {
##########  plot an MTM spectrum
    
###  plt.MTM0(frange,prange, plxy, M, freqs, amp )
    if(missing(dof)) { dof=NULL }
    if(missing(Fv)) { Fv=NULL }
    if(missing(COL)) { COL=1:M }

    paxis = unlist(strsplit(plxy, split=""))
    
    if(any(paxis=="D"))
      {
        pmax = max(prange)
        prange  = 10*log10(prange/pmax )
        
        if(any(paxis=="x")) { plxy="x" } else { plxy="" }
        plot(frange,prange,type='n',log=plxy, axes=FALSE, xlab="Hz", ylab="Db")
        for(i in 1:M)
          {
            amp[[i]]  = 20*log10(amp[[i]]/pmax )
          }
      }
    else
      {
        plot(frange,prange,type='n',log=plxy, axes=FALSE, xlab="Hz", ylab="Spec")
      }
    
    if(!is.null(dof))
      {
        
        for(i in 1:M)
          {
            why   = RPMG::RESCALE( dof[[i]], prange[1]  ,prange[2] , min(dof[[i]], na.rm = TRUE)  , max(dof[[i]], na.rm = TRUE)  )
            lines(freqs[[i]], why, col=COL[i], lty=2)
          }  
        
        
      }
    
    if(!is.null(Fv))
      {
        ppoints  =  c( 97.5)
        q1 = qf(ppoints/100, 2, 8)
        
        
        for(i in 1:M)
          {
            
            w = which(Fv[[i]]>q1)
            points(freqs[[i]][w], amp[[i]][w], col="brown", pch=8)
########  why   = RPMG::RESCALE( Fv[[i]], prange[1]  ,  prange[1] + (prange[2]-prange[1])/2 , min(Fv[[i]], na.rm = TRUE)  , max(Fv[[i]], na.rm = TRUE)  )
########  lines(freqs, why, col=i, lty=3)
          }  
        
        
      }
    
    
    for(i in 1:M)
      {
        lines(freqs[[i]], amp[[i]], col=COL[i], lty=1)
        text(freqs[[i]][1], amp[[i]][1], i, col=COL[i], pos=2)
      }
    if(!is.null(a$stamps))
      {
        legend("topright", lty=1, col=COL[1:M], legend=a$stamps)
      }
    
    axis(2, las=2)
    axis(1)
    box()
    
    
    
  }

