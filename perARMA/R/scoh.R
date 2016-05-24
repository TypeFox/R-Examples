 scoh <-function(x,m,win,...) 
  {

  scoh_full <-
  function(x,m,win,pfa,plflg,bfflg,ix,iy,nx,ny,datastr) 
  {

      if (length(win)!=m)
      {stop(" 'win' not of length 'm' ")}
 

      if (ix<0 || iy < 0)
      {stop("'ix' or 'iy' is negative")}

       nx0=length(x) 
       if (nx>nx0 || ny>nx0)
       {stop("'nx' or 'ny' > length(x)")}


      xt=fft(x-mean(x))
      mby2=floor(m/2)
      xt<-c(xt[1:mby2],xt,xt[(nx0-mby2+1):nx0])

      nxt=length(xt)

      ddiag=xt*Conj(xt)
      sdiag=signal::conv(win,ddiag)         
      sdiag=sdiag[m:length(sdiag)]   

      if (ix==iy & nx==ny)      
      { coh=diag(nx)     
        for (i in 1:(nx-1))
         {for (j in (i+1):ny)
          {
            indx=seq(0,(m-1))+ix+i                            
            indy=seq(0,(m-1))+iy+j-1                         
            corr=base::sum((xt[indx]*Conj(xt[indy]))*win)
            coh[i,j]=corr%*%Conj(corr)/(sdiag[ix+i]%*%sdiag[iy+j-1])
            coh[j,i]=coh[i,j]}
         }
     }
     if ( bfflg ) {                         
            pfa_corr=pfa*floor(m/2)/nx
            thrs = 1 - exp((log(pfa_corr)/(m-1)))
            } else {
            thrs = 1 - exp((log(pfa)/(m-1)))
           }
       ind<-which(Re(coh[,])<thrs)
       coh[ind]=0



     if (plflg) {
   frequency_p <- seq(ix,nx)
   frequency_q<- seq(iy,ny)
   image(frequency_p, frequency_q, Re(coh), col=gray((64:0)/64), axes = FALSE)
   axis(1, at=seq(ix,nx, by=10))
   axis(2, at=seq(ix,ny, by=10))
   box()
   title(main="Squared Coherence Statistic ",sub=paste("dataname=", datastr," for m= ",m))
  }
}
L<-modifyList(list(ix=0,iy=0,nx=length(x)/2,ny=length(x)/2,pfa=1,plflg=1,bfflg=1,datastr='data'), list(x=x,m=m,win=win,...))

 do.call(scoh_full,L)

 }

