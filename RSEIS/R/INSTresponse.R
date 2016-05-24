`INSTresponse` <-
function(Kal,key,ff,tt=tt,plotkey=NULL)
  {
    
    if(missing(plotkey)) { plotkey=NULL }
    
    if(missing(tt)) {   tt=c(20,0.008) }

    verb = FALSE

    nams=names(Kal)
    nam =  nams[key]
    Calib = Kal[[key]]
    norm  =Calib$Knorm;
    gain = Calib$Gain;
    meands=Calib$Sense;
    npole =Calib$np;
    nzero =Calib$nz;
    poles =Calib$poles;
    zeroes=Calib$zeros;

    if(verb)
      {
        print(paste(sep=" ","RESPONSE:", nam, npole, nzero, norm, gain))
        print('poles')
        print(poles)
        print('zeros')
        print(zeroes)
      }

    if(is.null(zeroes))
      {
       bb = 1
      }
    else
      {
    bb    =gpoly(zeroes); ##  convert zeros to polynomial coefficients
  }
    aa    =gpoly(poles);  ##  convert poles to polynomial coefficients


    if(length(ff)==2)
      {
        f=logspace(ff[1],ff[2]);
      }
    else
      {
        f = ff
      }
    w=2*pi*f;
    transfer=INSTFREQS(bb,aa,w)*norm;


if(length(tt)>0)
  {
  sintr=tt[2]                   ## % sample interval
  n=floor(tt[1]/sintr)          ## % number of time points
  n=n+(n%%2)                ##  % force n to be even
  n2=n/2
  tstart=sintr*(1-n2)           ## % time of first sample
  d=rep(0, length=n)                  ## % d is a delta function at time=0
  d[n2]=1/sintr
  
  DRET=FRWDft(d,n,tstart,sintr);      ##% fourier transform of delta function

  f1 = DRET$f
  
  RESP =(INVRft(DRET$G*INSTFREQS(bb,aa,2*pi*f1)*norm,n,tstart,sintr)); ##% inv trans

  resp=Re(RESP$g)
  
}


    

    if(length(plotkey)>0)
      {
        par(mfrow=c(length(plotkey) ,1))

        plotkey = tolower(plotkey)



        m1 = match(c("amp", "mag") , plotkey)
        m1 = m1[!is.na (m1)]
        
        if(  length(m1)>0)
          {
            
            mag=Mod(transfer)
                                        # plot(ff, mag, log='xy', axes=FALSE)

            fx = f>0
            
            ## plot(f[fx], mag[fx], type='l', log='xy', axes=FALSE,   xlab="", ylab="")

            plot(f[fx], 10*log10(mag[fx]/max(mag[fx])), type='l', log='x', axes=FALSE,   xlab="", ylab="")

            abline(h=(-3), lty=2, col=grey(.8))

            
            axis(1)
            axis(2)
            title(main=paste(sep=" ", nam, 'Amplitude transfer function'), xlab="frequency (Hz)", ylab='amplitude' )
            box()

          }

        m2 = match(c("phase") , plotkey)
        m2 = m2[!is.na (m2)]
        

        if(  length(m2)>0)
          {
            
            phase=Arg(transfer);
            plot(f[fx], LocalUnwrap(phase[fx]), type='l', log='x', axes=FALSE,   xlab="", ylab="")
            axis(1)
            axis(2)
            title(main=paste(sep=" ", nam, 'Phase transfer function'), xlab="frequency (Hz)", ylab='phase' )
            box()
          }


        m3 = match(c("resp") , plotkey)
        m3 = m3[!is.na (m3)]
        
        
        if(  length(m3)>0)
          {
            
            ##  locator(1)
            ##   par(mfrow=c(1,1))
            plot(RESP$t, resp, type='l', main=paste(sep=" ", nam,'impulse response function'), xlab="time (s)", ylab='amplitude')
          }
      }
    
    
    
    return(list(transfer=transfer  ,aa=aa, bb=bb, resp=resp))
    
    ##transfer=[]; resp=[]; t=[]; f=[];
  }

