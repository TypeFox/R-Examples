`MTMdisp` <-
function(a, f1=f1, f2=f2, len2=1024, PLOT=FALSE)
  {

    ###  calculate and plot an MTM spectrum
#   a = list(y=ampv, dt=0.008)

  if(missing(PLOT)) { PLOT=TRUE }
  if(missing(f1)) { f1 = 0.01 }
  if(missing(f2)) { f2 = 10 }
  len = length(a$y)
  if(missing(len2))
    {
      len2 = 2*next2(len)
    }
  if(len2<len)
    {
      len2 = 2*next2(len)
    }
  Mspec =   mtapspec(a$y,a$dt, klen=len2,  MTP=list(kind=1,nwin=5, npi=3,inorm=0)  )

  f=Mspec$freq

  amp = Mspec$spec[1:length(f)]

  #  sam = lowess(f,amp, f=10/length(f));
  #  sam$y[sam$y<=0] = amp[sam$y<=0];
  
  #  ma = cbind(amp, sam$y);
  ma = amp;
  flag = f>=f1 & f <= f2;
  displ = ma/(2*pi*f);

  if(PLOT==TRUE)
    {
                                        # 	matplot(f[flag],displ[flag,],type='l',log='xy',axes=FALSE, xlab="Hz")
  
      plot(range(f[flag]),range(displ[flag]),type='n',log='xy',axes=FALSE, xlab="Hz", ylab="Disp Spec")
      lines(f[flag], displ[flag], col=1, lty=1)       
      axis(2, las=2)
      axis(1)
      box()
    }
  invisible( list(len2=len2, f=f, f1=f1, f2=f2, displ=displ, ampsp=amp, flag=flag ) )
  }

