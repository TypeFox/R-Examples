`brune.doom` <-
function(amp, dt=1, f1=0.01, f2=15,  PLOTB=FALSE ,  tit="")
{
  ######  simpler version of brune.fit?  this is what David used I believe:  brune.doom(amp,dt0,1,20)

  ##  fit a brune model to a windowed trace
  ##  brune.doom( amp, dt, f1=0.01, f2=15, PLOT=FALSE  )
#######   function depends on:next2,  MTMdisp ,get.corner ,  brune.search
    ######   returns list:  WARN, tstar0 , gamma,omega0,  fc, alpha


  if(missing(dt)) { dt = 1 }
  if(missing(f1))
    {  f1 = 0.01 }
  if(missing(f2))
    { f2 = 14.0 }

  if(missing(PLOTB))
    { PLOTB=FALSE  }


  if(missing(tit))
    { tit = ""}

  ###  return this default structure, WARN is a warning for failure
  xc = list(SUCCESS=TRUE, WARN = "OK",tstar0 = 0, gamma = 0,
    omega0 = 0,
    fc = 0,
    alpha= 0)

  #######   X-axis coordinates
  ex = seq(0,length(amp)-1)*dt

  xv  = ex

  ####   remove the mean
  ampv = amp-mean(amp)
  a = list(y=ampv, dt=dt)

  if(length(a$y)<2)
    {
      return(list(WARN=FALSE, corn=0,ave=0,slope=0,interc=0,tstar0=0,omega0=0 ))
    }
  ##  ta = ts(a$y, start=0, deltat=dt)
  len2 = 2*next2(length(ampv))
  if(len2<1024)
    { len2 = 1024 }

#######  calculat the MTM spectrum for displacement (integrate the velocity seismogram)
  
  Spec =MTMdisp(a, f1=f1, f2=f2, len2=len2, PLOT=FALSE )
                                        #   Spec = spec(a, f1=f1, f2=f2, len2=len2, PLOT=FALSE )

  lspec = Spec$displ
  

  # print(paste(sep=' ', "brune.TEST", length(Spec$f), length( lspec), dt, f1, f2))
  
  #####  use grid search to find an initial corner frequency
  xc = get.corner(  Spec$f , lspec, dt, f1, f2, PLOT=FALSE)
  
  ## print(paste(sep=' ', "BF post", xc$omega0, xc$corn, xc$tstar0))

  fnyq = 1/(2*dt)

  #########  checks
  if(xc$corn>fnyq|xc$corn<=0|xc$tstar0<0)
    {
      xc$SUCCESS = FALSE;
      xc$WARN = "xc$corn>fnyq|xc$corn<=0|xc$tstar0<0";

      if(PLOTB==TRUE)
        {
          
          par(mfrow=c(1,1))
          flag= Spec$f>=f1 & Spec$f<=f2
          
          print(paste(sep=' ', "Brune DOOM PLOTB==TRUE", length(Spec$f[flag]), length( lspec[flag]), dt, f1, f2))
          
          LY = log10(lspec[flag])
          LF = log10(Spec$f[flag])
          print(paste(sep=' ', "Brune DOOM", length(LF), length(LY) ))
          plot(LF, LY ,  lty=1 , type='l' , xlab="Log Freq", ylab="Log Amp Spec")
          
          ## jbrune = brune.func(Spec$f[flag], jmod$omega0, jmod$tstar0 , jmod$fc,  jmod$alpha, jmod$gamma )
          ##   lines(LF, log10(jbrune), col=3)
          title(main="xc$corn>fnyq|xc$corn<=0|xc$tstar0<0")
          
        }
      return(xc)

    }

  if(is.numeric(xc$omega0)==TRUE & is.numeric(xc$corn)==TRUE  & is.numeric(xc$tstar0)==TRUE )
    {
      ########   do search for best fitting brune model
      jmod = brune.search(Spec$f, lspec, f1, f2,  xc$omega0, xc$corn, xc$tstar0, 2.0)
    }
  else
    {
      xc$SUCCESS = FALSE;
      xc$WARN = "Non-numeric";
      return(xc)
    }


  if(PLOTB==TRUE)
    {
      
      par(mfrow=c(1,1))
      flag= Spec$f>=f1 & Spec$f<=f2
      
      print(paste(sep=' ', "Brune DOOM PLOTB==TRUE", length(Spec$f[flag]), length( lspec[flag]), dt, f1, f2))

      LY = log10(lspec[flag])
      LF = log10(Spec$f[flag])
      
      print(paste(sep=' ', "Brune DOOM", length(LF), length(LY) ))
      plot(LF, LY ,  lty=1 , type='l' , xlab="Log Freq", ylab="Log Amp Spec")
      jbrune = brune.func(Spec$f[flag], jmod$omega0, jmod$tstar0 , jmod$fc,  jmod$alpha, jmod$gamma )
      lines(LF, log10(jbrune), col=3)


###### add.brune(list(y=Spec$f, dt=dt)  ,  xc$tstar0,   xc$corn,   xc$omega0, 4, f1=f1, f2=f2, plog=TRUE)
      abline(h=xc$ave, col=4)
      abline(xc$interc, xc$slope, col=2)
      
###  format.default(TIM$jd[1], digits=2, trim=FALSE)
      tit2 = paste(sep=' ',
        format.default(jmod$omega0, digits=4, trim=FALSE),
        format.default(jmod$tstar0, digits=4, trim=FALSE),
        format.default(jmod$fc,     digits=3, trim=FALSE),
        format.default(jmod$alpha,  digits=3, trim=FALSE),
        format.default(jmod$gamma,  digits=3, trim=FALSE))
      
      title(main=tit2)
     
     
      
      
    }

  ###  all went well, return the model parameters
  xc$SUCCESS = TRUE;
  xc$WARN = "OK";
  xc$tstar0 = jmod$tstar0;
  xc$gamma = jmod$gamma;
  xc$omega0 = jmod$omega0;
  xc$fc = jmod$fc;
  xc$alpha= jmod$alpha
  xc$chisqrd = jmod$chisqrd
   
  return(xc)
  
}

