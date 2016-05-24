`MTMplot` <-
function(a, f1=f1, f2=f2, len2=1024, PLOT=FALSE, PADDLAB=NULL, GUI=TRUE)
  {

    ###  calculate and plot an MTM spectrum
#   a = list(y=ampv, dt=0.008)

  if(missing(PLOT)) { PLOT=TRUE }
  if(missing(GUI)) { GUI=TRUE }
  
  if(missing(f1)) { f1 = 0.01 }
  if(missing(f2)) { f2 = 10 }

  if(missing(PADDLAB)) { PADDLAB=NULL}
  
  
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
  displ = ma ;

  if(PLOT==TRUE)
    {
      stdlab = c("DONE", "X-LOG", "Y-LOG", "VALS")
      labs = c(stdlab, PADDLAB)
      NLABS = length(labs)
      NOLAB = NLABS +1000
      colabs = (1:length(labs))
      pchlabs = rep(0,length(labs))
   
      plot(range(f[flag]),range(displ[flag]),type='n',log='',axes=FALSE, xlab="Hz", ylab="Spec")
      lines(f[flag], displ[flag], col=1, lty=1)       
      axis(2, las=2)
      axis(1)
      box()

      if(GUI==FALSE) { return( list(len2=len2, f=f, f1=f1, f2=f2, displ=displ, ampsp=amp, flag=flag ) ) }
      
      u = par("usr")
      sloc = list(x=c(u[1],u[2]))
      
      
      buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
      
      zloc = zlocator(COL=rgb(1,0.8, 0.8), ID=TRUE, NUM=FALSE , YN=1, style=0)
      Nclick = length(zloc$x)
      if(is.null(zloc$x)) { return(NULL) }
      K = RPMG::whichbutt(zloc ,buttons)

      sloc = zloc
      plogx=''
      plogy=''
      
      
      
      while(Nclick>0)
        {
          
          if(K[Nclick]==1)
            {
              break;
            }
          
          
     if(Nclick==1 & K[Nclick]==0)
       {
         plxy = NULL
       }
          
          
          
          if(K[Nclick]==2)
        {
          if( (plogx=='x')==TRUE ) { plogx = '' }
          else { plogx = "x" }
        }
          if(K[Nclick]==3)
            {
              if( (plogy=='y')==TRUE) { plogy = '' }
          else { plogy = "y" }
              
              
            }
          if(K[Nclick]==4)
            {
          alabs = format.default(zloc$x[1:(Nclick-1)], digits=3)
          print(paste( paste( sep=' ',"Frequencies:", paste(alabs, collapse=" "))))
          ## abline(v=zloc$x[1:(Nclick-1)])
          ## mtext(labs, at=zloc$x[1:(Nclick-1)], side=3, line=0)
        
        }

      plxy = paste(sep='', plogx , plogy)

      plot(range(f[flag]),range(displ[flag]),type='n',log=plxy,axes=FALSE, xlab="Hz", ylab="Spec")
      lines(f[flag], displ[flag], col=1, lty=1)       
      axis(2, las=2)
      axis(1)
      box()
        buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
 
      zloc = zlocator(COL=rgb(1,0.8, 0.8), NUM=FALSE , ID=TRUE ,YN=1, style=0)
      Nclick = length(zloc$x)
      if(is.null(zloc$x)) { return(sloc) }
      K =  RPMG::whichbutt(zloc ,buttons) 
      
    }


   
    }
  invisible( list(len2=len2, f=f, f1=f1, f2=f2, displ=displ, ampsp=amp, flag=flag ) )
  }

