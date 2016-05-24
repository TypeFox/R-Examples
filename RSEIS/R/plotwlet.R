`plotwlet` <-
function(baha, Ysig, dt , zscale=1,  zbound=NULL, col=rainbow(100) , ygrid=FALSE, STAMP="", xlab="Time, s" , units="", scaleloc=c(0.4,0.95) )
  {
    
   

    if(missing(zscale)) { zscale = 1 }
    
    if(missing(col)) { col=rainbow(50) }
    
     if(missing(ygrid)) { ygrid=FALSE }
       if(missing(STAMP)) { STAMP=NULL }
       if(missing(zbound)) { zbound=NULL }


    perc = 0.85

    
    if(baha$flip==TRUE)
      {
        
###NO: yax = rev(baha$nvoice*((1:baha$noctave)-1)/(baha$nvoice*baha$noctave))
        yax = rev(baha$nvoice*((1:baha$noctave))/(baha$nvoice*baha$noctave))
        
        
      }
    else
      {
        
        yax =baha$nvoice*((1:baha$noctave)-1)/(baha$nvoice*baha$noctave)
      }
    
    

    a = Ysig

    ##   DSPEC = baha$img
    numfreqs = nrow(baha$img)
    y =  yax
    x = dt*(1:nrow(baha$img))
     
    yflag = rep(TRUE, length(y))

    tim = dt*seq(0, length=length(a))
   
    
  ##   image(x,why,t(DSPEC[1:(numfreqs/2),]), add=TRUE, col = col,xlab='time', ylab='freq', axes=FALSE)

   ##   image(x,why,log10(t(DSPEC[1:(numfreqs/2),])), add=TRUE, col = col,xlab='time', ylab='freq', axes=FALSE)

    IMAT =baha$img

 ##print(y)
    why   = sort( RPMG::RESCALE( 1:ncol(baha$img) , 0 , perc , 1, ncol(baha$img) ))
 ##print(why)
    

    if(zscale<=1){ ImPlot = IMAT; Iunits="Amp" } 
    if(zscale==2){ ImPlot =     log10(IMAT) ; Iunits="Log Amp"}
    if(zscale==3){ ImPlot = sqrt(IMAT) ; Iunits="SQRT Amp"}
    if(zscale==4){ ImPlot = 20*log10( IMAT/ max( IMAT, na.rm=TRUE) ) ; Iunits="DB"}
       if(zscale>4){ ImPlot = IMAT; Iunits="Amp" }  

   ##   par(mfrow=c(1,1))

  ##   old.par <- par(no.readonly = TRUE)
  ##   on.exit(par(old.par))

    MAI = par("mai")
    MAI[4] = MAI[2]
    par(mai=MAI)

    
    par(xaxs='i', yaxs='i')

    plot(range(tim), c(0,1), axes=FALSE, type='n', xlab='', ylab='log(scale)')

    if(!is.null(zbound))
      {
        image.default(x=x, y=why , z=ImPlot , add=TRUE, col = col, zlim=zbound , xlab='time', ylab='log(scale)', axes=FALSE)
  }
    else
      {
        image.default(x=x, y=why , z=ImPlot , add=TRUE, col = col,  xlab='time', ylab='log(scale)', axes=FALSE)

      }

    if(!is.null(baha$ridge))
      {

         image(x=x, y=why , z=baha$ridge , add=TRUE, col = tomo.colors(100))


      }
    
    trace = RPMG::RESCALE( a, perc , 1.0  , min(a), max(a) )

    
   lines(tim, trace)

    ##  sy = RPMG::RESCALE( a, perc  , 1.0  , min(a), max(a) )
    Tdiff = max(tim)-min(tim)
    
    ##   segments(max(tim)-Tdiff*.04-DEVOL$wpars$Ns*dt, 0.76, max(tim)-Tdiff*.04, 0.76, lwd=2)
    

                                        # axis(1)
                                        #  axis(3)
    
    xtix = pretty(x, n=10)

    xedge = (x[2]-x[1])
    xleft = min(x)-xedge
    xright = max(x)+xedge
     xtix = xtix[xtix>=xleft & xtix<=xright]

    ### print(xtix)
    
 #    xtix = c(floor(min(x)),xtix,  floor(max(x)))
    
    axis(3,tck=.01,at=xtix,labels=FALSE)
    mtext( side=3,    at=xtix, text=xtix, line=.5)


    axis(1,tck=.01,at=xtix,labels=FALSE)
    mtext( side=1,    at=xtix, text=xtix, line=.25)

   
                                        #  title(xlab="Time, s")
    mtext(side=1, at=max(x), text=xlab , line=1.5, adj=1)


    if(!is.null(STAMP))
      {
        
        ###   mtext(side=1, at=max(x), text=STAMP , line=2.5, adj=1)
        mtext(side=3, at=0, text=STAMP , line=1.5, adj=0)
 
      }
    
      

    raxspec= RPMG::RESCALE(yax , 0 , perc , 0, 1 )
      
    axis(2, at=raxspec, labels=2^(1:baha$noctave))
    
   ##    axis(2, at=raxspec, labels=format.default(axspec, digits=3), pos=min(x))

    if(ygrid==TRUE)
      {
      
        ##  abline(h=raxspec)
        segments(rep(min(x), length(raxspec)), raxspec,  rep(max(x), length(raxspec)) , raxspec)

      }

    ###  right axis plotted on the time series
    ##   need a new version here that looks like the matlab axis


    
    axtrace = range(a)
    ## print(axtrace)
    
    raxtrace= RPMG::RESCALE( axtrace, perc , 1.0 , min(a), max(a) )

    Elabs =  RPMG::endSCALE(axtrace)

    exp = parse(text = Elabs[3])
 ##   mtext(Elabs[1], side = 4,at = perc , line=0.5, cex=0.8,las=2)
 ##   mtext(Elabs[2], side = 4,at = 1, line=0.5, , cex=0.8,las=2)
 
    
    axis(4, at=raxtrace , labels=Elabs[1:2] , pos=max(tim), tick=TRUE , line=0.5, cex.axis=0.8,las=2)

       mtext(exp, side = 3, at = max(tim), line=0.5, adj=-1  , cex=0.8)

    mtext(units, side = 4, at =mean(raxtrace) , line=0.5   , cex=0.8 ,las=1 )


    
   # axis(4, at=mean(c(perc, 1)) , labels=wiglab , pos=max(tim), tick=FALSE , line=-0.5, cex.axis=0.8)

####


    
    if(!is.null(zbound))
      {
        RPMG::HOZscale(zbound, col, units=Iunits, s1=scaleloc[1], s2=scaleloc[2])
      }
    else
      {
        RPMG::HOZscale( ImPlot, col, units=Iunits, s1=scaleloc[1], s2=scaleloc[2])
        
      }
    invisible(list(y=y[yflag], why=why, yBounds=c(0,perc), x=x, yat=raxspec))
    
  }

