`contwlet` <-
function(baha, Ysig, dt , clev=0.75,  NLEV=12, zscale=1,  zbound=NULL, col=col, ygrid=FALSE,  WUNITS="Volts",  PEAX=NULL)
  {
    if(missing(zscale)) { zscale = 1 }
    if(missing(clev)) { clev=.75 }
    if(missing(NLEV)) { NLEV=12 }
    
    if(missing(col)) { col=rainbow(50) }
    
     if(missing(ygrid)) { ygrid=FALSE }
    
    if(missing(zbound)) { zbound=NULL }
    if(missing(PEAX))  { PEAX=NULL }
         if(missing(WUNITS)) {   WUNITS="Volts"  }
     

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

    DSPEC = baha$img
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
    

    if(zscale<=1){ ImPlot = IMAT; units="Amp" } 
    if(zscale==2){ ImPlot =     log10(IMAT) ; units="Log Amp"}
    if(zscale==3){ ImPlot = sqrt(IMAT) ; units="SQRT Amp"}
    if(zscale==4){ ImPlot = 20*log10( IMAT/ max( IMAT, na.rm=TRUE) ) ; units="DB"}
       if(zscale>4){ ImPlot = IMAT; units="Amp" }  

   ##   par(mfrow=c(1,1))
    par(xaxs='i', yaxs='i')

    plot(range(tim), c(0,1), axes=FALSE, type='n', xlab='', ylab='log(scale)')
    print("PEAX:")
    print(PEAX)
    if(!is.null(PEAX$x))
      {
        abline(v=PEAX$x, col=grey(.6))
      }

    
    ###   just contour the high values
    ranz = range(ImPlot, na.rm=TRUE)

    if(!is.null(zbound))
      {
        print("in conwlet: zbound")
        print(zbound)
 
        if(zscale<=1){ ranz = zbound } 
        if(zscale==2){ ranz =     log10(zbound)}
        if(zscale==3){ ranz = sqrt(zbound) }
        if(zscale==4){ ranz = 20*log10(zbound / max(zbound , na.rm=TRUE) )}
        if(zscale>4){ ranz = zbound }  
        

        print("in conwlet: ranz")
        print(ranz)
        
        nlevs = pretty(seq(from=ranz[1], to=ranz[2] , length=NLEV))
      }
    else
      {

    ####  nlevs = pretty(seq(from=ranz[1]+clev*(diff(ranz)), to=ranz[2] , length=12))

    if(length(clev)==2)
      {
        nlevs = pretty(seq(from=ranz[1]+clev[1]*(diff(ranz)), to=ranz[1]+clev[2]*(diff(ranz)) , length=NLEV))
      }
    else
      {
        nlevs = pretty(seq(from=ranz[1]+clev*(diff(ranz)), to=ranz[2] , length=NLEV))
      }

  }

    print("nlevs")
    print(nlevs)
    contour(x=x, y=why , z=ImPlot , add=TRUE, levels=nlevs, xlab='time', ylab='log(scale)', axes=FALSE, drawlabels=FALSE)
    
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
    xtix = xtix[xtix<=max(x)]

 #    xtix = c(floor(min(x)),xtix,  floor(max(x)))
    
    axis(3,tck=.01,at=xtix,labels=FALSE)
    mtext( side=3,    at=xtix, text=xtix, line=.5)


    axis(1,tck=.01,at=xtix,labels=FALSE)
    mtext( side=1,    at=xtix, text=xtix, line=.25)

   
                                        #  title(xlab="Time, s")
    mtext(side=1, at=max(x), text="Time, s" , line=1.5, adj=1)



      

    raxspec= RPMG::RESCALE(yax , 0 , perc , 0, 1 )
      
    axis(2, at=raxspec, labels=2^(1:baha$noctave))
    
   ##    axis(2, at=raxspec, labels=format.default(axspec, digits=3), pos=min(x))

    if(ygrid==TRUE)
      {
      
        ##  abline(h=raxspec)
        segments(rep(min(x), length(raxspec)), raxspec,  rep(max(x), length(raxspec)) , raxspec)

      }


  

    
    axtrace = range(a)
    raxtrace= RPMG::RESCALE( axtrace, perc , 1.0 , min(a), max(a) )


    if( !is.na(WUNITS)  )
          {
            axis(4, at=raxtrace, labels=format.default(axtrace, digits=3), pos=max(tim, na.rm=TRUE))
            mtext(side=4, at=mean(raxtrace), text=WUNITS)
          }
    
    
   ###  RPMG::HOZscale( ImPlot, col, units=units, s1=0.4, s2=0.95)

    invisible(list(y=y[yflag], why=why, yBounds=c(0,perc), x=x, yat=raxspec))

  }

