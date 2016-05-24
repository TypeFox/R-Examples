`complex.hodo` <-
function(nbaz, dt=dt, labs=c("Vertical", "North", "East"), COL=rainbow(100) , STAMP="")
  {
    if(missing(labs))  { labs=c("Vertical", "North", "East") }
    if(missing(COL))  {  COL=rainbow(100)  }
    if(missing(dt))  {  dt=1  }
    if(missing(STAMP))  {  STAMP=""  }

    
######    par(mfrow=c(1,3))
######
######


    V = nbaz[,1]
    N = nbaz[,2]
    E = nbaz[,3]
    xx = range(E, na.rm =TRUE)
    yy = range(N, na.rm =TRUE)
    zz = range(V, na.rm =TRUE)


    
    sx = range(c(xx, yy))



    n = length(N)
    ncol = length(COL)

    if( ncol >3 )
      {
        
        KR = n/(ncol-1)
        
        KX = floor(seq(from=1, length=n)/(KR))+1
        
        cols=COL[KX]
      }
    else
      {

        cols = 1

      }
    
    ####

    par(yaxs='i', xaxs='i')
    
    plot(c(0, 3), c(0,2), type='n', asp=TRUE, axes=FALSE, ann=FALSE, yaxs='i', xaxs='i', xlab="", ylab="")
   #### title(sub=STAMP)

   ver  = RPMG::RESCALE(V, 1.66, 2, zz[1], zz[2])
   nor  = RPMG::RESCALE(N, 1.33, 1.66, yy[1], yy[2])
    eas  = RPMG::RESCALE(E, 1, 1.33, xx[1], xx[2])
    timvec = seq(from=0, length=length(V), by=dt)

    pax = pretty(timvec, n=10)
    
    xt = RPMG::RESCALE(timvec, 0, 3, min(timvec), max(timvec) )
    axt = RPMG::RESCALE(pax, 0, 3, min(timvec), max(timvec) )

   ####  axis(side=1, pos=2, at=axt, labels=pax)
    axis(side=3, tck=0.01, at=axt, labels=pax)


    
    segments(xt[1:(n-1)], ver[1:(n-1)],   xt[2:(n)], ver[2:(n)] , col=cols)
    segments(xt[1:(n-1)], nor[1:(n-1)],   xt[2:(n)], nor[2:(n)],  col=cols)
    segments(xt[1:(n-1)], eas[1:(n-1)],   xt[2:(n)], eas[2:(n)],  col=cols)
    
text(0, 1.66+0.5*.33, labels=labs[1], pos=4)
text(0, 1.33+0.5*.33, labels=labs[2], pos=4)
text(0, 1+0.5*.33, labels=labs[3], pos=4)

text(0, 2-0.1*.33, labels=STAMP, pos=4)

    
    segments(0, 0.5, 3, 0.5 , lty=2, col=rgb(.8, .8, .8))
    segments(c(.5, 1.5, 2.5), rep(0, 3), c(.5, 1.5, 2.5), rep(1, 3)  , lty=2, col=rgb(.8, .8, .8))
    
 ##  plot(c(0, 3), c(0,1), type='p', asp=TRUE, axes=FALSE, ann=FALSE, xlab="", ylab="")


    tics = pretty(sx)

    atics = RPMG::RESCALE(tics, 0, 1, sx[1], sx[2])
    w = which(atics>0.0 & atics<1)

    tics = tics[w]
    atics = atics[w]

    rect(0, 0, 1, 1)

   x  = RPMG::RESCALE(E, 0, 1, sx[1], sx[2])
   y  = RPMG::RESCALE(N, 0, 1, sx[1], sx[2])

   ## lines(x,y)
    points(x[1], y[1], pch=6) 
segments(x[1:(n-1)],y[1:(n-1)],x[2:n],y[2:n], col=cols)
DO.PMOT.ARR(x,y)

    
   #########   axis(2, labels=tics, at=atics)
    
    gaddtix(side=2, pos=0,   tck=-0.01, at=atics, labels=tics, col=1)
    gaddtix(side=1, pos=0,   tck=-0.01, at=atics, labels=tics, col=1)
    
    text(0, 0.8, labels=labs[2], srt=90, pos=4, offset=.75)
    
    text(0.8, 0, labels=labs[3],  pos=3)

    
      rect(1, 0, 2, 1)

   x  = RPMG::RESCALE(E, 1, 2, sx[1], sx[2])
   y  = RPMG::RESCALE(V, 0, 1, sx[1], sx[2])

  #####  lines(x,y)
    points(x[1], y[1], pch=6)
    segments(x[1:(n-1)],y[1:(n-1)],x[2:n],y[2:n], col=cols)
DO.PMOT.ARR(x,y)



    atics = RPMG::RESCALE(tics, 1, 2, sx[1], sx[2])

        gaddtix(side=1, pos=0,   tck=-0.01, at=atics, labels=tics, col=1)
    text(1, 0.8, labels=labs[1], srt=90,  pos=4, offset=.75)
    text(1.8, 0, labels=labs[3],  pos=3)

     rect(2, 0, 3, 1)
  x  = RPMG::RESCALE(N, 2, 3, sx[1], sx[2])
   y  = RPMG::RESCALE(V, 0, 1, sx[1], sx[2])
     points(x[1], y[1], pch=6)
segments(x[1:(n-1)],y[1:(n-1)],x[2:n],y[2:n], col=cols)
DO.PMOT.ARR(x,y)

   atics = RPMG::RESCALE(tics, 2, 3, sx[1], sx[2])

        gaddtix(side=1, pos=0,   tck=-0.01, at=atics, labels=tics, col=1)
      text(2, 0.8, labels=labs[1], srt=90, pos=4, offset=.75)
    text(2.8, 0, labels=labs[2],  pos=3)

    invisible(sx)
  }

