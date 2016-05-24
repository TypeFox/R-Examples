`plotMEC` <-
function(x, detail=0, up=FALSE)
  {
    ############  given strike dip and rake, plot focal mechanism with
    ############   special vectors shown
    if(missing(up)) up = FALSE
    if(missing(detail)) detail=0
   
    xmec = MRake(x$M)
    xmec$icol =  foc.icolor(xmec$rake1)
    xmec$ileg =  focleg(xmec$icol)
    xmec$fcol =   foc.color(xmec$icol)
    xmec$UP = up
    Beachfoc(xmec, fcol= xmec$fcol )

 if(detail>0)
      {

 if(detail>=1)
      {

        
    
    text(0,1, labels="N", pos=3, xpd=TRUE)
    text(1,0, labels="E", pos=4, xpd=TRUE)
    text(-1,0, labels="W", pos=2, xpd=TRUE)
    text(0,-1, labels="S", pos=1, xpd=TRUE)

    
    F = focpoint(xmec$F$az, xmec$F$dip, col=6,  lab="F", UP=xmec$UP)
    G = focpoint(xmec$G$az, xmec$G$dip, col=6,  lab="G", UP=xmec$UP)  
    P = focpoint(xmec$P$az, xmec$P$dip, col=6,  lab="P", UP=xmec$UP)
    T = focpoint(xmec$T$az, xmec$T$dip, col=6,  lab="T", UP=xmec$UP)
    V = focpoint(xmec$V$az, xmec$V$dip, col=6,  lab="V", UP=xmec$UP)
    U = focpoint(xmec$U$az, xmec$U$dip, col=6,  lab="U", UP=xmec$UP)
    segments(c(-.02, 0), c(0, -0.02), c(0.02, 0), c(0, 0.02), col='black' )
    PLNS = PlotPlanes(xmec, col1="blue", col2="green")
  }
   



  if(detail>=2)
      {

    naz = 100

    AZIM = min(RPMG::fmod(xmec$az1-180, 360), xmec$az1)
    
    
    alph  = 90-seq(from=0,to=AZIM , length=naz)
    erad = 1.01

    ex  = erad*cos(pi*alph/180)
    why = erad*sin(pi*alph/180)
    lines(ex,why, lwd=2, lty=2)
    arrows(ex[naz-1],why[naz-1], ex[naz],why[naz])
    iaz = floor(naz/2 )
    text(ex[iaz],why[iaz], labels=paste("Az=", AZIM)  , pos=4, font=2, xpd=TRUE)


    NIP = nipXY( xmec,  0, 0, fcol = "blue" , nipcol='black',  size = c(1,1), cex=1 )

    thick = 0.01; headlength = 0.2; headthick = 0.1
    fancyarrows(0,0, U$x, U$y, thick =thick , headlength =  headlength, headthick =headthick)
    fancyarrows(0,0, V$x, V$y, thick =thick , headlength =  headlength, headthick =headthick)
    fancyarrows(0,0, NIP$Q$x, NIP$Q$y, thick =thick , headlength =  headlength, headthick =headthick)

    
    rx = sqrt(F$x^2+F$y^2)
    ax = F$x/rx
    ay = F$y/rx

    arrows(ax, ay, rx*ax, rx*ay, lty=2, lwd=2, col='blue')
    text(mean(c(ax, rx*ax)), mean(c(ay,rx*ay) ), font=2, labels="Dip", pos=4)
    

    wm = which.min( (PLNS$LP2$x-NIP$Q$x)^2+(PLNS$LP2$y-NIP$Q$y)^2 )
   
    lines(PLNS$LP2$x[1:wm], PLNS$LP2$y[1:wm], col='purple', lwd=2)
    arrows(PLNS$LP2$x[wm-1], PLNS$LP2$y[wm-1], PLNS$LP2$x[wm], PLNS$LP2$y[wm],  col='purple', lwd=2)
    text((PLNS$LP2$x[round(wm/2)]), (PLNS$LP2$y[round(wm/2)]), labels="Rake", font=2, pos=4)
  }
  }

    
    upperlower = "Lower Hemisphere"
    if(xmec$UP)  upperlower = "Upper Hemisphere"
    
    title(main=paste(sep=" ", "Strike=", x$strike, "Dip=", x$dip, "Rake=", x$rake), xlab= upperlower)
    
    
  }

