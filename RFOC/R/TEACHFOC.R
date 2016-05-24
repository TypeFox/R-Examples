`TEACHFOC` <-
function(s,d,r, up=FALSE)
  {
    ############  given strike dip and rake, plot focal mechanism with
    ############   special vectors shown
    if(missing(up)) up = FALSE
    mc = CONVERTSDR(s,d,r )
    MEC = MRake(mc$M)
    MEC$icol =  foc.icolor(MEC$rake1)
    MEC$ileg =  focleg(MEC$icol)
    MEC$fcol =   foc.color(MEC$icol)
    MEC$UP = up
    Beachfoc(MEC, fcol= MEC$fcol )
    
    text(0,1, labels="N", pos=3, xpd=TRUE)
    text(1,0, labels="E", pos=4, xpd=TRUE)
    text(-1,0, labels="W", pos=2, xpd=TRUE)
    text(0,-1, labels="S", pos=1, xpd=TRUE)

    F = focpoint(MEC$F$az, MEC$F$dip, col=6,  lab="F", UP=MEC$UP)
    G = focpoint(MEC$G$az, MEC$G$dip, col=6,  lab="G", UP=MEC$UP)  
    P = focpoint(MEC$P$az, MEC$P$dip, col=6,  lab="P", UP=MEC$UP)
    T = focpoint(MEC$T$az, MEC$T$dip, col=6,  lab="T", UP=MEC$UP)
    V = focpoint(MEC$V$az, MEC$V$dip, col=6,  lab="V", UP=MEC$UP)
    U = focpoint(MEC$U$az, MEC$U$dip, col=6,  lab="U", UP=MEC$UP)
    segments(c(-.02, 0), c(0, -0.02), c(0.02, 0), c(0, 0.02), col='black' )

    
    PLNS = PlotPlanes(MEC, col1="blue", col2="green")

    naz = 100

    AZIM = min(RPMG::fmod(MEC$az1-180, 360), MEC$az1)
    
    
    alph  = 90-seq(from=0,to=AZIM , length=naz)
    erad = 1.01

    ex  = erad*cos(pi*alph/180)
    why = erad*sin(pi*alph/180)
    lines(ex,why, lwd=2, lty=2)
    arrows(ex[naz-1],why[naz-1], ex[naz],why[naz])
    iaz = floor(naz/2 )
    text(ex[iaz],why[iaz], labels=paste("Az=", formatC(AZIM, digits=6))  , pos=4, font=2, xpd=TRUE)


    NIP = nipXY( MEC,  0, 0, fcol = "blue" , nipcol='black',  size = c(1,1), cex=1 )

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
    
    upperlower = "Lower Hemisphere"
    if(MEC$UP)  upperlower = "Upper Hemisphere"
    
    title(main=paste(sep=" ", "Strike=", formatC(s, digits=6) , "Dip=", formatC(d, digits=6), "Rake=", formatC(r, digits=6)), xlab= upperlower)
    
    
  }

