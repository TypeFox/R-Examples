`horseshoe` <-
function(x,y, r1= 1, r2= 1.2, h1= .5, h2= .5, rot=list(cs=1, sn=0), col='black',  lwd=lwd, fill=FALSE  )
  {
    if(missing(col)) { col='black' }
    if(missing(lwd)) { lwd=1 }
     if(missing(r1)) { r1 = 1 }
     if(missing(r2)) { r2 = 1.2 }
     if(missing(h1)) { h1 = .5 }
     if(missing(h2)) { h2 = .5 }
     if(missing(rot)) { rot=list(cs=1, sn=0) }
     if(missing(fill)) { fill=FALSE }

    
 pin = par("pin")
  u = par('usr')
##  h = 1 mm
##  r =0.875 mm
## uinch =   (u[4]-u[3])/pin[2]
###  2.54 cn per inch
  umm =   (u[4]-u[3])/pin[2]/25.4
###  2.54 cn per inch
###25.4 mm/inch
  radx = r1*umm
  rady = r2*umm

  leg1  = h1*umm
  leg2  = h2*umm

    
       
    ########  r1= 1; r2= 1.2; leg1= .5; leg2= .5
    n = 24
    iang  = pi * seq(from = 0, to = 180, length = n)/180
    
    cx = radx*cos(iang)
    cy = rady*sin(iang)

    px=c(cx[1], cx, cx[n])
    py=c(-leg1, cy, -leg2)

    cs = rot$cs
    sn = rot$sn

    HORSE  = list()
      for(i in 1:length(x))
      {
        
        ax = x[i]+(px*cs[i]-py*sn[i])
        ay = y[i]+(px*sn[i]+py*cs[i])
        lines(ax,ay, col=col, lwd=lwd )
        if(fill) {
          polygon(ax,ay, col=col, border=col)
        }
        else
          {
            lines(ax,ay, col=col, lwd=lwd )

          }
        HORSE[[i]] = list(x=ax, y=ay)
      }

    invisible(HORSE)
  }

