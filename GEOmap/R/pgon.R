`pgon` <-
function(x,  y, siz=siz, col="black", border=NULL, K=5, startalph = -45, ... )
  {
    if(missing(col)) { col=rep(1,length=length(x)) }
    if(missing(siz)) { siz=rep(.3,length=length(x)) }
    if(missing(border)) { border= rep(1,length=length(x))}
    if(missing(K)) { K=5 }

    if(K<3) { K = 3 }

    if(length(siz)==1) {siz=rep(siz,length=length(x))  }
    if(length(col)==1) {col=rep(col,length=length(x))  }
    if(length(border)==1) {border=rep(border,length=length(x))  }


    
    if(missing(startalph))   { startalph = -45}

    up = par("usr")
    ui = par("pin")

    ratx = (up[2]-up[1])/ui[1]
    raty=  (up[4]-up[3])/ui[2]

    p.x = vector()
    p.y = vector()

    alph = (360/K)*pi/180
    stalph = startalph*pi/180

    angs = stalph+seq(from=0, to=2*pi, by=alph)
    phi = stalph+seq(from=0, to=2*pi, by=alph)
    
    for(i in 1:length(x))
      {
        x0 = x[i]
        y0 = y[i]

        usizx = siz[i]*ratx
        usizy = siz[i]*raty

            p.x =  x0+(usizx*cos(phi)+usizx*sin(phi))
            p.y =  y0+(-usizy*sin(phi)+usizy*cos(phi))

### text(p.x, p.y, labels=1:5)
        if(!is.null(col[i])) {  polygon(p.x, p.y, col=col[i], ...) }
        if(!is.null(border)) { lines(p.x, p.y, col=border[i], ...) }

      }

  }
