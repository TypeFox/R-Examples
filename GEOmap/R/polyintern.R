polyintern<-function(P, n=10, PLOT=FALSE)
  {

    if(missing(n)) n=10
    if(missing(PLOT) )  PLOT=FALSE
    
    
    x = P$x
    y = P$y

    xo = seq(from=min(x), to=max(x), length=n)
    
    yo = seq(from=min(y), to=max(y), length=n)

    mg = RPMG::meshgrid(xo, yo)

    
    
   ###   I = inpoly(mg$x, mg$y , list(x=c(P$x, P$x[1]), y=c(P$y, P$y[1])))

    kin = splancs::inout(cbind(as.vector(mg$x), as.vector(mg$y) ) ,cbind(x, y), bound=TRUE)
 
    nx = mg$x[kin]
    ny =  mg$y[kin]

    ef = rep(0, length=length(nx))
    
    for(j in 1:length(nx))
      {
        ## points(nx[j],ny[j],col='blue')
        dis = rep(NA, length=(length(P$x)-1))
        pdx = rep(NA, length=(length(P$x)-1))
        pdy =  rep(NA, length=(length(P$x)-1))


        if(FALSE)
          {
            for(k in 1:(length(P$x)-1))
              {
                Pdis = pline(P$x[k], P$y[k], P$x[k+1], P$y[k+1], nx[j], ny[j])
                ##  points(Pdis[4],Pdis[5], pch=3,col='brown')
                dis[k] = Pdis[1]
                pdx[k]=Pdis[4];
                pdy[k] =Pdis[5] 
                
              }
            wminA  = which.min(dis)
            
          }

        ax = P$x
        ay = P$y
        nel = length(ax)
        ex = nx[j]
        ey = ny[j]
        dis = 0;
        
        
        DDOUT = .C("CALL_polydistpoint",PACKAGE = "GEOmap",
          as.double(ax),
          as.double(ay),
          as.integer(nel),
          as.double(ex), as.double(ey), as.double(dis))
        

        dis = DDOUT[[6]]

        ##  print(dis)
        
        ##  dis = sqrt( (nx[j]-P$x)^2 + (ny[j]-P$y)^2 )

        
        ef[j] = dis
        ##ef[j] = min(dis, na.rm=TRUE)

        ##  segments(nx[j],ny[j], pdx[wminA], pdy[wminA], col='red')
      }

    
    zi = which.max(ef)

    if(PLOT)
      {
        plot(x,y, type='n', asp=1)
        polygon(x,y, col=rgb(.7,.7,1))
        
        points(mg$x[kin] , mg$y[kin] )

        points(nx[zi], ny[zi], col='red', pch=6)
      }
    
    invisible(list(x=nx[zi], y=ny[zi] , zi=zi, nx=nx, ny=ny, ef=ef ))
  }
