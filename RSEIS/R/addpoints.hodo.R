`addpoints.hodo` <-
function(nbaz, dt, sx, flag=1:10, pch=3, col=1)
  {
    if(missing(flag)) { flag = 1:length(nbaz[,2]) }
    if(missing(pch)) { pch = 3 }
    if(missing(col)) { col = "brown" }

    V = nbaz[,1]
    N = nbaz[,2]
    E = nbaz[,3]
    xx = range(E, na.rm =TRUE)
    yy = range(N, na.rm =TRUE)
    zz = range(V, na.rm =TRUE)

    timvec = seq(from=0, length=length(V), by=dt)
    xt = RPMG::RESCALE(timvec, 0, 3, min(timvec), max(timvec) )
 
    V = nbaz[flag,1]
    N = nbaz[flag,2]
    E = nbaz[flag,3]
    xt  = xt[flag]
    
    ver  = RPMG::RESCALE(V, 1.66, 2, zz[1], zz[2])
    nor  = RPMG::RESCALE(N, 1.33, 1.66, yy[1], yy[2])
    eas  = RPMG::RESCALE(E, 1, 1.33, xx[1], xx[2])
    
    points(xt, ver, pch=pch, col=col, cex=.5) 
    points(xt, nor, pch=pch, col=col, cex=.5) 
    points(xt, eas, pch=pch, col=col, cex=.5) 
    
    
    x  = RPMG::RESCALE(E, 0, 1, sx[1], sx[2])
    y  = RPMG::RESCALE(N, 0, 1, sx[1], sx[2])

    ## lines(x,y)
    points(x, y, pch=pch, col=col, cex=.5) 

    x  = RPMG::RESCALE(E, 1, 2, sx[1], sx[2])
    y  = RPMG::RESCALE(V, 0, 1, sx[1], sx[2])

#####  lines(x,y)
    points(x, y, pch=pch, col=col, cex=.5)

    x  = RPMG::RESCALE(N, 2, 3, sx[1], sx[2])
    y  = RPMG::RESCALE(V, 0, 1, sx[1], sx[2])
    points(x, y, pch=pch, col=col, cex=.5)


  }

