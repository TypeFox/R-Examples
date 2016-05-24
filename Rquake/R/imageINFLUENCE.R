imageINFLUENCE<-function(B, sta, proj)
  {
    NEX = 50
    col=RSEIS::tomo.colors(50)
     Mz = match(B$names, sta$name)
     
     gxy = GEOmap::GLOB.XY(  sta$lat[Mz]   ,sta$lon[Mz] , proj)
     
     ex = seq(from=min(gxy$x), to=max(gxy$x), length=NEX)
    why = seq(from=min(gxy$y), to=max(gxy$y), length=NEX)

     mval = as.vector(B$stats[3, ])

    DF = cbind(x=gxy$x , y=gxy$y ,  z=mval)
 zed  = MBA::mba.surf(DF, NEX, NEX, n = 8, m = 8, h = 8, extend=TRUE)$xyz.est
     ex = zed[[1]]
    why = zed[[2]]

    
     ## zed  = interp(x=gxy$x , y=gxy$y, z=mval, ex, why, duplicate="mean" )

     image(zed, col=col , add=TRUE)
     
     points(gxy$x, gxy$y)

     text(gxy$x, gxy$y, labels=B$names, pos=3)

   }
