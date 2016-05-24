rectPERIM<-function(x, y=1, pct=0)
  {
    if(missing(pct)) { pct=0}
    if(missing(y)) { y = x$y }
    if(is.list(x))
      {
        y = x$y
        x = x$x
      }

    rx = range(x, na.rm=TRUE)
    ry = range(y, na.rm=TRUE)

 if(pct>0 & pct>1)
      {
        pct = pct/100
      }
    if(pct>0)
      {
        rx = range( c(rx[1]-pct*(rx[2]-rx[1]), rx[2]+pct*(rx[2]-rx[1])))
        ry = range( c(ry[1]-pct*(ry[2]-ry[1]), ry[2]+pct*(ry[2]-ry[1])))

      }


    
    KX = c(rx[1], rx[2], rx[2], rx[1])
    KY = c(ry[1], ry[1], ry[2] , ry[2]  )

    return(list(x=KX, y=KY))
  }

