`pnet` <-
function(MN, add=FALSE, col = gray(.7), border='black', lwd=1 )
  {
    ###########
    if(missing(add))  { add=FALSE }
    if(missing(col))  { col = gray(.7) }  
    if(missing(lwd))  { lwd=1 }
    if(missing(border))  { border='black' }
######   pnet(MN)
######   pnet(MN, col = 'brown' , border='black', lwd=.5)

    if(add==FALSE)
      {
        plot(c(-1,1),c(-1,1), type='n', xlab='', ylab='', asp=1, axes=FALSE) 
      }
    pcirc(col)

    lines(MN$x1, MN$y1, col=col)
    lines(MN$x2, MN$y2, col=col)
    segments(c(-.02, 0), c(0, -0.02), c(0.02, 0), c(0, 0.02), col=1)
  }

