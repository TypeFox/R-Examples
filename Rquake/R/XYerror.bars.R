`XYerror.bars` <-
function(x, y, xlo=0, xhi=0, ylo=0, yhi=0, pch=1, col=1, barw=0.1, add=FALSE, ...)
  {
    if(missing(add)) { add=FALSE; }
    if(missing(pch)) { pch=1 }
    if(missing(col)) { col=1 }
    if(missing(barw)) {  barw=0.1 }



    if(add==FALSE)
      {
    plot(x,y, xlim=range(c(ylo[!is.na(ylo)], yhi[!is.na(yhi)])),
         ylim=range(c(xlo[!is.na(xlo)], xhi[!is.na(xhi)])),
         type='n', xlab="", ylab="", ...)
  }
    fi = par("fin")
    u = par("usr")
    xw = barw*((u[2]-u[1])/fi[1])
    yw = barw*((u[4]-u[3])/fi[2])
    
    points(x,y, pch=pch, col=col)
    
    segments(x,xlo, x, xhi, col=col)
    segments(x-xw,xlo,x+xw, xlo, col=col)
    segments(x-xw,xhi,x+xw, xhi, col=col)

    segments(ylo,y, yhi, y, col=col)
    segments(ylo, y-yw, ylo, y+yw, col=col)
    
    segments(yhi, y-yw ,  yhi, y+yw, col=col)
    
  }

