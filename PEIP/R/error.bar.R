error.bar <-
function(x, y, lo, hi, pch=1, col=1, barw=0.1, add=FALSE, ...)
  {
    if(missing(add)) { add=FALSE; }
    if(missing(pch)) { pch=1 }
    if(missing(col)) { col=1 }
    if(missing(barw)) {  barw=0.1 }
   

    
    if(add==FALSE)
      {
    plot(x,y, xlim=range(x[!is.na(x)]), ylim=range(c(lo[!is.na(lo)], hi[!is.na(hi)])), type='n', xlab="", ylab="", ...)
  }
    fi = par("fin")
    u = par("usr")
    w = barw*((u[2]-u[1])/fi[1])
    points(x,y, pch=pch, col=col)
    segments(x,lo, x, hi, col=col)
    segments(x-w,lo,x+w, lo, col=col)
    segments(x-w,hi,x+w, hi, col=col)
  }
