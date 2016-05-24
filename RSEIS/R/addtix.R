`addtix` <-
function(side=3, pos=0,   tck=0.005, at=c(0,1), labels=FALSE, col=2, ...)
  {
 ##X##   add tick marks to plot   
##X## ###     addtix(side=3, pos=y3+dy,   tck=0.005, at=c(0,1), labels=FALSE, col=2 )
    if(missing(side)) {  side = 3 }
    if(missing(pos)) {  pos = 0 }
    if(missing(tck)) {  tck = 0 }
    if(missing(at)) {  at = 0 }
    if(missing(labels)) {  labels = FALSE }
    if(missing(col)) {  col = "black" }
    

    n = length(at)
    u = par('usr')
    lines( c(at[1], at[n]), c(pos,pos) , col=col, ...)
    ###  x0, y0, x1, y1
    
    segments(at,pos, at, pos-tck , col=col)

  }

