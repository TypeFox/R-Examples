`PlotTernfoc` <-
function(h,v, x=0, y=0, siz=1, fcols='black', LABS=FALSE, add=FALSE)
  {
    if(missing(LABS)) { LABS=FALSE }
    if(missing(fcols)) { fcols=rep("black", length(h)) }
    if(missing(x)) { x=0 }
    if(missing(y)) { y=0 }
    if(missing(siz)) { siz=1 }
    if(missing(add)) { add=FALSE }
    
    
    V1 = ternfoc.point(90.0, 0.0 , 0.0)
    V2 = ternfoc.point(0.0, 90.0 , 0.0)
    V3 = ternfoc.point(0.0, 0.0 , 90.0)
    
    if(add==FALSE) plot( x+siz*c(V1$h, V2$h, V3$h),  y+siz*c(V1$v, V2$v, V3$v), type='n', ann=FALSE, axes=FALSE, asp=1  )
    polygon( x+siz*c(V1$h, V2$h, V3$h),  y+siz*c(V1$v, V2$v, V3$v), col='white' , border="black")
    
    points(x+siz*h,y+siz*v, pch=4, cex=.8, col=fcols, lwd=2)
    if(LABS==TRUE)
      {
        text( x+siz*c(V1$h, V2$h, V3$h),  y+siz*c(V1$v, V2$v, V3$v),
         labels=c("StrikeSlip", "Normal", "Thrust"), pos=c(3,1,1), xpd=TRUE)
      }
    
    
  }

