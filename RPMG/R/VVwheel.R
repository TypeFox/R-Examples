`VVwheel` <-
function( BIGMESH=NULL, v=1)
  {
    if(missing(v)) {  v=1 }
    if(missing(BIGMESH))
      {
        M = meshgrid( seq(from=-0.7071068, to=0.7071068, length=100) ,  seq(from=-0.7071068, to=0.7071068, length=100)   )

    ARE = sqrt((M$x)^2+(M$y)^2)
    
    ARE[ARE>1] = 1

    ANG = atan2(M$y,M$x)
    
    pANG = (ANG+pi)/(2*pi)
    
    dx = diff(M$x[1,])[1]
    dy = diff(M$y[,1])[1]
    
    
    RY =  range(M$y)
    RX =  range(M$x)

        BIGMESH=list(M=M, ARE=ARE, pANG = pANG, dx=dx, dy=dy, RY=RY, RX=RX)

      }

 
    plot(BIGMESH$RX, range(c(BIGMESH$RX  , BIGMESH$RX+diff(BIGMESH$RY))), type='n', xlab="", ylab="", axes=FALSE)
    mycols = hsv( BIGMESH$pANG,   BIGMESH$ARE, v )
    rect(BIGMESH$M$x, BIGMESH$M$y, BIGMESH$M$x+BIGMESH$dx, BIGMESH$M$y+BIGMESH$dy, col=mycols, border="transparent")
    mycols = hsv( BIGMESH$pANG,  v, BIGMESH$ARE )
    rect(BIGMESH$M$x, BIGMESH$M$y+diff(BIGMESH$RY), BIGMESH$M$x+BIGMESH$dx, BIGMESH$M$y+BIGMESH$dy+diff(BIGMESH$RY), col=mycols, border="transparent")


    invisible(BIGMESH)
    
  }

