`PLOT.TOMOXSEC` <-
function(XZSEC , depth=c(-25, 0) , COL=NULL, LIM=NULL, STA=NULL, ADD=FALSE)
  {
    ##require(GEOmap)
    
    if(missing(COL)) { COL=tomo.colors(100) }
    if(missing(LIM)) { LIM=NULL }
    if(missing(STA)) { STA=NULL }
    if(missing(depth)) { depth=range(XZSEC$y) }

    if(is.null(depth)) { depth=range(XZSEC$y) }
    if(missing(ADD)) {ADD=FALSE }

    
    if(ADD==FALSE)
      {
    plot(range(XZSEC$x) ,range(XZSEC$y) , ylim=depth, xlab="Km", ylab="Depth, Km", type='n', asp=TRUE )
  }
    
    ## image(XZSEC$x, XZSEC$y, XZSEC$z, col=COL, asp=TRUE, add=TRUE )

## abline(v=XZSEC$x, h=XZSEC$y)

    MXY1 = meshgrid(XZSEC$x[1:(length(XZSEC$x)-1)], XZSEC$y[1:(length(XZSEC$y)-1)])
    MXY2 = meshgrid(XZSEC$x[2:(length(XZSEC$x))], XZSEC$y[2:(length(XZSEC$y))])
##  rect( MXY1$x,  MXY1$y,  MXY2$x,  MXY2$y )

    if(is.null(LIM))
      {
        m1 = min(XZSEC$z, na.rm=TRUE)
        m2 = max(XZSEC$z, na.rm=TRUE)
        brs = seq(from=m1, to=m2, length=length(COL)   )
        
        
      }
    else
      {

        if(length(LIM)==2)
          { 
            m1 = LIM[1]
            m2 = LIM[2]
            brs = seq(from=m1, to=m2, length=length(COL)   )
            
          }
        else
          {

            brs = LIM
          }


        
      }
  

    
fs = findInterval( XZSEC$z, brs)
 mfs =   t(matrix(fs, ncol=ncol(XZSEC$z), nrow=nrow(XZSEC$z)))

    
##############    here the colors are flipped!

 ### rect( MXY1$x,  MXY1$y,  MXY2$x,  MXY2$y, border=NA, col= COL[100-mfs+1]  )       
   
    
   rect( MXY1$x,  MXY1$y,  MXY2$x,  MXY2$y, border=NA, col= COL[mfs]  )       
   
  }

