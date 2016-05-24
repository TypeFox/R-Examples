PostREQquake<-function(XQ, proj)
  {

    ZEDS = matrix(ncol=3, nrow=length(XQ))
    
    for(i in 1:length(XQ))
      {
        
        ZEDS[i, ] = c(XQ[[i]]$EQ$lat,XQ[[i]]$EQ$lon, XQ[[i]]$EQ$z ) 
      }

    
    ZEXY = GEOmap::GLOB.XY(ZEDS[,1], ZEDS[,2], proj)
    points(ZEXY, pch=8, col='red')

    for(i in 1:length(XQ))
      {
        KOV =  XQ[[i]]$ERR$cov[2:4, 2:4]
        ndf = XQ[[i]]$ERR$ndf
        
        eqlipse(ZEXY$x[i], ZEXY$y[i] , KOV,   wcols = c(1,2) , dof=ndf, border="blue"  )
      }
    
    

  }


